from __future__ import annotations

import re

import numpy as np
import pandas as pd

WT_KEY = "WT"
WT_SYN_PREFIX = "WTsyn_"

VARIANT_TYPES = ("missense", "stop", "wt_syn", "frameshift", "wt")


def load_barcode_variants(path: str) -> pd.DataFrame:
    """Load barcode_variants.csv produced by ``call_all_mutations``.

    Returns DataFrame with barcode as index and columns:
    nt_mutations, aa_mutations, read_count.
    """
    df = pd.read_csv(path, index_col="barcode")
    df["nt_mutations"] = df["nt_mutations"].fillna("")
    df["aa_mutations"] = df["aa_mutations"].fillna("")
    return df


_DIGITS = re.compile(r"\d+")


def _position_key(mutation: str) -> int:
    return int(_DIGITS.search(mutation).group())


def _canonical(mutations: str) -> str:
    """Join a ';'-separated mutation string in residue/nucleotide position order."""
    if ";" not in mutations:
        return mutations
    return ";".join(sorted(mutations.split(";"), key=_position_key))


def _count(mutations: pd.Series) -> pd.Series:
    return np.where(mutations == "", 0, mutations.str.count(";") + 1)


def classify_barcodes(barcode_variants: pd.DataFrame) -> pd.DataFrame:
    """Assign every barcode a ``variant`` key and ``variant_type``.

    Types, in order of precedence:

    - ``frameshift``: an amino-acid mutation ending in ``-`` (a codon with a
      deleted or uncalled base). Keyed by the full nucleotide mutation string.
    - ``stop``: contains a ``*`` mutation. Keyed by the amino-acid string.
    - ``missense``: one or more amino-acid mutations. Keyed by the amino-acid string.
    - ``wt_syn``: nucleotide changes that are all synonymous. Keyed by
      ``WTsyn_`` + the nucleotide string, because distinct synonymous
      alleles are independent replicates of wild-type activity.
    - ``wt``: no mutations at all, keyed ``WT``.

    Mutation strings are put in position order. Adds columns ``nt_mutations``,
    ``aa_mutations`` (canonicalized), ``nt_mutcount``, ``aa_mutcount``
    (0 for ``wt_syn`` and ``frameshift`` by definition), ``variant_type``
    and ``variant``.
    """
    df = barcode_variants[["nt_mutations", "aa_mutations"]].fillna("").copy()
    df["nt_mutations"] = df["nt_mutations"].map(_canonical)
    df["aa_mutations"] = df["aa_mutations"].map(_canonical)

    nt_count = _count(df["nt_mutations"])
    aa_count = _count(df["aa_mutations"])
    is_frameshift = df["aa_mutations"].str.contains("-", regex=False)
    is_stop = df["aa_mutations"].str.contains("*", regex=False) & ~is_frameshift

    variant_type = np.select(
        [is_frameshift, is_stop, aa_count >= 1, nt_count >= 1],
        ["frameshift", "stop", "missense", "wt_syn"],
        default="wt",
    )
    variant = np.select(
        [variant_type == "frameshift", variant_type == "wt_syn", variant_type == "wt"],
        [df["nt_mutations"], WT_SYN_PREFIX + df["nt_mutations"], WT_KEY],
        default=df["aa_mutations"],
    )

    df["nt_mutcount"] = nt_count
    df["aa_mutcount"] = np.where(np.isin(variant_type, ["frameshift", "wt_syn", "wt"]), 0, aa_count)
    df["variant_type"] = variant_type
    df["variant"] = variant
    return df


def _replicate_labels(enrichments: pd.DataFrame) -> list[str]:
    cols = [c for c in enrichments.columns if c.startswith("enrich_")]
    if not cols:
        raise ValueError("enrichments has no 'enrich_*' columns")
    return [c.removeprefix("enrich_") for c in cols]


def aggregate_enrichments(
    enrichments: pd.DataFrame,
    barcode_variants: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate per-barcode enrichments to per-variant statistics.

    Each replicate is summarised independently: a barcode only contributes to
    the replicates in which it has an enrichment (i.e. passed that replicate's
    NS filter), so ``count_A`` and ``count_B`` can differ. No minimum barcode
    count is applied; filter on the ``count_*`` columns afterwards.

    Parameters
    ----------
    enrichments:
        Per-barcode enrichments (``compute_replicate_enrichments``). Index =
        barcode, columns ``enrich_A``, ``enrich_B``, ...
    barcode_variants:
        ``barcode_variants.csv`` loaded via ``load_barcode_variants``, or the
        output of ``classify_barcodes``.

    Returns
    -------
    DataFrame indexed by ``variant`` with ``variant_type``, ``aa_mutcount`` and,
    for each replicate, ``count_*``, ``median_*``, ``mean_*``, ``std_*``, plus
    ``mean_enrich``: the mean of the per-replicate medians.
    """
    if "variant" not in barcode_variants.columns:
        barcode_variants = classify_barcodes(barcode_variants)

    labels = _replicate_labels(enrichments)
    classes = barcode_variants[["variant", "variant_type", "aa_mutcount"]]
    merged = enrichments.join(classes, how="inner")

    pieces = [merged.drop_duplicates("variant").set_index("variant")[["variant_type", "aa_mutcount"]]]
    for label in labels:
        values = merged.dropna(subset=[f"enrich_{label}"]).groupby("variant")[f"enrich_{label}"]
        stats = values.agg(["count", "median", "mean", "std"])
        stats.columns = [f"{name}_{label}" for name in stats.columns]
        pieces.append(stats)

    table = pd.concat(pieces, axis=1)
    for label in labels:
        table[f"count_{label}"] = table[f"count_{label}"].fillna(0).astype(int)

    table["mean_enrich"] = table[[f"median_{label}" for label in labels]].mean(axis=1)
    ordered = ["variant_type", "aa_mutcount"]
    for label in labels:
        ordered += [f"count_{label}", f"median_{label}", f"mean_{label}", f"std_{label}"]
    table = table[ordered + ["mean_enrich"]]
    table.index.name = "variant"
    return table.sort_index()
