from __future__ import annotations

import pandas as pd


def load_barcode_variants(path: str) -> pd.DataFrame:
    """Load barcode_variants.csv produced by ``call_all_mutations``.

    Returns DataFrame with barcode as index and columns:
    nt_mutations, aa_mutations, read_count.
    """
    df = pd.read_csv(path, index_col="barcode")
    df["nt_mutations"] = df["nt_mutations"].fillna("")
    df["aa_mutations"] = df["aa_mutations"].fillna("")
    return df


def aa_string(aa_mutations: str) -> str:
    """Normalize an aa_mutations string to a canonical variant key.

    Sorts individual mutations so 'E13N;D5A' and 'D5A;E13N' map to the same key.
    Empty string maps to 'WT'.
    """
    if not aa_mutations or aa_mutations == "":
        return "WT"
    parts = [p.strip() for p in aa_mutations.split(";") if p.strip()]
    return ";".join(sorted(parts))


def add_variant_column(barcode_variants: pd.DataFrame) -> pd.DataFrame:
    """Add a 'variant' column with the canonical AA variant string."""
    df = barcode_variants.copy()
    df["variant"] = df["aa_mutations"].apply(aa_string)
    return df


def aggregate_enrichments(
    enrichments: pd.DataFrame,
    barcode_variants: pd.DataFrame,
    min_barcodes: int = 3,
) -> pd.DataFrame:
    """Group per-barcode enrichments by AA variant and compute median per replicate.

    Parameters
    ----------
    enrichments:
        Per-barcode enrichment DataFrame (output of compute_replicate_enrichments).
        Index = barcode, columns = enrich_A, enrich_B, ...
    barcode_variants:
        barcode_variants.csv loaded via load_barcode_variants, with 'variant' column
        (add via add_variant_column if not present).
    min_barcodes:
        Variants represented by fewer barcodes than this are dropped.

    Returns
    -------
    DataFrame indexed by variant string, with columns:
    median_A, median_B, ..., count_A, count_B, ..., mean_enrich.
    """
    if "variant" not in barcode_variants.columns:
        barcode_variants = add_variant_column(barcode_variants)

    merged = enrichments.join(barcode_variants[["variant"]], how="inner")

    enrich_cols = [c for c in enrichments.columns if c.startswith("enrich_")]
    rep_labels = [c.replace("enrich_", "") for c in enrich_cols]

    records: list[dict] = []
    for variant, group in merged.groupby("variant"):
        row: dict = {"variant": variant}
        medians = []
        for col, label in zip(enrich_cols, rep_labels):
            vals = group[col].dropna()
            row[f"count_{label}"] = len(vals)
            row[f"median_{label}"] = vals.median() if len(vals) >= min_barcodes else float("nan")
            if len(vals) >= min_barcodes:
                medians.append(vals.median())
        row["mean_enrich"] = sum(medians) / len(medians) if medians else float("nan")
        records.append(row)

    result = pd.DataFrame(records).set_index("variant")
    return result[result[[c for c in result.columns if c.startswith("count_")]].max(axis=1) >= min_barcodes]
