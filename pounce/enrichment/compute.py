from __future__ import annotations

import numpy as np
import pandas as pd


def compute_enrichment(
    counts: pd.DataFrame,
    s_col: str,
    ns_col: str,
    pseudocount: float = 1.0,
    ns_min: int = 10,
) -> pd.Series:
    """Compute per-barcode log2 enrichment with library-size normalization.

    Formula:
        E = log2((S + pc) / (NS + pc)) + log2(NS_total / S_total)

    where pc = pseudocount / library_size, making ``pseudocount`` equivalent
    to that many raw counts.  The library-size correction term ensures that
    enrichment = 0 for a barcode present at the same frequency in both pools.

    Parameters
    ----------
    counts:
        Barcode count matrix.  Index = barcode, columns include s_col and ns_col.
    s_col:
        Column name for the selective sample.
    ns_col:
        Column name for the non-selective sample.
    pseudocount:
        Raw-count pseudocount added symmetrically to numerator and denominator.
    ns_min:
        Minimum NS count required per barcode; barcodes below this are excluded.

    Returns
    -------
    Series of enrichment values indexed by barcode (only barcodes passing ns_min).
    """
    s_total = counts[s_col].sum()
    ns_total = counts[ns_col].sum()

    if s_total == 0 or ns_total == 0:
        raise ValueError(f"Zero total reads in {'S' if s_total == 0 else 'NS'} column")

    mask = counts[ns_col] >= ns_min
    sub = counts.loc[mask, [s_col, ns_col]]

    enrich = (
        np.log2((sub[s_col] + pseudocount) / (sub[ns_col] + pseudocount))
        + np.log2(ns_total / s_total)
    )
    enrich.name = f"enrich_{s_col}_vs_{ns_col}"
    return enrich


def compute_replicate_enrichments(
    counts: pd.DataFrame,
    plate_info: pd.DataFrame,
    timepoint: str,
    pseudocount: float = 1.0,
    ns_min: int = 10,
) -> pd.DataFrame:
    """Compute enrichment for each replicate at a given timepoint.

    plate_info must have columns: sample_id, timepoint, replicate, condition.
    Returns a DataFrame with one column per replicate (enrich_A, enrich_B, ...).
    """
    tp = plate_info[plate_info["timepoint"] == timepoint]
    replicates = tp["replicate"].unique()

    results: dict[str, pd.Series] = {}
    for rep in sorted(replicates):
        rep_rows = tp[tp["replicate"] == rep]
        s_rows = rep_rows[rep_rows["condition"] == "S"]
        ns_rows = rep_rows[rep_rows["condition"] == "NS"]
        if s_rows.empty or ns_rows.empty:
            continue
        s_col = s_rows.iloc[0]["sample_id"]
        ns_col = ns_rows.iloc[0]["sample_id"]
        enrich = compute_enrichment(counts, s_col, ns_col, pseudocount, ns_min)
        results[f"enrich_{rep}"] = enrich

    return pd.DataFrame(results)
