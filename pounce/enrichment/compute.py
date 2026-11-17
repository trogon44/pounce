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

    where ``pc`` is a raw read count (default 1) added to both numerator and
    denominator, and S_total / NS_total are the column sums of ``counts``, so
    pass the full count matrix, not a subset.  The library-size term makes the
    enrichment 0 for a barcode present at the same frequency in both pools.

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
    harvest: str | None = None,
    pseudocount: float = 1.0,
    ns_min: int = 10,
) -> pd.DataFrame:
    """Compute enrichment for each replicate at a given timepoint.

    ``plate_info`` needs columns sample_id, timepoint, replicate and condition
    (``S`` or ``NS``), and optionally ``harvest``. Each selective sample is
    paired with the one non-selective sample that shares its timepoint,
    replicate and, if the column exists, harvest. Groups that have only an NS
    sample (for example an early-harvest control) are ignored. More than one
    S or NS sample in a group is an error rather than a guess, as is having
    several harvests with complete S/NS pairs and no ``harvest`` to choose one.

    Returns a DataFrame with one column per replicate (enrich_A, enrich_B, ...).
    """
    tp = plate_info[plate_info["timepoint"] == timepoint]
    has_harvest = "harvest" in plate_info.columns
    if harvest is not None:
        if not has_harvest:
            raise ValueError("harvest was given but plate_info has no 'harvest' column")
        tp = tp[tp["harvest"].astype(str) == str(harvest)]

    group_cols = ["replicate", "harvest"] if has_harvest else ["replicate"]
    pairs: list[tuple[str, str, str, str]] = []  # (replicate, harvest, s_col, ns_col)
    for key, rows in tp.groupby(group_cols, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        rep, hv = key[0], (key[1] if has_harvest else "")
        s_rows = rows[rows["condition"] == "S"]
        ns_rows = rows[rows["condition"] == "NS"]
        if s_rows.empty or ns_rows.empty:
            continue
        where = f"replicate {rep}" + (f", harvest {hv}" if has_harvest else "") + f" at {timepoint}"
        for kind, group in (("selective", s_rows), ("non-selective", ns_rows)):
            if len(group) > 1:
                raise ValueError(
                    f"{where} has {len(group)} {kind} samples ({', '.join(group['sample_id'])}); "
                    "expected exactly one. List only the S/NS pair you want to compare in plate_info."
                )
        pairs.append((str(rep), str(hv), s_rows.iloc[0]["sample_id"], ns_rows.iloc[0]["sample_id"]))

    if not pairs:
        raise ValueError(f"No complete S/NS pair found at timepoint {timepoint}")
    harvests = sorted({hv for _, hv, _, _ in pairs})
    if len(harvests) > 1:
        raise ValueError(
            f"S/NS pairs exist at several harvests ({', '.join(harvests)}) at {timepoint}; "
            "choose one with harvest=."
        )

    results = {
        f"enrich_{rep}": compute_enrichment(counts, s_col, ns_col, pseudocount, ns_min)
        for rep, _, s_col, ns_col in pairs
    }
    return pd.DataFrame(results)
