from __future__ import annotations

import pandas as pd


def normalize_to_anchors(
    variant_enrichments: pd.DataFrame,
    min_barcodes: int = 10,
    anchor_stat: str = "mean",
    enrich_col: str = "mean_enrich",
) -> pd.Series:
    """Rescale enrichments so the dead level is 0 and the wild-type level is 1.

    Anchors are the ``stop`` variants (dead) and the ``wt_syn`` variants
    (wild-type) that have at least ``min_barcodes`` barcodes in *every*
    replicate. Each anchor level is the ``anchor_stat`` ("mean" or "median") of
    ``enrich_col`` across its anchor variants, and every variant becomes
    ``(E - dead) / (wt - dead)``.

    Returns a Series named ``enrich_norm`` aligned to ``variant_enrichments``.
    """
    if anchor_stat not in ("mean", "median"):
        raise ValueError(f"anchor_stat must be 'mean' or 'median', not {anchor_stat!r}")

    count_cols = [c for c in variant_enrichments.columns if c.startswith("count_")]
    if not count_cols:
        raise ValueError("variant_enrichments has no 'count_*' columns")
    sufficient = (variant_enrichments[count_cols] >= min_barcodes).all(axis=1)

    levels = {}
    for variant_type, label in (("wt_syn", "WT-synonymous"), ("stop", "stop")):
        anchors = variant_enrichments.loc[
            (variant_enrichments["variant_type"] == variant_type) & sufficient, enrich_col
        ].dropna()
        if anchors.empty:
            raise ValueError(
                f"No {label} variants with >= {min_barcodes} barcodes in every replicate; "
                "cannot normalize."
            )
        levels[variant_type] = getattr(anchors, anchor_stat)()

    scale = levels["wt_syn"] - levels["stop"]
    if abs(scale) < 1e-9:
        raise ValueError("WT and stop enrichments are indistinguishable; normalization failed.")

    normalized = (variant_enrichments[enrich_col] - levels["stop"]) / scale
    normalized.name = "enrich_norm"
    return normalized
