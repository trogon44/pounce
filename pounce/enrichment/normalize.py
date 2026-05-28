from __future__ import annotations

import numpy as np
import pandas as pd


def _is_stop_variant(variant: str) -> bool:
    return "*" in variant and ";" not in variant  # single-mutation stop


def _is_wt_synonymous(variant: str) -> bool:
    return variant == "WT"


def normalize_to_wt_stop(
    variant_enrichments: pd.DataFrame,
    min_barcodes: int = 10,
    enrich_col: str = "mean_enrich",
) -> pd.Series:
    """Normalize enrichments so that WT synonymous = 1.0 and single stops = 0.0.

    Anchors are computed from the observed distributions of:
    - WT-synonymous variants (variant == 'WT') with >= min_barcodes barcodes
    - Single premature-stop variants (one mutation, ends in '*') with >= min_barcodes barcodes

    Returns a Series aligned to variant_enrichments.index.
    """
    count_cols = [c for c in variant_enrichments.columns if c.startswith("count_")]

    def _has_min_barcodes(row: pd.Series) -> bool:
        return any(row[c] >= min_barcodes for c in count_cols if c in row.index)

    sufficient = variant_enrichments.apply(_has_min_barcodes, axis=1)

    wt_mask = (
        variant_enrichments.index.map(_is_wt_synonymous)
        & sufficient
    )
    stop_mask = (
        variant_enrichments.index.map(_is_stop_variant)
        & sufficient
    )

    wt_vals = variant_enrichments.loc[wt_mask, enrich_col].dropna()
    stop_vals = variant_enrichments.loc[stop_mask, enrich_col].dropna()

    if wt_vals.empty:
        raise ValueError(
            "No WT-synonymous variants with sufficient barcode coverage for normalization. "
            "Check that barcode_variants.csv was generated from the same library."
        )
    if stop_vals.empty:
        raise ValueError(
            "No single-stop variants with sufficient barcode coverage for normalization."
        )

    wt_mean = wt_vals.median()
    stop_mean = stop_vals.median()
    scale = wt_mean - stop_mean

    if abs(scale) < 1e-9:
        raise ValueError("WT and stop enrichments are indistinguishable; normalization failed.")

    normalized = (variant_enrichments[enrich_col] - stop_mean) / scale
    normalized.name = "enrich_normalized"
    return normalized
