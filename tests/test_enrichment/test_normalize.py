"""Tests for pounce.enrichment.normalize."""

from __future__ import annotations

import pandas as pd
import pytest

from pounce.enrichment.normalize import normalize_to_wt_stop


def _make_variant_df(variants: dict[str, float], count: int = 15) -> pd.DataFrame:
    rows = {v: {"mean_enrich": e, "count_A": count} for v, e in variants.items()}
    return pd.DataFrame.from_dict(rows, orient="index")


def test_wt_maps_to_one():
    df = _make_variant_df({"WT": 2.0, "K5*": -1.0, "E13N": 1.5})
    normalized = normalize_to_wt_stop(df, min_barcodes=1)
    assert abs(normalized["WT"] - 1.0) < 1e-9


def test_stop_maps_to_zero():
    df = _make_variant_df({"WT": 2.0, "K5*": -1.0, "E13N": 1.5})
    normalized = normalize_to_wt_stop(df, min_barcodes=1)
    assert abs(normalized["K5*"] - 0.0) < 1e-9


def test_intermediate_variant_scaled():
    df = _make_variant_df({"WT": 4.0, "K5*": 0.0, "E13N": 2.0})
    normalized = normalize_to_wt_stop(df, min_barcodes=1)
    # E13N is midway between WT and stop → should be ≈ 0.5
    assert abs(normalized["E13N"] - 0.5) < 1e-9


def test_multiple_stops_uses_median():
    df = _make_variant_df({
        "WT": 2.0,
        "K5*": -1.0,
        "A10*": -0.8,  # slightly different stop enrichment
        "E13N": 1.0,
    })
    # median of stops: (-1.0 + -0.8) / 2 = -0.9
    normalized = normalize_to_wt_stop(df, min_barcodes=1)
    assert abs(normalized["WT"] - 1.0) < 1e-9


def test_no_wt_raises():
    df = _make_variant_df({"K5*": -1.0, "E13N": 1.5})
    with pytest.raises(ValueError, match="WT-synonymous"):
        normalize_to_wt_stop(df, min_barcodes=1)


def test_no_stops_raises():
    df = _make_variant_df({"WT": 2.0, "E13N": 1.5})
    with pytest.raises(ValueError, match="stop"):
        normalize_to_wt_stop(df, min_barcodes=1)


def test_min_barcodes_filters_anchors():
    """Anchors with fewer than min_barcodes should not be used."""
    df = _make_variant_df({"WT": 2.0, "K5*": -1.0, "E13N": 1.5}, count=5)
    with pytest.raises(ValueError):
        normalize_to_wt_stop(df, min_barcodes=10)
