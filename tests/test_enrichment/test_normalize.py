"""Tests for pounce.enrichment.normalize."""

from __future__ import annotations

import pandas as pd
import pytest

from pounce.enrichment.normalize import normalize_to_anchors


def _table(rows: list[tuple[str, str, float, int, int]]) -> pd.DataFrame:
    """rows: (variant, variant_type, mean_enrich, count_A, count_B)."""
    return pd.DataFrame(
        rows, columns=["variant", "variant_type", "mean_enrich", "count_A", "count_B"]
    ).set_index("variant")


BASE = [
    ("WTsyn_A1G", "wt_syn", 4.0, 20, 20),
    ("WTsyn_C2T", "wt_syn", 6.0, 20, 20),
    ("K5*", "stop", -2.0, 20, 20),
    ("A10*", "stop", 0.0, 20, 20),
    ("E13N", "missense", 3.0, 20, 20),
]


def test_anchor_means_map_to_one_and_zero():
    out = normalize_to_anchors(_table(BASE), min_barcodes=10)
    # wt level 5.0, dead level -1.0 -> scale 6
    assert out["E13N"] == pytest.approx((3.0 + 1.0) / 6.0)
    assert out["WTsyn_A1G"] == pytest.approx((4.0 + 1.0) / 6.0)
    assert out.name == "enrich_norm"


def test_anchor_stat_median_differs_from_mean():
    rows = BASE + [("WTsyn_G3A", "wt_syn", 20.0, 20, 20)]
    mean_out = normalize_to_anchors(_table(rows), anchor_stat="mean")
    median_out = normalize_to_anchors(_table(rows), anchor_stat="median")
    # wt level: mean(4,6,20)=10; median=6
    assert mean_out["E13N"] == pytest.approx((3.0 + 1.0) / 11.0)
    assert median_out["E13N"] == pytest.approx((3.0 + 1.0) / 7.0)


def test_anchor_needs_min_barcodes_in_every_replicate():
    rows = [r for r in BASE if r[0] != "WTsyn_C2T"] + [("WTsyn_C2T", "wt_syn", 100.0, 50, 3)]
    out = normalize_to_anchors(_table(rows), min_barcodes=10)
    # WTsyn_C2T has only 3 barcodes in replicate B so it is not an anchor: wt level = 4.0
    assert out["E13N"] == pytest.approx((3.0 + 1.0) / 5.0)


def test_non_anchor_types_ignored():
    rows = BASE + [("F1-", "frameshift", -9.0, 50, 50), ("WT", "wt", 99.0, 50, 50)]
    assert normalize_to_anchors(_table(rows))["E13N"] == pytest.approx((3.0 + 1.0) / 6.0)


def test_no_wt_syn_raises():
    rows = [r for r in BASE if r[1] != "wt_syn"]
    with pytest.raises(ValueError, match="WT-synonymous"):
        normalize_to_anchors(_table(rows))


def test_no_stops_raises():
    rows = [r for r in BASE if r[1] != "stop"]
    with pytest.raises(ValueError, match="stop"):
        normalize_to_anchors(_table(rows))


def test_identical_anchor_levels_raise():
    rows = [("WTsyn_A1G", "wt_syn", 1.0, 20, 20), ("K5*", "stop", 1.0, 20, 20)]
    with pytest.raises(ValueError, match="indistinguishable"):
        normalize_to_anchors(_table(rows))


def test_bad_anchor_stat_raises():
    with pytest.raises(ValueError, match="anchor_stat"):
        normalize_to_anchors(_table(BASE), anchor_stat="mode")
