"""Tests for pounce.enrichment.compute."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from pounce.enrichment.compute import compute_enrichment, compute_replicate_enrichments


def _make_counts(s_vals: list[int], ns_vals: list[int]) -> pd.DataFrame:
    barcodes = [f"BC{i:02d}" for i in range(len(s_vals))]
    return pd.DataFrame({"S": s_vals, "NS": ns_vals}, index=barcodes)


def test_equal_frequencies_give_zero():
    """Barcode at same frequency in S and NS should give enrichment ≈ 0."""
    # All barcodes identical count in both pools
    counts = _make_counts([10] * 5, [10] * 5)
    enrich = compute_enrichment(counts, "S", "NS", pseudocount=1.0, ns_min=0)
    assert all(abs(v) < 0.01 for v in enrich)


def test_higher_s_gives_positive_enrichment():
    counts = _make_counts([100, 1], [10, 10])
    enrich = compute_enrichment(counts, "S", "NS", pseudocount=1.0, ns_min=0)
    assert enrich.iloc[0] > 0
    assert enrich.iloc[1] < 0


def test_ns_min_filter():
    counts = _make_counts([50, 50], [3, 20])
    enrich = compute_enrichment(counts, "S", "NS", pseudocount=1.0, ns_min=5)
    # BC00 has NS=3, below ns_min=5 → excluded
    assert "BC00" not in enrich.index
    assert "BC01" in enrich.index


def test_pseudocount_affects_low_count_barcodes():
    """With pc=1, a barcode with 0 S counts should still produce a finite value."""
    counts = _make_counts([0, 10], [10, 10])
    enrich = compute_enrichment(counts, "S", "NS", pseudocount=1.0, ns_min=0)
    assert math.isfinite(enrich.iloc[0])
    assert enrich.iloc[0] < enrich.iloc[1]


def test_zero_s_total_raises():
    counts = _make_counts([0, 0], [10, 10])
    with pytest.raises(ValueError, match="Zero total reads"):
        compute_enrichment(counts, "S", "NS", pseudocount=0.0, ns_min=0)


def test_replicate_enrichments(data_dir):
    plate = pd.read_csv(data_dir / "test_plate_info.csv")
    plate["count_file"] = plate["sample_id"].apply(
        lambda sid: str(data_dir / f"{sid}_barcode_counts.tsv")
    )
    from pounce.shortread.matrix import build_count_matrix
    counts = build_count_matrix(plate, ns_min_count=1)
    result = compute_replicate_enrichments(counts, plate, timepoint="T1")

    assert "enrich_A" in result.columns
    assert "enrich_B" in result.columns
    assert len(result) > 0
    # All values finite
    assert result.apply(lambda col: col.apply(math.isfinite)).all().all()


def _plate(rows):
    return pd.DataFrame(rows, columns=["sample_id", "timepoint", "replicate", "condition"])


def test_multiple_selective_samples_per_replicate_raises():
    """Silently taking the first S row once paired a replicate with the wrong sample."""
    counts = pd.DataFrame({"S1": [5, 5], "S2": [9, 9], "NS": [10, 10]}, index=["a", "b"])
    plate = _plate([("S1", "T2", "A", "S"), ("S2", "T2", "A", "S"), ("NS", "T2", "A", "NS")])
    with pytest.raises(ValueError, match="selective samples"):
        compute_replicate_enrichments(counts, plate, timepoint="T2", ns_min=0)


def test_multiple_ns_samples_per_replicate_raises():
    counts = pd.DataFrame({"S": [5, 5], "NS1": [10, 10], "NS2": [9, 9]}, index=["a", "b"])
    plate = _plate([("S", "T2", "A", "S"), ("NS1", "T2", "A", "NS"), ("NS2", "T2", "A", "NS")])
    with pytest.raises(ValueError, match="non-selective samples"):
        compute_replicate_enrichments(counts, plate, timepoint="T2", ns_min=0)


def test_replicates_are_filtered_independently_on_their_own_ns_column():
    counts = pd.DataFrame(
        {"SA": [20, 20], "NSA": [30, 2], "SB": [20, 20], "NSB": [2, 30]}, index=["x", "y"]
    )
    plate = _plate([
        ("SA", "T2", "A", "S"), ("NSA", "T2", "A", "NS"),
        ("SB", "T2", "B", "S"), ("NSB", "T2", "B", "NS"),
    ])
    out = compute_replicate_enrichments(counts, plate, timepoint="T2", ns_min=10)
    assert out["enrich_A"].notna().to_dict() == {"x": True, "y": False}
    assert out["enrich_B"].notna().to_dict() == {"x": False, "y": True}


STACKED_PLATE = pd.DataFrame(
    [
        ("SP01", "T0", "A", "NS", "00h"), ("SP03", "T0", "A", "NS", "24h"),
        ("SP11", "T0", "A", "S", "72h"), ("SP12", "T0", "A", "NS", "72h"),
        ("SP02", "T0", "B", "NS", "00h"),
        ("SP13", "T0", "B", "S", "72h"), ("SP14", "T0", "B", "NS", "72h"),
        ("SP07", "T2", "A", "S", "48h"), ("SP08", "T2", "A", "NS", "48h"),
    ],
    columns=["sample_id", "timepoint", "replicate", "condition", "harvest"],
)


def _stacked_counts():
    rng = np.random.default_rng(0)
    cols = STACKED_PLATE["sample_id"].tolist()
    return pd.DataFrame(rng.integers(20, 200, size=(30, len(cols))), columns=cols,
                        index=[f"b{i}" for i in range(30)])


def test_s_is_paired_with_ns_of_the_same_harvest():
    counts = _stacked_counts()
    out = compute_replicate_enrichments(counts, STACKED_PLATE, "T0", ns_min=0)
    assert list(out.columns) == ["enrich_A", "enrich_B"]
    expected_a = compute_enrichment(counts, "SP11", "SP12", ns_min=0)   # 72h with 72h
    wrong_a = compute_enrichment(counts, "SP11", "SP01", ns_min=0)      # 72h with 00h
    pd.testing.assert_series_equal(out["enrich_A"], expected_a, check_names=False)
    assert not np.allclose(out["enrich_A"], wrong_a)


def test_early_ns_only_harvests_are_ignored():
    out = compute_replicate_enrichments(_stacked_counts(), STACKED_PLATE, "T0", ns_min=0)
    assert out.shape[1] == 2  # 00h/24h NS-only groups produce no column


def test_several_paired_harvests_require_a_choice():
    plate = pd.concat([STACKED_PLATE, pd.DataFrame(
        [("SP20", "T0", "A", "S", "48h"), ("SP21", "T0", "A", "NS", "48h")],
        columns=STACKED_PLATE.columns)], ignore_index=True)
    counts = pd.DataFrame(1, index=["a", "b"], columns=plate["sample_id"])
    with pytest.raises(ValueError, match="several harvests"):
        compute_replicate_enrichments(counts, plate, "T0", ns_min=0)
    out = compute_replicate_enrichments(counts + 5, plate, "T0", harvest="48h", ns_min=0)
    assert list(out.columns) == ["enrich_A"]


def test_no_complete_pair_raises():
    with pytest.raises(ValueError, match="No complete S/NS pair"):
        compute_replicate_enrichments(_stacked_counts(), STACKED_PLATE, "T0", harvest="00h", ns_min=0)


def test_harvest_without_column_raises():
    counts = pd.DataFrame(1, index=["a"], columns=["S", "NS"])
    plate = _plate([("S", "T2", "A", "S"), ("NS", "T2", "A", "NS")])
    with pytest.raises(ValueError, match="no 'harvest' column"):
        compute_replicate_enrichments(counts, plate, "T2", harvest="24h", ns_min=0)
