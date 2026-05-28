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
