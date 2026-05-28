"""Tests for pounce.shortread.matrix."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pounce.shortread.matrix import build_count_matrix, load_count_tsv


def _write_tsv(path: Path, counts: dict[str, int]) -> None:
    with open(path, "w") as f:
        for bc, n in sorted(counts.items()):
            f.write(f"{bc}\t{n}\n")


def test_load_count_tsv(tmp_path):
    _write_tsv(tmp_path / "s.tsv", {"AAAA": 5, "CCCC": 10})
    s = load_count_tsv(tmp_path / "s.tsv")
    assert s["AAAA"] == 5
    assert s["CCCC"] == 10


def test_build_matrix_zeros_absent_barcodes(tmp_path, data_dir):
    """Barcodes absent in a sample should get count 0, not NaN."""
    plate = pd.read_csv(data_dir / "test_plate_info.csv")
    plate["count_file"] = plate["sample_id"].apply(
        lambda sid: str(data_dir / f"{sid}_barcode_counts.tsv")
    )
    mat = build_count_matrix(plate, ns_min_count=1)

    # Every cell should be non-negative integer
    assert (mat >= 0).all().all()
    assert mat.dtypes.apply(lambda d: d in [int, "int64", "int32"]).all() or \
           (mat.fillna(0) == mat.astype(int)).all().all()


def test_build_matrix_ns_filter(tmp_path):
    """Barcodes not meeting ns_min_count in any NS sample are excluded."""
    ns_path = tmp_path / "NS.tsv"
    s_path = tmp_path / "S.tsv"
    _write_tsv(ns_path, {"AAAA": 1, "CCCC": 10})
    _write_tsv(s_path, {"AAAA": 50, "CCCC": 5})

    plate = pd.DataFrame([
        {"sample_id": "NS", "fastq_r2": "", "timepoint": "T1",
         "replicate": "A", "condition": "NS", "count_file": str(ns_path)},
        {"sample_id": "S", "fastq_r2": "", "timepoint": "T1",
         "replicate": "A", "condition": "S", "count_file": str(s_path)},
    ])

    mat = build_count_matrix(plate, ns_min_count=5)
    # "AAAA" has only 1 NS count < 5, so should be excluded
    assert "AAAA" not in mat.index
    assert "CCCC" in mat.index


def test_build_matrix_uses_all_ns_samples(tmp_path):
    """A barcode should be included if it meets ns_min_count in ANY NS sample."""
    ns1 = tmp_path / "NS1.tsv"
    ns2 = tmp_path / "NS2.tsv"
    _write_tsv(ns1, {"AAAA": 1})
    _write_tsv(ns2, {"AAAA": 10})  # meets threshold in NS2

    plate = pd.DataFrame([
        {"sample_id": "NS1", "fastq_r2": "", "timepoint": "T1",
         "replicate": "A", "condition": "NS", "count_file": str(ns1)},
        {"sample_id": "NS2", "fastq_r2": "", "timepoint": "T1",
         "replicate": "B", "condition": "NS", "count_file": str(ns2)},
    ])

    mat = build_count_matrix(plate, ns_min_count=5)
    assert "AAAA" in mat.index
