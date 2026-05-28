"""Tests for pounce.shortread.count (pure-Python path)."""

from __future__ import annotations

import gzip
from pathlib import Path

from pounce.shortread.count import count_barcodes_python

BARCODE_LEN = 30


def _write_fastq(path: Path, barcode_counts: dict[str, int]) -> None:
    adapter = "ATCG" * 5  # dummy flanks, not included in this test
    with gzip.open(path, "wt") as f:
        read_i = 0
        for bc, n in barcode_counts.items():
            seq = bc  # pure barcode reads (already trimmed)
            qual = "I" * len(seq)
            for _ in range(n):
                f.write(f"@read{read_i}\n{seq}\n+\n{qual}\n")
                read_i += 1


def test_count_matches_input(tmp_path):
    barcodes = {
        "A" * BARCODE_LEN: 10,
        "C" * BARCODE_LEN: 5,
        "G" * BARCODE_LEN: 1,
    }
    fastq = tmp_path / "sample_barcodes.fastq.gz"
    _write_fastq(fastq, barcodes)

    tsv = tmp_path / "counts.tsv"
    count_barcodes_python(fastq, tsv, BARCODE_LEN)

    result = {}
    for line in tsv.read_text().strip().splitlines():
        bc, cnt = line.split("\t")
        result[bc] = int(cnt)

    assert result == barcodes


def test_wrong_length_excluded(tmp_path):
    """Barcodes of wrong length must not appear in output."""
    fastq = tmp_path / "sample.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        for seq, n in [("A" * BARCODE_LEN, 5), ("A" * (BARCODE_LEN - 1), 3)]:
            for i in range(n):
                f.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")

    tsv = tmp_path / "counts.tsv"
    count_barcodes_python(fastq, tsv, BARCODE_LEN)
    result = {l.split("\t")[0] for l in tsv.read_text().strip().splitlines()}
    assert "A" * (BARCODE_LEN - 1) not in result


def test_count_empty_fastq(tmp_path):
    fastq = tmp_path / "empty.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        pass
    tsv = tmp_path / "counts.tsv"
    count_barcodes_python(fastq, tsv, BARCODE_LEN)
    assert tsv.read_text().strip() == ""
