"""Tests for pounce.longread.filter."""

from __future__ import annotations

import gzip

from pounce.longread.filter import _passes_filter, filter_consensus

DEFAULTS = dict(
    min_reads=2, min_mean_quality=10.0, max_n=5,
    low_support_reads=2, low_support_min_quality=20.0, low_support_max_n=0,
)


def passes(size, n_count=0, quality=40.0, **overrides):
    return _passes_filter(size, n_count, quality, **{**DEFAULTS, **overrides})


def test_single_read_barcodes_rejected():
    assert not passes(1)


def test_two_read_barcodes_need_agreement():
    assert passes(2, n_count=0, quality=40.0)
    assert not passes(2, n_count=1, quality=40.0)
    assert not passes(2, n_count=0, quality=19.9)


def test_three_or_more_reads_use_relaxed_thresholds():
    assert passes(3, n_count=5, quality=10.0)
    assert not passes(3, n_count=6, quality=40.0)
    assert not passes(3, n_count=0, quality=9.9)


def test_thresholds_are_parameters():
    assert not passes(3, min_reads=4)
    assert passes(2, n_count=3, quality=15.0, low_support_reads=1)


def test_filter_consensus_end_to_end(tmp_path):
    fastq = tmp_path / "consensus.fastq.gz"
    records = [
        ("bc_one_read", 1, "ACGT", "IIII"),
        ("bc_two_clean", 2, "ACGT", "IIII"),
        ("bc_two_n", 2, "ACNT", "IIII"),
        ("bc_many_ok", 8, "ANNT", "IIII"),
        ("bc_many_noisy", 8, "ACGT", "++++"),  # Q10 mean passes; checked below
        ("bc_many_lowq", 8, "ACGT", "&&&&"),    # Q5 fails
    ]
    with gzip.open(fastq, "wt") as f:
        for name, size, seq, qual in records:
            f.write(f"@{name} size={size}\n{seq}\n+\n{qual}\n")

    out = tmp_path / "filtered.fasta.gz"
    total, passed = filter_consensus(fastq, out)
    kept = [line[1:].split()[0] for line in gzip.open(out, "rt") if line.startswith(">")]

    assert total == 6
    assert kept == ["bc_two_clean", "bc_many_ok", "bc_many_noisy"]
    assert passed == 3
