"""Tests for pounce.longread.extract interval extraction."""

from __future__ import annotations

import random

import pysam
import pytest

from pounce.longread.extract import _extract_interval

REF_LEN = 400


def make_read(cigar: list[tuple[int, int]], ref_start: int, seq: str) -> pysam.AlignedSegment:
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref", "LN": REF_LEN}]}
    )
    read = pysam.AlignedSegment(header)
    read.query_name = "r"
    read.query_sequence = seq
    read.reference_id = 0
    read.reference_start = ref_start
    read.cigartuples = cigar
    read.mapping_quality = 60
    return read


def reference_walk(read: pysam.AlignedSegment, ref_start: int, ref_end: int) -> str | None:
    """Slow but obviously correct version built on get_aligned_pairs."""
    by_ref = {r: q for q, r in read.get_aligned_pairs(matches_only=False) if r is not None}
    out = []
    for pos in range(ref_start, ref_end):
        if pos not in by_ref:
            return None
        q = by_ref[pos]
        out.append("N" if q is None else read.query_sequence[q])
    return "".join(out)


def random_alignment(rng: random.Random):
    ops: list[tuple[int, int]] = []
    q_len = 0
    if rng.random() < 0.5:
        n = rng.randint(1, 20)
        ops.append((4, n))
        q_len += n
    for _ in range(rng.randint(2, 12)):
        kind = rng.choices([0, 1, 2], weights=[6, 1, 1])[0]
        n = rng.randint(1, 30) if kind == 0 else rng.randint(1, 4)
        ops.append((kind, n))
        if kind in (0, 1):
            q_len += n
    if ops[-1][0] != 0:
        ops.append((0, 5))
        q_len += 5
    if ops[0][0] in (1, 2):
        ops.insert(0, (0, 5))
        q_len += 5
    # collapse adjacent identical ops is unnecessary; the aligner would not emit I next to D but pysam does not care
    seq = "".join(rng.choice("ACGT") for _ in range(q_len))
    return ops, seq


def test_simple_match():
    read = make_read([(0, 100)], 50, "ACGT" * 25)
    assert _extract_interval(read, 60, 70) == ("ACGT" * 25)[10:20]


def test_deletion_filled_with_n():
    read = make_read([(0, 10), (2, 3), (0, 10)], 0, "A" * 10 + "C" * 10)
    assert _extract_interval(read, 8, 15) == "AANNNCC"


def test_insertion_dropped():
    read = make_read([(0, 5), (1, 2), (0, 5)], 0, "AAAAATTCCCCC")
    assert _extract_interval(read, 3, 8) == "AACCC"


def test_soft_clip_skipped():
    read = make_read([(4, 4), (0, 6)], 10, "TTTTACGTAC")
    assert _extract_interval(read, 10, 14) == "ACGT"


@pytest.mark.parametrize("start,end", [(0, 10), (45, 70), (95, 110)])
def test_read_not_spanning_interval_is_none(start, end):
    read = make_read([(0, 50)], 50, "A" * 50)
    assert _extract_interval(read, start, end) is None


def test_partial_interval_padded_with_n_when_allowed():
    read = make_read([(0, 10)], 5, "ACGTACGTAC")        # covers ref 5..15
    assert _extract_interval(read, 0, 20) is None
    assert _extract_interval(read, 0, 20, allow_partial=True) == "N" * 5 + "ACGTACGTAC" + "N" * 5
    assert _extract_interval(read, 10, 20, allow_partial=True) == "CGTAC" + "N" * 5
    assert _extract_interval(read, 15, 25, allow_partial=True) is None  # no overlap at all


def test_partial_interval_with_deletion_and_insertion():
    read = make_read([(0, 4), (2, 2), (1, 3), (0, 4)], 10, "ACGT" + "TTT" + "ACGT")  # ref 10..20
    assert _extract_interval(read, 8, 22, allow_partial=True) == "NN" + "ACGT" + "NN" + "ACGT" + "NN"


def test_matches_get_aligned_pairs_on_random_alignments():
    rng = random.Random(7)
    checked = 0
    for _ in range(500):
        ops, seq = random_alignment(rng)
        ref_start = rng.randint(0, 50)
        read = make_read(ops, ref_start, seq)
        if read.reference_end > REF_LEN:
            continue
        a = rng.randint(ref_start - 5, read.reference_end)
        b = a + rng.randint(1, 40)
        assert _extract_interval(read, a, b) == reference_walk(read, a, b), (ops, a, b)
        checked += 1
    assert checked > 400
