"""Tests for pounce.reference."""

from __future__ import annotations

import pytest
from pounce.reference import parse_reference, ReferenceConfig


def test_parse_finds_n_run(plasmid_fasta, wt_orf_fasta, orf_coords):
    start, end = orf_coords
    ref = parse_reference(plasmid_fasta, wt_orf_fasta, start, end)
    assert isinstance(ref, ReferenceConfig)
    assert ref.barcode_len == 30


def test_adapter_sequences_flank_barcode(plasmid_fasta, wt_orf_fasta, orf_coords):
    from Bio import SeqIO
    start, end = orf_coords
    ref = parse_reference(plasmid_fasta, wt_orf_fasta, start, end)

    plasmid_seq = str(next(SeqIO.parse(plasmid_fasta, "fasta")).seq).upper()
    assert plasmid_seq[ref.barcode_start - 20:ref.barcode_start] == ref.adapter_5p
    assert plasmid_seq[ref.barcode_end:ref.barcode_end + 20] == ref.adapter_3p


def test_orf_length_matches(plasmid_fasta, wt_orf_fasta, orf_coords):
    from Bio import SeqIO
    start, end = orf_coords
    ref = parse_reference(plasmid_fasta, wt_orf_fasta, start, end)
    wt_len = len(next(SeqIO.parse(wt_orf_fasta, "fasta")).seq)
    assert ref.orf_len == wt_len


def test_multiple_n_runs_raises(tmp_path, wt_orf_fasta, orf_coords):
    bad = tmp_path / "bad_plasmid.fa"
    bad.write_text(">bad\nATCGNNNATCGNNNATCG\n")
    start, end = orf_coords
    with pytest.raises(ValueError, match="N-run"):
        parse_reference(bad, wt_orf_fasta, start, end)


def test_orf_length_mismatch_raises(tmp_path, plasmid_fasta, orf_coords):
    wrong_orf = tmp_path / "wrong_orf.fa"
    wrong_orf.write_text(">WT\nATCG\n")
    start, end = orf_coords
    with pytest.raises(ValueError, match="length"):
        parse_reference(plasmid_fasta, wrong_orf, start, end)
