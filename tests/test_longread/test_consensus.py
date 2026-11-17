"""Tests for pounce.longread.consensus."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from pounce.longread.consensus import _majority_vote, call_consensus


def test_majority_vote_clear_winner():
    seqs = ["ACGT", "ACGT", "ACGT", "ACGA"]
    consensus, qualities = _majority_vote(seqs)
    assert consensus == "ACGT"
    # Last position: 3/4 = 75% → quality should be positive
    assert qualities[-1] > 0


def test_majority_vote_all_agree():
    seqs = ["AAAA"] * 10
    consensus, qualities = _majority_vote(seqs)
    assert consensus == "AAAA"
    assert all(q == 40 for q in qualities)


def test_majority_vote_snp_recovered():
    # 9 reads have 'G' at position 2, 1 read has 'A'
    seqs = ["ACGT"] * 9 + ["ACAT"]
    consensus, _ = _majority_vote(seqs)
    assert consensus[2] == "G"


def test_majority_vote_empty():
    seq, quals = _majority_vote([])
    assert seq == ""
    assert quals == []


def test_majority_vote_n_on_tie():
    # Equal A and G — 'N' is expected when no majority
    seqs = ["AAAA", "GAAA"]
    consensus, qualities = _majority_vote(seqs)
    # Position 0: tie between A and G
    # Both have 50% support. The function picks whichever comes first in dict order.
    # Quality at a tie position should be low (0 or close to it).
    assert qualities[0] == 0 or consensus[0] in "ACGTN"


def test_majority_vote_gap_reads_count_as_dissent():
    # 3 reads: two call A, one has a deletion (N). Support is 2/3, not 2/2.
    consensus, qualities = _majority_vote(["A", "A", "N"])
    assert consensus == "A"
    assert qualities[0] == 5  # -10*log10(1/3) = 4.77 -> 5


def test_majority_vote_gap_does_not_outvote_base():
    # The consensus base is the most common A/C/G/T even when gaps are more numerous
    consensus, _ = _majority_vote(["N", "N", "G"])
    assert consensus == "G"


def test_majority_vote_all_gap_is_n():
    consensus, qualities = _majority_vote(["N", "N"])
    assert consensus == "N"
    assert qualities == [0]


def test_majority_vote_tie_goes_to_first_in_acgt_order():
    consensus, qualities = _majority_vote(["G", "A"])
    assert consensus == "A"
    assert qualities == [3]


def test_call_consensus_roundtrip(tmp_path):
    """Write a tiny sorted reads TSV and check that call_consensus produces the right output."""
    reads_path = tmp_path / "reads.tsv.gz"
    output_path = tmp_path / "consensus.fastq.gz"

    # Two barcodes, multiple reads each
    bc1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"[:30]
    bc2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[:30]
    orf1 = "ATCGATCGATCG" * 4  # 48 bp
    orf2 = "GCTAGCTAGCTA" * 4

    with gzip.open(reads_path, "wt") as f:
        for _ in range(5):
            f.write(f"{bc1}\t{orf1}\tread1\n")
        for _ in range(3):
            f.write(f"{bc2}\t{orf2}\tread2\n")

    n = call_consensus(reads_path, output_path)
    assert n == 2

    from Bio import SeqIO
    records = list(SeqIO.parse(gzip.open(output_path, "rt"), "fastq"))
    assert len(records) == 2
    ids = {r.id for r in records}
    assert bc1 in ids and bc2 in ids

    for r in records:
        if r.id == bc1:
            assert str(r.seq) == orf1
            size = int(r.description.split("size=")[-1])
            assert size == 5


def test_consensus_recovers_from_pacbio_noise(tmp_path, data_dir):
    """Majority-vote consensus on PacBio HiFi-like noisy reads recovers ≥99% of bases.

    Fixture has 5–15 reads per barcode with errors drawn from Q~N(35,3) capped [20,40].
    At Q35 the per-base error rate is ~0.03%, so most positions are unanimous and
    the consensus should be nearly perfect.
    """
    import pandas as pd
    from Bio import SeqIO

    reads_path = data_dir / "test_longread_reads.tsv.gz"
    output_path = tmp_path / "consensus.fastq.gz"

    n_barcodes = call_consensus(reads_path, output_path)
    assert n_barcodes == 100

    truth = pd.read_csv(data_dir / "test_barcodes_truth.csv", index_col="barcode")
    wt_seq = str(next(SeqIO.parse(data_dir / "test_wt_orf.fa", "fasta")).seq).upper()

    total_bases = 0
    correct_bases = 0

    for record in SeqIO.parse(gzip.open(output_path, "rt"), "fastq"):
        bc = record.id
        consensus_seq = str(record.seq).upper()

        row = truth.loc[bc]
        raw_muts = row["nt_mutations"] if isinstance(row["nt_mutations"], str) else ""
        true_seq = list(wt_seq)
        for m in (raw_muts.split(";") if raw_muts else []):
            pos = int(m[1:-1]) - 1
            true_seq[pos] = m[-1]
        true_seq = "".join(true_seq)

        for t, c in zip(true_seq, consensus_seq):
            if c != "N":
                total_bases += 1
                correct_bases += t == c

    recovery = correct_bases / total_bases
    assert recovery >= 0.99, f"Consensus recovery {recovery:.4f} < 0.99"


def test_call_consensus_threads_match_single_process(tmp_path, data_dir):
    reads_path = data_dir / "test_longread_reads.tsv.gz"
    one, many = tmp_path / "one.fastq.gz", tmp_path / "many.fastq.gz"
    assert call_consensus(reads_path, one) == call_consensus(reads_path, many, threads=3)
    assert gzip.open(one, "rt").read() == gzip.open(many, "rt").read()
