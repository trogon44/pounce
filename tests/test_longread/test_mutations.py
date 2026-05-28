"""Tests for pounce.longread.mutations."""

from __future__ import annotations

import csv
import gzip
import tempfile
from pathlib import Path

from pounce.longread.mutations import (
    call_nt_mutations,
    nuc_to_aa_mutations,
    call_all_mutations,
)

WT = "ATGAAAGTTGCTATCGACTGGGCCTTTGCTGAAGAGAAGTTGGCTGTC"  # 48 bp, 16 codons


def test_call_nt_mutations_none():
    assert call_nt_mutations(WT, WT) == []


def test_call_nt_mutations_single():
    mut = list(WT)
    mut[3] = "C"  # A4C (1-based position 4)
    muts = call_nt_mutations("".join(mut), WT)
    assert muts == ["A4C"]


def test_call_nt_mutations_multiple():
    mut = list(WT)
    mut[0] = "C"  # A1C
    mut[9] = "T"  # A10T (if WT[9]='A', else different)
    wt_base_9 = WT[9]
    muts = call_nt_mutations("".join(mut), WT)
    assert "A1C" in muts
    assert any(m.endswith("T") and m[1:-1] == "10" for m in muts)


def test_call_nt_mutations_n_ignored():
    # 'N' in consensus should not be called as a mutation
    mut = list(WT)
    mut[3] = "N"
    muts = call_nt_mutations("".join(mut), WT)
    assert not any("4" in m for m in muts)


def test_nuc_to_aa_synonymous():
    # Change the 3rd position of codon 1 (ATG → ATG synonymous not possible; use codon 2)
    # Codon 2: positions 3,4,5 = AAA (K); change pos 5 (A6G) → AAG (still K)
    assert WT[3:6] == "AAA"
    syn_mut = list(WT)
    syn_mut[5] = "G"  # A6G: AAA → AAG, both = Lys
    aa_muts = nuc_to_aa_mutations(["A6G"], "".join(syn_mut)[:6] + WT[6:])
    # Should be empty (synonymous)
    assert aa_muts == []


def test_nuc_to_aa_missense():
    # Codon 2: AAA (K) → change pos 4 (A4C) → CAA (Q) — K→Q missense
    mut_seq = list(WT)
    mut_seq[3] = "C"  # A4C
    aa_muts = nuc_to_aa_mutations(["A4C"], WT)
    assert len(aa_muts) == 1
    assert aa_muts[0].startswith("K")  # WT AA at codon 2


def test_nuc_to_aa_stop():
    # Introduce a stop codon at codon 2: AAA → TAA
    mut_seq = list(WT)
    mut_seq[3] = "T"  # A4T → TAA = stop
    aa_muts = nuc_to_aa_mutations(["A4T"], WT)
    assert any("*" in m for m in aa_muts)


def test_call_all_mutations_roundtrip(tmp_path, data_dir):
    """call_all_mutations on truth barcodes recovers expected AA mutations."""
    import pandas as pd

    filtered_fasta = tmp_path / "filtered.fasta.gz"
    wt_orf = data_dir / "test_wt_orf.fa"

    from Bio import SeqIO
    wt_seq = str(next(SeqIO.parse(wt_orf, "fasta")).seq).upper()

    truth = pd.read_csv(data_dir / "test_barcodes_truth.csv", index_col="barcode")

    # Write a tiny filtered fasta from truth
    with gzip.open(filtered_fasta, "wt") as f:
        for bc, row in truth.head(10).iterrows():
            raw = row["nt_mutations"] if isinstance(row["nt_mutations"], str) else ""
            nt_muts = [m for m in raw.split(";") if m]
            seq = list(wt_seq)
            for m in nt_muts:
                pos = int(m[1:-1]) - 1
                seq[pos] = m[-1]
            f.write(f">{bc} size=5\n{''.join(seq)}\n")

    output_csv = tmp_path / "barcode_variants.csv"
    n = call_all_mutations(filtered_fasta, wt_orf, output_csv)
    assert n == 10

    result = pd.read_csv(output_csv, index_col="barcode")
    for bc, row in truth.head(10).iterrows():
        expected_nt = "" if pd.isna(row["nt_mutations"]) else row["nt_mutations"]
        got_nt = "" if pd.isna(result.loc[bc, "nt_mutations"]) else result.loc[bc, "nt_mutations"]
        assert got_nt == expected_nt, f"Barcode {bc}: expected {expected_nt!r}, got {got_nt!r}"
