#!/usr/bin/env python3
"""Generate synthetic test fixtures for pounce unit tests.

Run once: python tests/fixtures/generate.py
Outputs go to tests/fixtures/data/ and are checked in to the repo.
"""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

# Plasmid layout
FLANK_5P = "ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # 50 bp
FLANK_MID = "TCGATCGATCGATCGATCG"  # 19 bp between barcode and ORF
FLANK_3P = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"  # 50 bp
BARCODE_LEN = 30
ORF_LEN = 300  # 100 codons

BASES = list("ACGT")
random.seed(42)


# ---------------------------------------------------------------------------
# Sequencing error simulation
# ---------------------------------------------------------------------------

def illumina_q(pos: int, length: int) -> int:
    """Sample a Q-score for Illumina R2 at read position pos (0-based).

    Models the gradual quality decline toward the 3' end typical of R2 reads:
    mean declines linearly from Q35 at pos=0 to Q28 at pos=length-1.
    """
    mean_q = 35.0 - 7.0 * (pos / max(length - 1, 1))
    q = round(random.gauss(mean_q, 3))
    return max(10, min(40, q))


def pacbio_hifi_q(pos: int, length: int) -> int:
    """Sample a Q-score for PacBio HiFi.

    HiFi reads are high and roughly uniform along the read (Q30–Q40).
    """
    q = round(random.gauss(35, 3))
    return max(20, min(40, q))


def add_sequencing_errors(seq: str, q_fn) -> tuple[str, str]:
    """Add per-base errors at the rate implied by sampled Q-scores.

    For each base, a Q-score is drawn from q_fn(position, length).
    An error is introduced with probability 10^(-Q/10); the substituted
    base is chosen uniformly from the other three.

    Returns (noisy_sequence, quality_string) where quality_string is
    Phred+33 ASCII encoding of the sampled Q-scores.
    """
    n = len(seq)
    noisy = list(seq)
    quals: list[int] = []
    for i, base in enumerate(seq):
        q = q_fn(i, n)
        quals.append(q)
        if random.random() < 10 ** (-q / 10):
            noisy[i] = random.choice([b for b in BASES if b != base])
    qual_str = "".join(chr(q + 33) for q in quals)
    return "".join(noisy), qual_str


# ---------------------------------------------------------------------------
# Reference / ORF helpers
# ---------------------------------------------------------------------------

def _build_codon_table() -> dict[str, str]:
    from Bio.Data import CodonTable
    table = CodonTable.standard_dna_table
    d = dict(table.forward_table)
    for stop in table.stop_codons:
        d[stop] = "*"
    return d


def random_seq(n: int) -> str:
    return "".join(random.choices(BASES, k=n))


def make_wt_orf(length: int = ORF_LEN) -> str:
    """Random ORF with no internal stop codons."""
    codon_table = _build_codon_table()
    sense_codons = [c for c, aa in codon_table.items() if aa != "*"]
    codons = [random.choice(sense_codons) for _ in range(length // 3)]
    return "".join(codons)


def mutate_seq(seq: str, n_muts: int) -> tuple[str, list[str]]:
    """Apply n_muts random SNPs; return (mutated_seq, ['REFposALT', ...])."""
    positions = random.sample(range(len(seq)), n_muts)
    seq_list = list(seq)
    mutations = []
    for pos in sorted(positions):
        wt = seq_list[pos]
        alt = random.choice([b for b in BASES if b != wt])
        mutations.append(f"{wt}{pos + 1}{alt}")
        seq_list[pos] = alt
    return "".join(seq_list), mutations


def nuc_muts_to_aa(nt_muts: list[str], wt_seq: str) -> list[str]:
    from Bio.Seq import Seq

    if not nt_muts:
        return []
    wt = list(wt_seq)
    mut = list(wt_seq)
    parsed = [(int(m[1:-1]) - 1, m[0], m[-1]) for m in nt_muts]
    for pos, _, alt in parsed:
        mut[pos] = alt
    affected = sorted({pos // 3 for pos, _, _ in parsed})
    aa_muts = []
    for ci in affected:
        wt_codon = "".join(wt[ci * 3:(ci + 1) * 3])
        mut_codon = "".join(mut[ci * 3:(ci + 1) * 3])
        wt_aa = str(Seq(wt_codon).translate())
        mut_aa = str(Seq(mut_codon).translate())
        if wt_aa != mut_aa:
            aa_muts.append(f"{wt_aa}{ci + 1}{mut_aa}")
    return aa_muts


# ---------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------

def make_fixtures() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # 1. WT ORF and plasmid
    wt_orf = make_wt_orf(ORF_LEN)
    (DATA_DIR / "test_wt_orf.fa").write_text(f">WT_ORF\n{wt_orf}\n")

    plasmid_seq = FLANK_5P + ("N" * BARCODE_LEN) + FLANK_MID + wt_orf + FLANK_3P
    (DATA_DIR / "test_plasmid.fa").write_text(f">test_plasmid\n{plasmid_seq}\n")

    orf_start = len(FLANK_5P) + BARCODE_LEN + len(FLANK_MID)
    orf_end = orf_start + ORF_LEN
    (DATA_DIR / "test_coords.txt").write_text(f"{orf_start}\n{orf_end}\n")

    adapter_5p = FLANK_5P[-20:]
    adapter_3p = FLANK_MID[:20]

    # 2. 100 barcodes with assigned mutations (ground truth)
    barcodes: list[str] = []
    seen: set[str] = set()
    while len(barcodes) < 100:
        bc = random_seq(BARCODE_LEN)
        if bc not in seen:
            seen.add(bc)
            barcodes.append(bc)

    truth_rows = []
    bc_to_orf: dict[str, str] = {}
    for bc in barcodes:
        n_muts = random.choices([0, 1, 2, 3], weights=[30, 50, 15, 5])[0]
        if n_muts == 0:
            bc_to_orf[bc] = wt_orf
            truth_rows.append({"barcode": bc, "nt_mutations": "", "aa_mutations": ""})
        else:
            mut_orf, nt_muts = mutate_seq(wt_orf, n_muts)
            aa_muts = nuc_muts_to_aa(nt_muts, wt_orf)
            bc_to_orf[bc] = mut_orf
            truth_rows.append({
                "barcode": bc,
                "nt_mutations": ";".join(nt_muts),
                "aa_mutations": ";".join(aa_muts),
            })

    with open(DATA_DIR / "test_barcodes_truth.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["barcode", "nt_mutations", "aa_mutations"])
        w.writeheader()
        w.writerows(truth_rows)

    # 3. Short-read R2 FASTQs with realistic Illumina errors
    ns_counts = {bc: random.randint(20, 60) for bc in barcodes}
    s_counts = {
        bc: max(1, ns_counts[bc] * (5 if i < 20 else 1) + random.randint(-2, 2))
        for i, bc in enumerate(barcodes)
    }

    def _make_sample_reads(bc_counts: dict[str, int], fastq_path: Path) -> None:
        with gzip.open(fastq_path, "wt") as f:
            read_idx = 0
            for bc, n in bc_counts.items():
                true_seq = adapter_5p + bc + adapter_3p
                for _ in range(n):
                    noisy_seq, qual_str = add_sequencing_errors(true_seq, illumina_q)
                    f.write(f"@read{read_idx}\n{noisy_seq}\n+\n{qual_str}\n")
                    read_idx += 1

    _make_sample_reads(ns_counts, DATA_DIR / "test_NS_A_R2.fastq.gz")
    _make_sample_reads(
        {bc: max(1, c + random.randint(-5, 5)) for bc, c in ns_counts.items()},
        DATA_DIR / "test_NS_B_R2.fastq.gz",
    )
    _make_sample_reads(s_counts, DATA_DIR / "test_S_A_R2.fastq.gz")
    _make_sample_reads(
        {bc: max(1, c + random.randint(-5, 5)) for bc, c in s_counts.items()},
        DATA_DIR / "test_S_B_R2.fastq.gz",
    )

    # 4. Pre-built barcode count TSVs (bypass trim+count for unit tests)
    for sid, counts_dict in [
        ("SP01", ns_counts), ("SP02", ns_counts),
        ("SP03", s_counts),  ("SP04", s_counts),
    ]:
        with open(DATA_DIR / f"{sid}_barcode_counts.tsv", "w") as f:
            for bc, cnt in sorted(counts_dict.items()):
                f.write(f"{bc}\t{cnt}\n")

    # 5. Plate info CSV
    with open(DATA_DIR / "test_plate_info.csv", "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["sample_id", "fastq_r2", "timepoint", "replicate", "condition"]
        )
        w.writeheader()
        for sid, tp, rep, cond, fname in [
            ("SP01", "T1", "A", "NS", "test_NS_A_R2.fastq.gz"),
            ("SP02", "T1", "B", "NS", "test_NS_B_R2.fastq.gz"),
            ("SP03", "T1", "A", "S",  "test_S_A_R2.fastq.gz"),
            ("SP04", "T1", "B", "S",  "test_S_B_R2.fastq.gz"),
        ]:
            w.writerow({
                "sample_id": sid,
                "fastq_r2": str(DATA_DIR / fname),
                "timepoint": tp,
                "replicate": rep,
                "condition": cond,
            })

    # 6. Long-read sorted TSV for consensus testing (PacBio HiFi error model)
    #    Format: barcode<TAB>orf_seq, sorted by barcode, 5-15 reads per barcode
    reads_per_bc = {bc: random.randint(5, 15) for bc in barcodes}
    rows: list[tuple[str, str]] = []
    for bc in barcodes:
        true_orf = bc_to_orf[bc]
        for _ in range(reads_per_bc[bc]):
            noisy_orf, _ = add_sequencing_errors(true_orf, pacbio_hifi_q)
            rows.append((bc, noisy_orf))

    rows.sort(key=lambda x: x[0])
    with gzip.open(DATA_DIR / "test_longread_reads.tsv.gz", "wt") as f:
        for bc, orf_seq in rows:
            f.write(f"{bc}\t{orf_seq}\n")

    total_reads = sum(reads_per_bc.values())
    print(f"Generated fixtures in {DATA_DIR}")
    print(f"  orf_start={orf_start}, orf_end={orf_end}")
    print(f"  adapter_5p={adapter_5p}")
    print(f"  adapter_3p={adapter_3p}")
    print(f"  100 barcodes | {sum(ns_counts.values())} NS reads | {sum(s_counts.values())} S reads")
    print(f"  {total_reads} long reads across 100 barcodes (5–15 per barcode, PacBio HiFi errors)")


if __name__ == "__main__":
    make_fixtures()
