from __future__ import annotations

import gzip
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _passes_filter(
    size: int,
    sequence: str,
    mean_quality: float,
    min_reads: int,
    max_n_frac: float,
    min_mean_quality: float,
) -> bool:
    if size < min_reads:
        return False
    n_frac = sequence.count("N") / max(len(sequence), 1)
    if n_frac > max_n_frac:
        return False
    if mean_quality < min_mean_quality:
        return False
    return True


def filter_consensus(
    consensus_fastq: str | Path,
    output_fasta: str | Path,
    min_reads: int = 3,
    max_n_frac: float = 0.01,
    min_mean_quality: float = 20.0,
) -> tuple[int, int]:
    """Filter consensus FASTQ and write passing records as FASTA.

    Returns (n_total, n_passed).
    """
    consensus_fastq = Path(consensus_fastq)
    output_fasta = Path(output_fasta)

    n_total = 0
    n_passed = 0

    with gzip.open(output_fasta, "wt") as out:
        for record in SeqIO.parse(gzip.open(consensus_fastq, "rt"), "fastq"):
            n_total += 1

            # Header format: "@{barcode} size={N}"
            size_str = record.description.split("size=")[-1]
            try:
                size = int(size_str)
            except ValueError:
                continue

            qualities = record.letter_annotations["phred_quality"]
            mean_q = sum(qualities) / len(qualities) if qualities else 0.0
            sequence = str(record.seq)

            if _passes_filter(size, sequence, mean_q, min_reads, max_n_frac, min_mean_quality):
                out.write(f">{record.id} size={size}\n{sequence}\n")
                n_passed += 1

    return n_total, n_passed
