from __future__ import annotations

import gzip
from pathlib import Path

from Bio import SeqIO


def _passes_filter(
    size: int,
    n_count: int,
    mean_quality: float,
    min_reads: int,
    min_mean_quality: float,
    max_n: int,
    low_support_reads: int,
    low_support_min_quality: float,
    low_support_max_n: int,
) -> bool:
    if size < min_reads:
        return False
    if size <= low_support_reads:
        if mean_quality < low_support_min_quality or n_count > low_support_max_n:
            return False
    return mean_quality >= min_mean_quality and n_count <= max_n


def filter_consensus(
    consensus_fastq: str | Path,
    output_fasta: str | Path,
    min_reads: int = 2,
    min_mean_quality: float = 10.0,
    max_n: int = 5,
    low_support_reads: int = 2,
    low_support_min_quality: float = 20.0,
    low_support_max_n: int = 0,
) -> tuple[int, int]:
    """Filter consensus FASTQ and write passing records as FASTA.

    A barcode is kept if it has at least ``min_reads`` reads, mean consensus
    quality >= ``min_mean_quality`` and at most ``max_n`` ``N`` bases. Barcodes
    supported by ``low_support_reads`` reads or fewer must additionally meet
    the stricter ``low_support_min_quality`` and ``low_support_max_n``, since
    with so few reads a consensus is only trustworthy if the reads agree
    completely. The defaults are the thresholds used for the TasR ePCR library.

    Returns (n_total, n_passed).
    """
    consensus_fastq = Path(consensus_fastq)
    output_fasta = Path(output_fasta)

    n_total = 0
    n_passed = 0

    with gzip.open(consensus_fastq, "rt") as handle, gzip.open(output_fasta, "wt") as out:
        for record in SeqIO.parse(handle, "fastq"):
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

            if _passes_filter(
                size, sequence.count("N"), mean_q, min_reads, min_mean_quality, max_n,
                low_support_reads, low_support_min_quality, low_support_max_n,
            ):
                out.write(f">{record.id} size={size}\n{sequence}\n")
                n_passed += 1

    return n_total, n_passed
