from __future__ import annotations

import gzip
from itertools import groupby
from math import log10
from multiprocessing import Pool
from operator import itemgetter
from pathlib import Path
from typing import Iterator


def _majority_vote(sequences: list[str]) -> tuple[str, list[int]]:
    """Call consensus by per-position majority vote over A/C/G/T.

    Reads carrying a gap or ``N`` at a position do not vote, but they do count
    in the denominator, so the Phred-scaled quality reflects the fraction of
    *all* reads supporting the consensus base. Ties go to the first of A, C, G,
    T. A position where no read has an A/C/G/T becomes ``N`` with quality 0.
    """
    if not sequences:
        return "", []

    n_reads = len(sequences)
    consensus: list[str] = []
    qualities: list[int] = []

    for column in zip(*sequences):
        best_base, best_count = "N", 0
        for base in "ACGT":
            count = column.count(base)
            if count > best_count:
                best_base, best_count = base, count
        frac_support = best_count / n_reads
        q = int(round(-10 * log10(max(1e-6, 1.0 - frac_support))))
        consensus.append(best_base)
        qualities.append(min(40, max(0, q)))

    return "".join(consensus), qualities


def _iter_sorted_reads(reads_path: Path) -> Iterator[tuple[str, str]]:
    with gzip.open(reads_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            yield parts[0], parts[1]  # barcode, orf_seq


def _iter_barcode_groups(reads_path: Path) -> Iterator[tuple[str, list[str]]]:
    for barcode, group in groupby(_iter_sorted_reads(reads_path), key=itemgetter(0)):
        yield barcode, [orf_seq for _, orf_seq in group]


def _fastq_record(item: tuple[str, list[str]]) -> str:
    barcode, seqs = item
    consensus_seq, qualities = _majority_vote(seqs)
    q_str = "".join(chr(q + 33) for q in qualities)
    return f"@{barcode} size={len(seqs)}\n{consensus_seq}\n+\n{q_str}\n"


def call_consensus(
    reads_path: str | Path,
    output_path: str | Path,
    threads: int = 1,
) -> int:
    """Call per-barcode consensus from a sorted (barcode, orf_seq) TSV.

    Reads must be sorted by barcode (col 1) so that identical barcodes are
    consecutive — produced by ``extract_reads_sorted``. With ``threads`` > 1
    the per-barcode votes are computed in a process pool; output order and
    content are identical to the single-process run.

    Returns the number of consensus sequences written.
    """
    reads_path = Path(reads_path)
    output_path = Path(output_path)
    groups = _iter_barcode_groups(reads_path)

    n_written = 0
    with gzip.open(output_path, "wt") as out:
        if threads > 1:
            with Pool(threads) as pool:
                for record in pool.imap(_fastq_record, groups, chunksize=256):
                    out.write(record)
                    n_written += 1
        else:
            for item in groups:
                out.write(_fastq_record(item))
                n_written += 1

    return n_written
