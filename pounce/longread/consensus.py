from __future__ import annotations

import gzip
from collections import defaultdict
from itertools import groupby
from math import log10
from pathlib import Path
from typing import Iterator


def _majority_vote(sequences: list[str]) -> tuple[str, list[int]]:
    """Call consensus by per-position majority vote.

    Returns (consensus_sequence, phred_quality_scores).
    Quality reflects the fraction of reads supporting the consensus base.
    """
    if not sequences:
        return "", []

    length = len(sequences[0])
    consensus: list[str] = []
    qualities: list[int] = []

    for pos in range(length):
        counts: dict[str, int] = defaultdict(int)
        for seq in sequences:
            if pos < len(seq):
                base = seq[pos]
                if base in "ACGT":
                    counts[base] += 1

        if not counts:
            consensus.append("N")
            qualities.append(0)
            continue

        best_base = max(counts, key=counts.__getitem__)
        best_count = counts[best_base]
        total = sum(counts.values())

        frac_support = best_count / total
        q = min(40, max(0, int(-10 * log10(max(1e-6, 1.0 - frac_support)))))

        consensus.append(best_base)
        qualities.append(q)

    return "".join(consensus), qualities


def _iter_sorted_reads(reads_path: Path) -> Iterator[tuple[str, str]]:
    with gzip.open(reads_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            yield parts[0], parts[1]  # barcode, orf_seq


def call_consensus(
    reads_path: str | Path,
    output_path: str | Path,
) -> int:
    """Call per-barcode consensus from a sorted (barcode, orf_seq) TSV.

    Reads must be sorted by barcode (col 1) so that identical barcodes are
    consecutive — produced by ``extract_reads_sorted``.

    Returns the number of consensus sequences written.
    """
    reads_path = Path(reads_path)
    output_path = Path(output_path)

    n_written = 0
    with gzip.open(output_path, "wt") as out:
        for barcode, group_iter in groupby(_iter_sorted_reads(reads_path), key=lambda x: x[0]):
            seqs = [orf_seq for _, orf_seq in group_iter]
            consensus_seq, qualities = _majority_vote(seqs)
            q_str = "".join(chr(q + 33) for q in qualities)
            out.write(f"@{barcode} size={len(seqs)}\n{consensus_seq}\n+\n{q_str}\n")
            n_written += 1

    return n_written
