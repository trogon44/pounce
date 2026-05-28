from __future__ import annotations

import gzip
import subprocess
import tempfile
from pathlib import Path
from typing import Iterator

import pysam

from pounce.reference import ReferenceConfig


def _extract_interval(
    read: pysam.AlignedSegment,
    ref_start: int,
    ref_end: int,
) -> str | None:
    """Extract read bases covering [ref_start, ref_end) in reference coordinates.

    Returns None if the read does not fully span the interval.
    Deletions relative to the reference are filled with 'N'.
    """
    pairs = read.get_aligned_pairs(matches_only=False)
    result: list[str] = []
    covered_ref_positions: set[int] = set()

    for ref_pos, read_pos in pairs:
        if ref_pos is None or ref_pos < ref_start:
            continue
        if ref_pos >= ref_end:
            break
        covered_ref_positions.add(ref_pos)
        result.append("N" if read_pos is None else read.query_sequence[read_pos])

    expected = set(range(ref_start, ref_end))
    if covered_ref_positions != expected:
        return None

    return "".join(result)


def _iter_records(
    bam_path: Path,
    ref_config: ReferenceConfig,
    mapq_min: int,
) -> Iterator[tuple[str, str, str]]:
    """Yield (barcode, orf_seq, read_id) for each usable read."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if read.mapping_quality < mapq_min:
                continue
            if read.query_sequence is None:
                continue

            barcode = _extract_interval(
                read, ref_config.barcode_start, ref_config.barcode_end
            )
            if barcode is None or "N" in barcode:
                continue

            orf_seq = _extract_interval(
                read, ref_config.orf_start, ref_config.orf_end
            )
            if orf_seq is None:
                continue

            yield barcode, orf_seq, read.query_name


def extract_reads_sorted(
    bam_path: str | Path,
    ref_config: ReferenceConfig,
    output_path: str | Path,
    mapq_min: int = 0,
    sort_buffer: str = "2G",
) -> None:
    """Extract (barcode, orf_seq, read_id) from BAM and write sorted by barcode.

    Output is a gzip-compressed TSV with columns: barcode, orf_seq, read_id.
    Sorting by barcode groups all reads for the same barcode consecutively,
    enabling single-pass consensus calling.
    """
    bam_path = Path(bam_path)
    output_path = Path(output_path)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        n_written = 0
        for barcode, orf_seq, read_id in _iter_records(bam_path, ref_config, mapq_min):
            tmp.write(f"{barcode}\t{orf_seq}\t{read_id}\n")
            n_written += 1

    if n_written == 0:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(
            "No reads passed extraction filters. Check BAM alignment and reference coordinates."
        )

    result = subprocess.run(
        f"LC_ALL=C sort -k1,1 --buffer-size={sort_buffer} {tmp_path} | gzip > {output_path}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    tmp_path.unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"sort failed: {result.stderr}")
