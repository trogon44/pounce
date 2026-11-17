from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterator

import pysam

from pounce.reference import ReferenceConfig


def _is_aligned(bam_path: Path) -> bool:
    """Return True if the BAM header declares at least one reference sequence."""
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        return len(bam.header.references) > 0


def _align_to_plasmid(bam_path: Path, plasmid_path: Path, work_dir: Path, threads: int) -> Path:
    """Align an unaligned PacBio HiFi/CCS BAM to the plasmid reference.

    Raw instrument BAMs (e.g. hifi_reads/*.bam) have no @SQ header and are not
    aligned to anything. Interval extraction needs reference-coordinate
    alignment, so we align with minimap2 (map-hifi preset) and produce a
    coordinate-sorted, indexed BAM in work_dir.
    """
    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"'{tool}' not found on PATH; it is required to align the unaligned "
                "input BAM to the plasmid reference before extraction."
            )

    aligned_bam = work_dir / f"{bam_path.stem}.aligned_to_plasmid.bam"
    cmd = (
        f"samtools fastq -@ {threads} {bam_path} | "
        f"minimap2 -t {threads} -ax map-hifi {plasmid_path} - | "
        f"samtools sort -@ {threads} -o {aligned_bam} -"
    )
    result = subprocess.run(
        f"set -o pipefail; {cmd}", shell=True, executable="/bin/bash",
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Alignment to plasmid reference failed: {result.stderr}")

    subprocess.run(["samtools", "index", "-@", str(threads), str(aligned_bam)], check=True)
    return aligned_bam


def _ensure_aligned_bam(bam_path: Path, plasmid_path: Path, work_dir: Path, threads: int) -> Path:
    """Return a path to an indexed BAM aligned to plasmid_path.

    If bam_path is already aligned, index it in place if needed and return it
    unchanged. Otherwise, align it to the plasmid reference first.
    """
    if not _is_aligned(bam_path):
        return _align_to_plasmid(bam_path, plasmid_path, work_dir, threads)

    has_index = bam_path.with_suffix(bam_path.suffix + ".bai").exists() or Path(
        str(bam_path) + ".bai"
    ).exists()
    if not has_index:
        pysam.index(str(bam_path))
    return bam_path


# CIGAR operations that consume the query / the reference (pysam codes).
_ALIGNED = (0, 7, 8)   # M, =, X
_INSERTION = 1
_DELETION = 2
_REF_SKIP = 3
_SOFT_CLIP = 4


def _extract_interval(
    read: pysam.AlignedSegment,
    ref_start: int,
    ref_end: int,
    allow_partial: bool = False,
) -> str | None:
    """Extract read bases covering [ref_start, ref_end) in reference coordinates.

    Returns None if the read does not fully span the interval, unless
    ``allow_partial`` is set, in which case the part of the interval outside
    the read's alignment is filled with 'N' (None is still returned if the read
    does not overlap the interval at all). Deletions relative to the reference
    are filled with 'N'; insertions are dropped. Walks the CIGAR instead of
    building per-base aligned pairs, which is much faster for the multi-kb
    HiFi reads this is run on millions of times.
    """
    lo_bound = max(read.reference_start, ref_start)
    hi_bound = min(read.reference_end, ref_end)
    if lo_bound >= hi_bound:
        return None
    if not allow_partial and (lo_bound != ref_start or hi_bound != ref_end):
        return None

    seq = read.query_sequence
    q_pos = 0
    r_pos = read.reference_start
    pieces: list[str] = []

    for op, length in read.cigartuples:
        if op in _ALIGNED:
            lo, hi = max(r_pos, lo_bound), min(r_pos + length, hi_bound)
            if lo < hi:
                pieces.append(seq[q_pos + lo - r_pos:q_pos + hi - r_pos])
            q_pos += length
            r_pos += length
        elif op == _DELETION:
            lo, hi = max(r_pos, lo_bound), min(r_pos + length, hi_bound)
            if lo < hi:
                pieces.append("N" * (hi - lo))
            r_pos += length
        elif op == _REF_SKIP:
            if max(r_pos, lo_bound) < min(r_pos + length, hi_bound):
                return None
            r_pos += length
        elif op in (_INSERTION, _SOFT_CLIP):
            q_pos += length
        if r_pos >= hi_bound:
            break

    body = "".join(pieces)
    if len(body) != hi_bound - lo_bound:
        return None
    return "N" * (lo_bound - ref_start) + body + "N" * (ref_end - hi_bound)


def _iter_records(
    bam_path: Path,
    ref_config: ReferenceConfig,
    mapq_min: int,
    threads: int = 1,
    allow_partial_orf: bool = True,
) -> Iterator[tuple[str, str, str]]:
    """Yield (barcode, orf_seq, read_id) for each usable read.

    The barcode must be fully spanned with no deletions. With
    ``allow_partial_orf`` a read that covers only part of the ORF still
    contributes, with ``N`` for the uncovered bases.
    """
    with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
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
                read, ref_config.orf_start, ref_config.orf_end, allow_partial_orf
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
    threads: int = 4,
    allow_partial_orf: bool = True,
) -> None:
    """Extract (barcode, orf_seq, read_id) from BAM and write sorted by barcode.

    Output is a gzip-compressed TSV with columns: barcode, orf_seq, read_id.
    Sorting by barcode groups all reads for the same barcode consecutively,
    enabling single-pass consensus calling.
    """
    bam_path = Path(bam_path)
    output_path = Path(output_path)

    # The aligned BAM of a full PacBio cell is tens of GB: keep it beside the output, not in /tmp.
    with tempfile.TemporaryDirectory(dir=output_path.parent) as align_dir:
        aligned_bam_path = _ensure_aligned_bam(bam_path, ref_config.plasmid_path, Path(align_dir), threads)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, dir=output_path.parent
        ) as tmp:
            tmp_path = Path(tmp.name)
            n_written = 0
            for barcode, orf_seq, read_id in _iter_records(
                aligned_bam_path, ref_config, mapq_min, threads, allow_partial_orf
            ):
                tmp.write(f"{barcode}\t{orf_seq}\t{read_id}\n")
                n_written += 1

    if n_written == 0:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(
            "No reads passed extraction filters. Check BAM alignment and reference coordinates."
        )

    result = subprocess.run(
        f"LC_ALL=C sort -k1,1 --buffer-size={sort_buffer} -T {output_path.parent} {tmp_path} | gzip > {output_path}",
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    tmp_path.unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"sort failed: {result.stderr}")
