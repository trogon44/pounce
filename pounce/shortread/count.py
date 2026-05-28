from __future__ import annotations

import subprocess
from collections import Counter
from pathlib import Path


def count_barcodes_shell(
    barcodes_fastq: str | Path,
    output_tsv: str | Path,
    barcode_len: int,
    sort_buffer: str = "4G",
    sort_parallel: int = 0,
) -> None:
    """Count barcodes using awk/sort/uniq pipeline.

    Preferred for large files (>500 MB). Output TSV: barcode<TAB>count, no header.
    """
    barcodes_fastq = Path(barcodes_fastq)

    decompressor = "zcat" if str(barcodes_fastq).endswith(".gz") else "cat"
    sort_opts = f"LC_ALL=C sort --buffer-size={sort_buffer}"
    if sort_parallel > 0:
        sort_opts += f" --parallel={sort_parallel}"

    cmd = (
        f"{decompressor} {barcodes_fastq} | "
        f"awk 'NR%4==2 && length($0)=={barcode_len} {{print}}' | "
        f"{sort_opts} | "
        f"uniq -c | "
        f"awk '{{print $2\"\\t\"$1}}' "
        f"> {output_tsv}"
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"barcode counting failed for {barcodes_fastq}:\n{result.stderr}")


def count_barcodes_python(
    barcodes_fastq: str | Path,
    output_tsv: str | Path,
    barcode_len: int,
) -> None:
    """Count barcodes in pure Python (used in tests; prefer count_barcodes_shell in production)."""
    import gzip as _gzip

    barcodes_fastq = Path(barcodes_fastq)
    counts: Counter[str] = Counter()

    opener = _gzip.open if str(barcodes_fastq).endswith(".gz") else open
    with opener(barcodes_fastq, "rt") as f:
        while True:
            f.readline()  # @header
            seq = f.readline().rstrip("\n")
            f.readline()  # +
            f.readline()  # quality
            if not seq:
                break
            if len(seq) == barcode_len and set(seq) <= {"A", "C", "G", "T"}:
                counts[seq] += 1

    with open(output_tsv, "w") as out:
        for barcode, count in sorted(counts.items()):
            out.write(f"{barcode}\t{count}\n")
