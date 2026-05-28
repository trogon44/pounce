from __future__ import annotations

import subprocess
from pathlib import Path

from pounce.reference import ReferenceConfig


def trim_barcodes(
    fastq_r2: str | Path,
    output_path: str | Path,
    ref_config: ReferenceConfig,
    error_rate: float = 0.05,
    min_overlap: int = 10,
    threads: int = 1,
) -> None:
    """Run cutadapt to extract barcodes from R2 reads using linked adapters.

    The adapter sequences are derived automatically from the plasmid FASTA
    flanking the N-run barcode region (stored in ref_config).
    """
    adapter = (
        f"{ref_config.adapter_5p};min_overlap={min_overlap}"
        f"...{ref_config.adapter_3p};min_overlap={min_overlap}"
    )
    barcode_len = ref_config.barcode_len

    cmd = [
        "cutadapt",
        "-e", str(error_rate),
        "-a", adapter,
        "--discard-untrimmed",
        "-m", str(barcode_len - 1),
        "-M", str(barcode_len + 1),
        "-j", str(threads),
        "-o", str(output_path),
        str(fastq_r2),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"cutadapt failed for {fastq_r2}:\n{result.stderr}"
        )
