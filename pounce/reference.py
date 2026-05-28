from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True)
class ReferenceConfig:
    plasmid_path: Path
    wt_orf_path: Path
    barcode_start: int    # 0-based inclusive
    barcode_end: int      # 0-based exclusive
    barcode_len: int
    orf_start: int        # 0-based inclusive, within plasmid
    orf_end: int          # 0-based exclusive, within plasmid
    adapter_5p: str       # plasmid sequence upstream of barcode, used by cutadapt
    adapter_3p: str       # plasmid sequence downstream of barcode, used by cutadapt

    @property
    def orf_len(self) -> int:
        return self.orf_end - self.orf_start


def parse_reference(
    plasmid_fasta: str | Path,
    wt_orf_fasta: str | Path,
    orf_start: int,
    orf_end: int,
    adapter_flank: int = 20,
) -> ReferenceConfig:
    """Build a ReferenceConfig from a plasmid FASTA with Ns marking the barcode.

    Parameters
    ----------
    plasmid_fasta:
        Path to FASTA with exactly one contiguous run of N's marking the barcode region.
    wt_orf_fasta:
        Path to FASTA of the wild-type ORF sequence (used for mutation calling).
    orf_start:
        0-based start of the ORF within the plasmid sequence.
    orf_end:
        0-based exclusive end of the ORF within the plasmid sequence.
    adapter_flank:
        Number of plasmid bases flanking the barcode to use as cutadapt adapter sequences.
    """
    plasmid_path = Path(plasmid_fasta)
    wt_orf_path = Path(wt_orf_fasta)

    records = list(SeqIO.parse(plasmid_path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Plasmid FASTA must contain exactly one record, found {len(records)}")
    plasmid_seq = str(records[0].seq).upper()

    n_runs = list(re.finditer(r"N+", plasmid_seq))
    if len(n_runs) != 1:
        raise ValueError(
            f"Expected exactly one N-run in plasmid FASTA to mark the barcode, "
            f"found {len(n_runs)}"
        )

    m = n_runs[0]
    barcode_start = m.start()
    barcode_end = m.end()

    adapter_5p = plasmid_seq[max(0, barcode_start - adapter_flank):barcode_start]
    adapter_3p = plasmid_seq[barcode_end:barcode_end + adapter_flank]

    if len(adapter_5p) < adapter_flank or len(adapter_3p) < adapter_flank:
        raise ValueError(
            f"Barcode is too close to a plasmid end to extract {adapter_flank}-bp adapters; "
            "reduce --adapter-flank or check your plasmid FASTA"
        )

    wt_records = list(SeqIO.parse(wt_orf_path, "fasta"))
    if len(wt_records) != 1:
        raise ValueError(f"WT ORF FASTA must contain exactly one record, found {len(wt_records)}")
    wt_orf_len = len(wt_records[0].seq)
    if wt_orf_len != orf_end - orf_start:
        raise ValueError(
            f"WT ORF FASTA length ({wt_orf_len}) does not match orf_end - orf_start "
            f"({orf_end - orf_start})"
        )

    return ReferenceConfig(
        plasmid_path=plasmid_path,
        wt_orf_path=wt_orf_path,
        barcode_start=barcode_start,
        barcode_end=barcode_end,
        barcode_len=barcode_end - barcode_start,
        orf_start=orf_start,
        orf_end=orf_end,
        adapter_5p=adapter_5p,
        adapter_3p=adapter_3p,
    )
