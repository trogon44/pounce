from __future__ import annotations

import csv
import gzip
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def _parse_nuc_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse 'A334G' → ('A', 334, 'G'). Position is 1-based."""
    return mut_str[0], int(mut_str[1:-1]), mut_str[-1]


def call_nt_mutations(consensus_seq: str, wt_seq: str) -> list[str]:
    """Return list of nucleotide mutations as 'REF{pos}ALT' strings (1-based).

    A consensus ``N`` (every read has a deletion or no call there) is reported
    as a mutation to ``N`` rather than ignored, so such barcodes are not
    mistaken for wild type.
    """
    if len(consensus_seq) != len(wt_seq):
        raise ValueError(
            f"Length mismatch: consensus={len(consensus_seq)}, WT={len(wt_seq)}"
        )
    return [
        f"{wt}{i + 1}{cons}"
        for i, (wt, cons) in enumerate(zip(wt_seq, consensus_seq))
        if cons != wt
    ]


def nuc_to_aa_mutations(nt_mutations: list[str], wt_orf_seq: str) -> list[str]:
    """Convert nucleotide mutation list to amino acid mutation list.

    Synonymous mutations produce no entries. Stop codons are represented as '*'.
    Multiple mutations in the same codon are resolved together. A codon that
    contains an ``N`` is reported as ``{wt_aa}{pos}-``. Entries are ordered by
    residue position.
    """
    if not nt_mutations:
        return []

    wt_bases = list(wt_orf_seq)
    mut_bases = list(wt_orf_seq)

    parsed = [_parse_nuc_mutation(m) for m in nt_mutations]
    for _, pos_1b, alt in parsed:
        mut_bases[pos_1b - 1] = alt

    affected_codons = sorted({(pos_1b - 1) // 3 for _, pos_1b, _ in parsed})

    aa_muts: list[str] = []
    for codon_idx in affected_codons:
        wt_codon = "".join(wt_bases[codon_idx * 3:(codon_idx + 1) * 3])
        mut_codon = "".join(mut_bases[codon_idx * 3:(codon_idx + 1) * 3])
        wt_aa = str(Seq(wt_codon).translate())
        mut_aa = "-" if "N" in mut_codon else str(Seq(mut_codon).translate())
        if wt_aa != mut_aa:
            aa_muts.append(f"{wt_aa}{codon_idx + 1}{mut_aa}")

    return aa_muts


def call_all_mutations(
    filtered_fasta: str | Path,
    wt_orf_fasta: str | Path,
    output_csv: str | Path,
) -> int:
    """Call nt and AA mutations for every barcode in the filtered consensus FASTA.

    Output CSV columns: barcode, nt_mutations, aa_mutations, read_count
    Mutations are semicolon-separated strings, empty string when WT.

    Returns number of barcodes written.
    """
    filtered_fasta = Path(filtered_fasta)
    wt_orf_fasta = Path(wt_orf_fasta)
    output_csv = Path(output_csv)

    wt_records = list(SeqIO.parse(wt_orf_fasta, "fasta"))
    if len(wt_records) != 1:
        raise ValueError(f"Expected one WT ORF record, found {len(wt_records)}")
    wt_seq = str(wt_records[0].seq).upper()

    n_written = 0
    opener = gzip.open if str(filtered_fasta).endswith(".gz") else open
    with opener(filtered_fasta, "rt") as fasta_in, open(output_csv, "w", newline="") as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(["barcode", "nt_mutations", "aa_mutations", "read_count"])

        header = None
        seq_buf: list[str] = []

        def _flush():
            nonlocal n_written
            if header is None:
                return
            # Header format: "{barcode} size={N}" or just "{barcode}"
            parts = header.split()
            barcode = parts[0]
            size_str = parts[1].split("=")[-1] if len(parts) > 1 else "0"
            read_count = int(size_str) if size_str.isdigit() else 0

            seq = "".join(seq_buf).upper()
            nt_muts = call_nt_mutations(seq, wt_seq)
            aa_muts = nuc_to_aa_mutations(nt_muts, wt_seq)
            writer.writerow([
                barcode,
                ";".join(nt_muts),
                ";".join(aa_muts),
                read_count,
            ])
            n_written += 1

        for line in fasta_in:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                header = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line)
        _flush()

    return n_written
