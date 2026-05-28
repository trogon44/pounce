from __future__ import annotations

from pathlib import Path

import click
import pandas as pd


@click.group()
def shortread():
    """Short-read (Illumina) barcode counting and matrix assembly."""


@shortread.command()
@click.argument("plate_info_csv")
@click.argument("plasmid_fasta")
@click.argument("wt_orf_fasta")
@click.option("--orf-start", required=True, type=int)
@click.option("--orf-end", required=True, type=int)
@click.option("--out-dir", "-o", required=True)
@click.option("--error-rate", default=0.05, show_default=True)
@click.option("--threads", default=1, show_default=True)
@click.option("--adapter-flank", default=20, show_default=True)
def trim(plate_info_csv, plasmid_fasta, wt_orf_fasta, orf_start, orf_end,
         out_dir, error_rate, threads, adapter_flank):
    """Trim barcodes from R2 FASTQ files using adapter sequences from plasmid."""
    from pounce.reference import parse_reference
    from pounce.shortread.trim import trim_barcodes

    ref = parse_reference(plasmid_fasta, wt_orf_fasta, orf_start, orf_end, adapter_flank)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    plate = pd.read_csv(plate_info_csv)
    for _, row in plate.iterrows():
        sid = row["sample_id"]
        output = out / f"{sid}_barcodes.fastq.gz"
        click.echo(f"Trimming {sid}...")
        trim_barcodes(row["fastq_r2"], output, ref,
                      error_rate=error_rate, threads=threads)
    click.echo(f"Trimmed {len(plate)} samples → {out_dir}")


@shortread.command()
@click.argument("plate_info_csv")
@click.option("--barcode-len", default=30, show_default=True)
@click.option("--barcodes-dir", default=".", show_default=True,
              help="Directory containing {sample_id}_barcodes.fastq.gz files")
@click.option("--out-dir", "-o", required=True)
@click.option("--sort-buffer", default="4G", show_default=True)
@click.option("--sort-parallel", default=0, show_default=True)
def count(plate_info_csv, barcode_len, barcodes_dir, out_dir, sort_buffer, sort_parallel):
    """Count barcodes per sample using the awk/sort/uniq pipeline."""
    from pounce.shortread.count import count_barcodes_shell

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    barcodes_dir = Path(barcodes_dir)

    plate = pd.read_csv(plate_info_csv)
    for _, row in plate.iterrows():
        sid = row["sample_id"]
        fastq = barcodes_dir / f"{sid}_barcodes.fastq.gz"
        tsv = out / f"{sid}_barcode_counts.tsv"
        click.echo(f"Counting {sid}...")
        count_barcodes_shell(fastq, tsv, barcode_len,
                             sort_buffer=sort_buffer, sort_parallel=sort_parallel)
    click.echo(f"Counted {len(plate)} samples → {out_dir}")


@shortread.command()
@click.argument("plate_info_csv")
@click.option("--counts-dir", required=True,
              help="Directory containing {sample_id}_barcode_counts.tsv files")
@click.option("--output", "-o", required=True, help="Output barcode_counts.csv")
@click.option("--ns-min-count", default=1, show_default=True,
              help="Min count in any NS sample for barcode to be included")
def matrix(plate_info_csv, counts_dir, output, ns_min_count):
    """Build barcode × sample count matrix from per-sample TSVs."""
    from pounce.shortread.matrix import build_count_matrix

    plate = pd.read_csv(plate_info_csv)
    counts_dir = Path(counts_dir)
    plate["count_file"] = plate["sample_id"].apply(
        lambda sid: str(counts_dir / f"{sid}_barcode_counts.tsv")
    )

    mat = build_count_matrix(plate, ns_min_count=ns_min_count)
    mat.to_csv(output)
    click.echo(f"Matrix: {mat.shape[0]:,} barcodes × {mat.shape[1]} samples → {output}")


@shortread.command()
@click.argument("plate_info_csv")
@click.argument("plasmid_fasta")
@click.argument("wt_orf_fasta")
@click.option("--orf-start", required=True, type=int)
@click.option("--orf-end", required=True, type=int)
@click.option("--out-dir", "-o", required=True)
@click.option("--ns-min-count", default=1, show_default=True)
@click.option("--threads", default=1, show_default=True)
@click.pass_context
def run(ctx, plate_info_csv, plasmid_fasta, wt_orf_fasta, orf_start, orf_end,
        out_dir, ns_min_count, threads):
    """Run the full short-read pipeline: trim → count → matrix."""
    out = Path(out_dir)
    ctx.invoke(trim, plate_info_csv=plate_info_csv,
               plasmid_fasta=plasmid_fasta, wt_orf_fasta=wt_orf_fasta,
               orf_start=orf_start, orf_end=orf_end,
               out_dir=str(out / "barcodes"), threads=threads)
    ctx.invoke(count, plate_info_csv=plate_info_csv,
               barcodes_dir=str(out / "barcodes"),
               out_dir=str(out / "counts"))
    ctx.invoke(matrix, plate_info_csv=plate_info_csv,
               counts_dir=str(out / "counts"),
               output=str(out / "barcode_counts.csv"),
               ns_min_count=ns_min_count)
