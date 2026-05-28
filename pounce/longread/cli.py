from __future__ import annotations

from pathlib import Path

import click


@click.group()
def longread():
    """Long-read (PacBio) barcode extraction and variant mapping."""


@longread.command()
@click.argument("bam")
@click.argument("plasmid_fasta")
@click.argument("wt_orf_fasta")
@click.option("--orf-start", required=True, type=int, help="0-based ORF start in plasmid")
@click.option("--orf-end", required=True, type=int, help="0-based exclusive ORF end in plasmid")
@click.option("--output", "-o", required=True, help="Output reads.tsv.gz")
@click.option("--mapq-min", default=0, show_default=True)
@click.option("--adapter-flank", default=20, show_default=True)
@click.option("--sort-buffer", default="2G", show_default=True)
def extract(bam, plasmid_fasta, wt_orf_fasta, orf_start, orf_end,
            output, mapq_min, adapter_flank, sort_buffer):
    """Extract (barcode, orf_seq) pairs from BAM, sorted by barcode."""
    from pounce.reference import parse_reference
    from pounce.longread.extract import extract_reads_sorted

    ref = parse_reference(plasmid_fasta, wt_orf_fasta, orf_start, orf_end, adapter_flank)
    click.echo(f"Barcode region: {ref.barcode_start}–{ref.barcode_end} ({ref.barcode_len} bp)")
    extract_reads_sorted(bam, ref, output, mapq_min=mapq_min, sort_buffer=sort_buffer)
    click.echo(f"Written: {output}")


@longread.command()
@click.argument("reads_tsv_gz")
@click.option("--output", "-o", required=True, help="Output consensus.fastq.gz")
def consensus(reads_tsv_gz, output):
    """Call per-barcode consensus from sorted reads TSV."""
    from pounce.longread.consensus import call_consensus

    n = call_consensus(reads_tsv_gz, output)
    click.echo(f"Called consensus for {n:,} barcodes → {output}")


@longread.command("filter")
@click.argument("consensus_fastq_gz")
@click.option("--output", "-o", required=True, help="Output filtered.fasta.gz")
@click.option("--min-reads", default=3, show_default=True)
@click.option("--max-n-frac", default=0.01, show_default=True)
@click.option("--min-mean-quality", default=20.0, show_default=True)
def filter_cmd(consensus_fastq_gz, output, min_reads, max_n_frac, min_mean_quality):
    """Filter consensus FASTQ by read count and sequence quality."""
    from pounce.longread.filter import filter_consensus

    total, passed = filter_consensus(
        consensus_fastq_gz, output, min_reads, max_n_frac, min_mean_quality
    )
    click.echo(f"Passed: {passed:,} / {total:,} ({100 * passed / max(total, 1):.1f}%)")


@longread.command()
@click.argument("filtered_fasta_gz")
@click.argument("wt_orf_fasta")
@click.option("--output", "-o", required=True, help="Output barcode_variants.csv")
def mutations(filtered_fasta_gz, wt_orf_fasta, output):
    """Call nucleotide and amino-acid mutations for each consensus sequence."""
    from pounce.longread.mutations import call_all_mutations

    n = call_all_mutations(filtered_fasta_gz, wt_orf_fasta, output)
    click.echo(f"Called mutations for {n:,} barcodes → {output}")


@longread.command()
@click.argument("bam")
@click.argument("plasmid_fasta")
@click.argument("wt_orf_fasta")
@click.option("--orf-start", required=True, type=int)
@click.option("--orf-end", required=True, type=int)
@click.option("--out-dir", "-o", required=True, help="Output directory")
@click.option("--mapq-min", default=0, show_default=True)
@click.option("--adapter-flank", default=20, show_default=True)
@click.option("--min-reads", default=3, show_default=True)
@click.option("--max-n-frac", default=0.01, show_default=True)
@click.option("--min-mean-quality", default=20.0, show_default=True)
@click.pass_context
def run(ctx, bam, plasmid_fasta, wt_orf_fasta, orf_start, orf_end, out_dir,
        mapq_min, adapter_flank, min_reads, max_n_frac, min_mean_quality):
    """Run the full long-read pipeline: extract → consensus → filter → mutations."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    ctx.invoke(extract,
               bam=bam, plasmid_fasta=plasmid_fasta, wt_orf_fasta=wt_orf_fasta,
               orf_start=orf_start, orf_end=orf_end,
               output=str(out / "reads.tsv.gz"),
               mapq_min=mapq_min, adapter_flank=adapter_flank, sort_buffer="2G")

    ctx.invoke(consensus,
               reads_tsv_gz=str(out / "reads.tsv.gz"),
               output=str(out / "consensus.fastq.gz"))

    ctx.invoke(filter_cmd,
               consensus_fastq_gz=str(out / "consensus.fastq.gz"),
               output=str(out / "consensus_filtered.fasta.gz"),
               min_reads=min_reads, max_n_frac=max_n_frac,
               min_mean_quality=min_mean_quality)

    ctx.invoke(mutations,
               filtered_fasta_gz=str(out / "consensus_filtered.fasta.gz"),
               wt_orf_fasta=wt_orf_fasta,
               output=str(out / "barcode_variants.csv"))
