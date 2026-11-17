from __future__ import annotations

from pathlib import Path

import click


@click.group()
def longread():
    """Long-read (PacBio) barcode extraction and variant mapping."""


_orf_coverage_option = click.option(
    "--allow-partial-orf/--require-full-orf", default=True, show_default=True,
    help="Keep reads that span the barcode but cover only part of the ORF, "
         "calling uncovered bases as N (they add votes but never override covered bases)",
)


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
@click.option("--threads", default=4, show_default=True,
              help="Threads for BAM decompression and for aligning an unaligned input BAM")
@_orf_coverage_option
def extract(bam, plasmid_fasta, wt_orf_fasta, orf_start, orf_end,
            output, mapq_min, adapter_flank, sort_buffer, threads, allow_partial_orf):
    """Extract (barcode, orf_seq) pairs from BAM, sorted by barcode."""
    from pounce.reference import parse_reference
    from pounce.longread.extract import extract_reads_sorted

    ref = parse_reference(plasmid_fasta, wt_orf_fasta, orf_start, orf_end, adapter_flank)
    click.echo(f"Barcode region: {ref.barcode_start}–{ref.barcode_end} ({ref.barcode_len} bp)")
    extract_reads_sorted(
        bam, ref, output, mapq_min=mapq_min, sort_buffer=sort_buffer, threads=threads,
        allow_partial_orf=allow_partial_orf,
    )
    click.echo(f"Written: {output}")


@longread.command()
@click.argument("reads_tsv_gz")
@click.option("--output", "-o", required=True, help="Output consensus.fastq.gz")
@click.option("--threads", default=1, show_default=True, help="Processes for consensus calling")
def consensus(reads_tsv_gz, output, threads):
    """Call per-barcode consensus from sorted reads TSV."""
    from pounce.longread.consensus import call_consensus

    n = call_consensus(reads_tsv_gz, output, threads=threads)
    click.echo(f"Called consensus for {n:,} barcodes → {output}")


_FILTER_OPTIONS = [
    click.option("--min-reads", default=2, show_default=True,
                 help="Min PacBio reads per barcode"),
    click.option("--min-mean-quality", default=10.0, show_default=True,
                 help="Min mean consensus Phred quality"),
    click.option("--max-n", default=5, show_default=True,
                 help="Max N bases in a consensus"),
    click.option("--low-support-reads", default=2, show_default=True,
                 help="Barcodes with this many reads or fewer get the stricter thresholds below"),
    click.option("--low-support-min-quality", default=20.0, show_default=True),
    click.option("--low-support-max-n", default=0, show_default=True),
]


def _filter_options(func):
    for option in reversed(_FILTER_OPTIONS):
        func = option(func)
    return func


@longread.command("filter")
@click.argument("consensus_fastq_gz")
@click.option("--output", "-o", required=True, help="Output filtered.fasta.gz")
@_filter_options
def filter_cmd(consensus_fastq_gz, output, min_reads, min_mean_quality, max_n,
               low_support_reads, low_support_min_quality, low_support_max_n):
    """Filter consensus FASTQ by read count and sequence quality."""
    from pounce.longread.filter import filter_consensus

    total, passed = filter_consensus(
        consensus_fastq_gz, output, min_reads, min_mean_quality, max_n,
        low_support_reads, low_support_min_quality, low_support_max_n,
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
@click.option("--threads", default=4, show_default=True,
              help="Threads for BAM decompression and for aligning an unaligned input BAM")
@_orf_coverage_option
@_filter_options
@click.pass_context
def run(ctx, bam, plasmid_fasta, wt_orf_fasta, orf_start, orf_end, out_dir,
        mapq_min, adapter_flank, threads, allow_partial_orf, min_reads, min_mean_quality, max_n,
        low_support_reads, low_support_min_quality, low_support_max_n):
    """Run the full long-read pipeline: extract → consensus → filter → mutations."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    ctx.invoke(extract,
               bam=bam, plasmid_fasta=plasmid_fasta, wt_orf_fasta=wt_orf_fasta,
               orf_start=orf_start, orf_end=orf_end,
               output=str(out / "reads.tsv.gz"),
               mapq_min=mapq_min, adapter_flank=adapter_flank, sort_buffer="2G",
               threads=threads, allow_partial_orf=allow_partial_orf)

    ctx.invoke(consensus,
               reads_tsv_gz=str(out / "reads.tsv.gz"),
               output=str(out / "consensus.fastq.gz"),
               threads=threads)

    ctx.invoke(filter_cmd,
               consensus_fastq_gz=str(out / "consensus.fastq.gz"),
               output=str(out / "consensus_filtered.fasta.gz"),
               min_reads=min_reads, min_mean_quality=min_mean_quality, max_n=max_n,
               low_support_reads=low_support_reads,
               low_support_min_quality=low_support_min_quality,
               low_support_max_n=low_support_max_n)

    ctx.invoke(mutations,
               filtered_fasta_gz=str(out / "consensus_filtered.fasta.gz"),
               wt_orf_fasta=wt_orf_fasta,
               output=str(out / "barcode_variants.csv"))
