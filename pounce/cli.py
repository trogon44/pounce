import click

from pounce.longread.cli import longread
from pounce.shortread.cli import shortread
from pounce.enrichment.cli import enrichment


@click.group()
@click.version_option()
def cli():
    """pounce: pooled barcode screen analysis.

    Three independent pipelines:\n
      longread   — PacBio barcode extraction and variant mapping\n
      shortread  — Illumina barcode counting and count matrix assembly\n
      enrichment — Selective vs. non-selective enrichment scoring
    """


cli.add_command(longread)
cli.add_command(shortread)
cli.add_command(enrichment)
