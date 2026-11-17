from __future__ import annotations

import click
import pandas as pd


@click.group()
def enrichment():
    """Compute enrichment scores from barcode count matrices."""


@enrichment.command()
@click.argument("barcode_counts_csv")
@click.argument("plate_info_csv")
@click.option("--timepoint", required=True, help="Timepoint to use (e.g. T2)")
@click.option("--harvest", default=None,
              help="Harvest to use (needs a 'harvest' column in plate_info); required "
                   "when several harvests have both an S and an NS sample")
@click.option("--output", "-o", required=True, help="Output barcode_enrichments.csv")
@click.option("--pseudocount", default=1.0, show_default=True)
@click.option("--ns-min", default=10, show_default=True)
def barcodes(barcode_counts_csv, plate_info_csv, timepoint, harvest, output, pseudocount, ns_min):
    """Compute per-barcode log2 enrichments for each replicate at a given timepoint.

    Each S sample is paired with the NS sample of the same timepoint, replicate
    and harvest. Ambiguous pairings are an error, never a guess.
    """
    from pounce.enrichment.compute import compute_replicate_enrichments

    counts = pd.read_csv(barcode_counts_csv, index_col="barcode")
    plate = pd.read_csv(plate_info_csv)

    result = compute_replicate_enrichments(
        counts, plate, timepoint=timepoint, harvest=harvest,
        pseudocount=pseudocount, ns_min=ns_min,
    )
    result.to_csv(output)
    click.echo(
        f"{len(result):,} barcodes × {result.shape[1]} replicates → {output}"
    )


@enrichment.command()
@click.argument("barcode_enrichments_csv")
@click.argument("barcode_variants_csv")
@click.option("--output", "-o", required=True, help="Output variant_enrichments.csv")
def variants(barcode_enrichments_csv, barcode_variants_csv, output):
    """Aggregate per-barcode enrichments to per-variant statistics.

    Within each replicate the barcodes of a variant are summarised by their
    median (also mean and std); ``mean_enrich`` is the mean of the replicate
    medians. No barcode-count threshold is applied: filter the ``count_*``
    columns of the output as needed.
    """
    from pounce.enrichment.variants import (
        load_barcode_variants,
        aggregate_enrichments,
    )

    enrichments = pd.read_csv(barcode_enrichments_csv, index_col="barcode")
    bv = load_barcode_variants(barcode_variants_csv)

    result = aggregate_enrichments(enrichments, bv)
    result.to_csv(output)
    click.echo(f"{len(result):,} variants → {output}")
    click.echo(result["variant_type"].value_counts().to_string())


@enrichment.command()
@click.argument("variant_enrichments_csv")
@click.option("--output", "-o", required=True, help="Output with an added enrich_norm column")
@click.option("--min-barcodes", default=10, show_default=True,
              help="Min barcodes, in every replicate, for a variant to serve as an anchor")
@click.option("--anchor-stat", type=click.Choice(["mean", "median"]), default="mean",
              show_default=True, help="How the dead and WT levels are taken across anchor variants")
def normalize(variant_enrichments_csv, output, min_barcodes, anchor_stat):
    """Add enrich_norm: stop variants = 0, WT-synonymous variants = 1."""
    from pounce.enrichment.normalize import normalize_to_anchors

    result = pd.read_csv(variant_enrichments_csv, index_col="variant")
    result["enrich_norm"] = normalize_to_anchors(
        result, min_barcodes=min_barcodes, anchor_stat=anchor_stat
    )
    result.to_csv(output)
    click.echo(f"Normalized {len(result):,} variants → {output}")
