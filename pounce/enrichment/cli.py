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
@click.option("--output", "-o", required=True, help="Output barcode_enrichments.csv")
@click.option("--pseudocount", default=1.0, show_default=True)
@click.option("--ns-min", default=10, show_default=True)
def barcodes(barcode_counts_csv, plate_info_csv, timepoint, output, pseudocount, ns_min):
    """Compute per-barcode log2 enrichments for each replicate at a given timepoint."""
    from pounce.enrichment.compute import compute_replicate_enrichments

    counts = pd.read_csv(barcode_counts_csv, index_col="barcode")
    plate = pd.read_csv(plate_info_csv)

    result = compute_replicate_enrichments(
        counts, plate, timepoint=timepoint,
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
@click.option("--min-barcodes", default=3, show_default=True)
def variants(barcode_enrichments_csv, barcode_variants_csv, output, min_barcodes):
    """Aggregate per-barcode enrichments to per-AA-variant median enrichments."""
    from pounce.enrichment.variants import (
        load_barcode_variants,
        add_variant_column,
        aggregate_enrichments,
    )

    enrichments = pd.read_csv(barcode_enrichments_csv, index_col="barcode")
    bv = add_variant_column(load_barcode_variants(barcode_variants_csv))

    result = aggregate_enrichments(enrichments, bv, min_barcodes=min_barcodes)
    result.to_csv(output)
    click.echo(f"{len(result):,} variants → {output}")


@enrichment.command()
@click.argument("variant_enrichments_csv")
@click.option("--output", "-o", required=True, help="Output normalized_enrichments.csv")
@click.option("--min-barcodes", default=10, show_default=True,
              help="Min barcodes for WT/stop anchor variants")
def normalize(variant_enrichments_csv, output, min_barcodes):
    """Normalize variant enrichments to WT=1.0, single stops=0.0."""
    from pounce.enrichment.normalize import normalize_to_wt_stop

    result = pd.read_csv(variant_enrichments_csv, index_col="variant")
    result["enrich_normalized"] = normalize_to_wt_stop(result, min_barcodes=min_barcodes)
    result.to_csv(output)
    click.echo(f"Normalized {len(result):,} variants → {output}")
