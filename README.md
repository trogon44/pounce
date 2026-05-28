# pounce

Pooled barcode screen analysis for barcoded ePCR selection experiments.

Given a library of randomly mutagenized variants each tagged with a unique DNA barcode,
`pounce` handles three independent analysis stages:

1. **`longread`** — PacBio HiFi reads → barcode-to-variant mapping
2. **`shortread`** — Illumina reads → per-barcode count matrix
3. **`enrichment`** — count matrix → selective vs. non-selective enrichment scores

The stages are independent CLI commands. It is common to have Illumina data before PacBio
data, or to analyze barcode-level enrichments without variant mapping.

## Installation

```bash
pip install -e .
```

Requires `cutadapt` and GNU coreutils (`sort`, `uniq`) on `PATH` for the short-read pipeline.

## Inputs

Every run needs three reference files:

| File | Description |
|------|-------------|
| `plasmid.fa` | Full plasmid FASTA with **a single run of N's** marking the barcode position |
| `wt_orf.fa` | Wild-type ORF sequence (used for mutation calling) |
| `plate_info.csv` | Sample metadata (see schema below) |

The barcode coordinates and cutadapt adapter sequences are derived automatically from the
plasmid FASTA — you do not need to specify them manually.

**plate_info.csv schema:**

```
sample_id, fastq_r2, timepoint, replicate, condition
SP01, /path/to/SP01_R2.fastq.gz, T0, A, NS
SP02, /path/to/SP02_R2.fastq.gz, T0, B, NS
SP03, /path/to/SP03_R2.fastq.gz, T1, A, S
SP04, /path/to/SP04_R2.fastq.gz, T1, B, S
```

`condition` must be `S` (selective) or `NS` (non-selective).

## Usage

### Long-read pipeline (PacBio → variants)

```bash
# Run all steps at once
pounce longread run reads.bam plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 \
    --out-dir lr_out/

# Or step by step
pounce longread extract reads.bam plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 -o reads.tsv.gz
pounce longread consensus reads.tsv.gz -o consensus.fastq.gz
pounce longread filter consensus.fastq.gz -o consensus_filtered.fasta.gz
pounce longread mutations consensus_filtered.fasta.gz wt_orf.fa -o barcode_variants.csv
```

`barcode_variants.csv` columns: `barcode`, `nt_mutations`, `aa_mutations`, `read_count`

Mutations are semicolon-separated strings, e.g. `A334G;C501T` / `K112E;R201W`.
Empty string = wild type.

### Short-read pipeline (Illumina → counts)

```bash
pounce shortread run plate_info.csv plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 \
    --out-dir sr_out/
```

Produces `sr_out/barcode_counts.csv`: barcodes × samples count matrix.

### Enrichment

```bash
# Per-barcode enrichment (works without variant mapping)
pounce enrichment barcodes sr_out/barcode_counts.csv plate_info.csv \
    --timepoint T1 -o barcode_enrichments.csv

# Aggregate to per-AA-variant enrichments (requires barcode_variants.csv)
pounce enrichment variants barcode_enrichments.csv barcode_variants.csv \
    -o variant_enrichments.csv

# Normalize: WT synonymous = 1.0, single stop codons = 0.0
pounce enrichment normalize variant_enrichments.csv \
    -o normalized_enrichments.csv
```

### Key parameters

| Flag | Default | Notes |
|------|---------|-------|
| `--pseudocount` | `1` | Raw-count pseudocount (enrichment barcodes) |
| `--ns-min` | `10` | Min NS reads per barcode |
| `--min-barcodes` | `3` | Min barcodes per variant (applied post-hoc) |
| `--min-reads` | `3` | Min PacBio reads per barcode for consensus |

## Enrichment formula

```
E = log2((S + pc) / (NS + pc)) + log2(NS_total / S_total)
```

where `pc = pseudocount / library_size` (i.e., `pseudocount` is in raw-count units).
The library-size correction ensures E = 0 for a barcode at identical frequency in both pools.

Per-variant enrichment is the **median** across all barcodes mapping to that variant,
computed independently per replicate.

## Development

```bash
pip install -e ".[dev]"

# Generate synthetic test fixtures (done automatically by pytest if missing)
python tests/fixtures/generate.py

pytest tests/ -v
```

Test fixtures are fully synthetic: random ORF, 100 random barcodes, random SNP mutations,
and synthetic Illumina/PacBio-like reads with known ground-truth enrichments.
No real experimental data is required to run the test suite.
