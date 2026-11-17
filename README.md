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

The input BAM can be the raw instrument HiFi BAM (it is aligned to `plasmid.fa` with
minimap2 automatically) or an already coordinate-sorted alignment to the plasmid.

```bash
# Run all steps at once
pounce longread run reads.bam plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 \
    --out-dir lr_out/ --threads 12

# Or step by step
pounce longread extract reads.bam plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 -o reads.tsv.gz
pounce longread consensus reads.tsv.gz -o consensus.fastq.gz --threads 12
pounce longread filter consensus.fastq.gz -o consensus_filtered.fasta.gz
pounce longread mutations consensus_filtered.fasta.gz wt_orf.fa -o barcode_variants.csv
```

`barcode_variants.csv` columns: `barcode`, `nt_mutations`, `aa_mutations`, `read_count`

Mutations are semicolon-separated strings in position order, e.g. `A334G;C501T` /
`K112E;R201W`. Empty string = no mutation. A base that every read lacks (a deletion, or
positions no read covers) is reported as `N` (`G534N`), and a codon containing an `N` as
`L178-`; such barcodes are classed as `frameshift` downstream rather than mistaken for WT.

How reads are used:

- A read must span the whole barcode with no deletion in it. By default it only has to
  *overlap* the ORF (`--allow-partial-orf`); bases it does not cover are `N`, which never
  override a covered base but do count against the consensus quality. About 12% of reads in
  the Dec 2025 library start partway into the ORF. Use `--require-full-orf` to drop them.
- Consensus is a per-position majority vote over A/C/G/T, quality is the Phred-scaled
  fraction of *all* reads supporting the call.
- The consensus filter keeps barcodes with at least 2 reads, mean quality ≥ 10 and at most
  5 `N`s; barcodes with 2 reads or fewer must additionally have quality ≥ 20 and no `N`
  (`--min-reads`, `--min-mean-quality`, `--max-n`, `--low-support-*`).

### Short-read pipeline (Illumina → counts)

```bash
pounce shortread run plate_info.csv plasmid.fa wt_orf.fa \
    --orf-start 2931 --orf-end 3945 \
    --out-dir sr_out/ --ns-min-count 6 --threads 12
```

Produces `sr_out/barcode_counts.csv`: barcodes × samples count matrix. Barcodes are kept if
they reach `--ns-min-count` reads in at least one non-selective sample.

### Enrichment

```bash
# Per-barcode enrichment (works without variant mapping)
pounce enrichment barcodes sr_out/barcode_counts.csv plate_info.csv \
    --timepoint T2 --ns-min 10 -o barcode_enrichments.csv

# Aggregate to per-variant statistics (requires barcode_variants.csv)
pounce enrichment variants barcode_enrichments.csv barcode_variants.csv \
    -o variant_enrichments.csv

# Add enrich_norm: stop variants = 0, WT-synonymous variants = 1
pounce enrichment normalize variant_enrichments.csv -o normalized_enrichments.csv
```

Each `S` sample is paired with the one `NS` sample that has the same `timepoint`,
`replicate` and, if `plate_info.csv` has a `harvest` column, the same `harvest`. NS samples
with no matching S (for example 00 h or 24 h controls) are simply not used. pounce stops with
an error instead of guessing when a group has more than one S or NS sample, or when several
harvests have complete pairs and `--harvest` was not given:

```
sample_id, fastq_r2, timepoint, replicate, condition, harvest
SP11_72h, ..., T0, A, S,  72h
SP12_72h, ..., T0, A, NS, 72h
SP01_00h, ..., T0, A, NS, 00h      # ignored: no S at 00h
```

```bash
pounce enrichment barcodes counts.csv plate_info.csv --timepoint T0 --harvest 72h -o T0_72h.csv
```

Output columns are always `enrich_<replicate>`, whichever harvest was chosen.

`variant_enrichments.csv` has one row per variant with `variant_type`, `aa_mutcount` and,
for each replicate, `count_*`, `median_*`, `mean_*`, `std_*`, plus `mean_enrich`. Variants
are keyed and typed as follows:

| `variant_type` | Rule | Key |
|---|---|---|
| `missense` | amino-acid changes, no stop or `-` | `K2Q;W300R` |
| `stop` | contains a `*` | `K2*;E234K` |
| `wt_syn` | nucleotide changes, all synonymous | `WTsyn_A6G;C9T` |
| `frameshift` | a codon with a deleted/uncalled base | nucleotide string |
| `wt` | no mutations | `WT` |

No barcode-count threshold is applied when the table is built. Filter on `count_*`
afterwards (for example `count_A >= 10 and count_B >= 10`).

### Key parameters

| Flag | Default | Notes |
|------|---------|-------|
| `--pseudocount` | `1` | Raw-count pseudocount (enrichment barcodes) |
| `--ns-min` | `10` | Min NS reads per barcode, applied per replicate on its own NS column |
| `--ns-min-count` | `1` | Short-read matrix: min count in any NS sample (`6` reproduces the ePCR library) |
| `--min-barcodes` | `10` | `enrichment normalize` only: min barcodes (every replicate) for an anchor variant |
| `--anchor-stat` | `mean` | `enrichment normalize`: mean or median across anchor variants |
| `--min-reads` | `2` | Min PacBio reads per barcode for consensus |

## Enrichment formula

```
E = log2((S + pc) / (NS + pc)) + log2(NS_total / S_total)
```

`pc` is a raw read count (default 1), added to numerator and denominator. `S_total` and
`NS_total` are the column sums of the count matrix passed in, so the library-size term makes
E = 0 for a barcode at identical frequency in both pools.

Each replicate is summarised independently: the variant value is the **median** E over the
barcodes of that variant (limits the effect of outlier barcodes), and `mean_enrich` is the
**mean of the replicate medians**. `enrich_norm = (mean_enrich - dead) / (wt - dead)`, where
`dead` and `wt` are the mean of `mean_enrich` over the `stop` and `wt_syn` variants with
enough barcodes in every replicate.

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
