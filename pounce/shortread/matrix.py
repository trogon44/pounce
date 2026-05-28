from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_count_tsv(path: str | Path) -> pd.Series:
    """Load a barcode_counts.tsv into a Series indexed by barcode."""
    df = pd.read_csv(path, sep="\t", names=["barcode", "count"], header=None)
    return df.set_index("barcode")["count"]


def build_count_matrix(
    plate_info: pd.DataFrame,
    ns_min_count: int = 1,
) -> pd.DataFrame:
    """Aggregate per-sample TSVs into a barcode × sample count matrix.

    Parameters
    ----------
    plate_info:
        DataFrame with at least columns: sample_id, count_file, condition.
        ``count_file`` should point to the *_barcode_counts.tsv for each sample.
    ns_min_count:
        Minimum count in any NS (non-selective) sample for a barcode to be
        included. Barcodes never reaching this threshold in any NS sample are dropped.

    Returns
    -------
    DataFrame with barcodes as index and sample_ids as columns.
    """
    per_sample: dict[str, pd.Series] = {}
    for _, row in plate_info.iterrows():
        per_sample[row["sample_id"]] = load_count_tsv(row["count_file"])

    ns_ids = plate_info.loc[plate_info["condition"] == "NS", "sample_id"].tolist()
    if not ns_ids:
        raise ValueError("plate_info contains no NS (non-selective) samples")

    good_barcodes: set[str] = set()
    for sid in ns_ids:
        s = per_sample[sid]
        good_barcodes |= set(s[s >= ns_min_count].index)

    if not good_barcodes:
        raise ValueError(
            f"No barcodes meet ns_min_count={ns_min_count} in any NS sample"
        )

    index = sorted(good_barcodes)
    matrix = pd.DataFrame(0, index=index, columns=plate_info["sample_id"].tolist())
    matrix.index.name = "barcode"

    for sid, series in per_sample.items():
        shared = series.index.intersection(matrix.index)
        matrix.loc[shared, sid] = series.loc[shared]

    return matrix
