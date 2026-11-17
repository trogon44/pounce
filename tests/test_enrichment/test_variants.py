"""Tests for pounce.enrichment.variants."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pounce.enrichment.variants import aggregate_enrichments, classify_barcodes


def _barcode_variants(rows: dict[str, tuple[str, str]]) -> pd.DataFrame:
    df = pd.DataFrame.from_dict(rows, orient="index", columns=["nt_mutations", "aa_mutations"])
    df.index.name = "barcode"
    return df


def test_classification_rules():
    bv = _barcode_variants({
        "b_wt": ("", ""),
        "b_syn": ("A6G;C9T", ""),
        "b_mis": ("A4C", "K2Q"),
        "b_mis2": ("A400C;A4C", "K134Q;K2Q"),
        "b_stop": ("A4T;G700A", "K2*;E234K"),
        "b_fs": ("A4N;A20C", "K2-;E7A"),
        "b_fs_syn": ("A4N", "K2-"),
    })
    out = classify_barcodes(bv)
    assert out["variant_type"].to_dict() == {
        "b_wt": "wt", "b_syn": "wt_syn", "b_mis": "missense", "b_mis2": "missense",
        "b_stop": "stop", "b_fs": "frameshift", "b_fs_syn": "frameshift",
    }
    # missense/stop keyed by AA string in residue order; syn/frameshift keyed by nt string
    assert out["variant"].to_dict() == {
        "b_wt": "WT", "b_syn": "WTsyn_A6G;C9T", "b_mis": "K2Q", "b_mis2": "K2Q;K134Q",
        "b_stop": "K2*;E234K", "b_fs": "A4N;A20C", "b_fs_syn": "A4N",
    }
    assert out["aa_mutcount"].to_dict() == {
        "b_wt": 0, "b_syn": 0, "b_mis": 1, "b_mis2": 2, "b_stop": 2, "b_fs": 0, "b_fs_syn": 0,
    }
    assert out.loc["b_syn", "nt_mutcount"] == 2


def test_classification_sorts_by_position_not_alphabetically():
    bv = _barcode_variants({"b": ("G900A;A4C", "W300R;K2Q")})
    out = classify_barcodes(bv)
    assert out.loc["b", "variant"] == "K2Q;W300R"
    assert out.loc["b", "nt_mutations"] == "A4C;G900A"


def _enrichments() -> pd.DataFrame:
    # Five barcodes of variant K2Q, four with a replicate-A value, three with a replicate-B value
    return pd.DataFrame(
        {
            "enrich_A": [1.0, 2.0, 10.0, 3.0, np.nan],
            "enrich_B": [0.0, np.nan, 4.0, np.nan, 2.0],
        },
        index=[f"b{i}" for i in range(5)],
    )


def test_aggregate_known_values():
    bv = _barcode_variants({f"b{i}": ("A4C", "K2Q") for i in range(5)})
    table = aggregate_enrichments(_enrichments(), bv)

    row = table.loc["K2Q"]
    assert row["variant_type"] == "missense" and row["aa_mutcount"] == 1
    # independent per-replicate filtering: counts differ
    assert (row["count_A"], row["count_B"]) == (4, 3)
    assert row["median_A"] == pytest.approx(2.5)        # median of 1,2,3,10 resists the outlier
    assert row["mean_A"] == pytest.approx(4.0)
    assert row["std_A"] == pytest.approx(np.std([1, 2, 3, 10], ddof=1))
    assert row["median_B"] == pytest.approx(2.0)
    assert row["mean_enrich"] == pytest.approx((2.5 + 2.0) / 2)  # mean of replicate medians


def test_aggregate_applies_no_barcode_threshold():
    bv = _barcode_variants({f"b{i}": ("A4C", "K2Q") for i in range(4)} | {"b4": ("A7C", "K3Q")})
    table = aggregate_enrichments(_enrichments(), bv)
    # K3Q is supported by a single barcode, still present
    assert "K3Q" in table.index
    assert table.loc["K3Q", "count_A"] == 0 and np.isnan(table.loc["K3Q", "median_A"])
    assert table.loc["K3Q", "count_B"] == 1


def test_aggregate_column_order_and_syn_keys():
    bv = _barcode_variants({f"b{i}": ("A6G", "") for i in range(5)})
    table = aggregate_enrichments(_enrichments(), bv)
    assert list(table.index) == ["WTsyn_A6G"]
    assert list(table.columns) == [
        "variant_type", "aa_mutcount",
        "count_A", "median_A", "mean_A", "std_A",
        "count_B", "median_B", "mean_B", "std_B",
        "mean_enrich",
    ]
