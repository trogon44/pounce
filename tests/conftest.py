"""Shared fixtures for pounce tests."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "fixtures" / "data"


def pytest_configure(config):
    if not DATA.exists() or not (DATA / "test_plasmid.fa").exists():
        import subprocess, sys
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "fixtures" / "generate.py")],
            check=True,
        )


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA


@pytest.fixture(scope="session")
def plasmid_fasta(data_dir) -> Path:
    return data_dir / "test_plasmid.fa"


@pytest.fixture(scope="session")
def wt_orf_fasta(data_dir) -> Path:
    return data_dir / "test_wt_orf.fa"


@pytest.fixture(scope="session")
def orf_coords(data_dir) -> tuple[int, int]:
    lines = (data_dir / "test_coords.txt").read_text().strip().splitlines()
    return int(lines[0]), int(lines[1])


@pytest.fixture(scope="session")
def ref_config(plasmid_fasta, wt_orf_fasta, orf_coords):
    from pounce.reference import parse_reference
    start, end = orf_coords
    return parse_reference(plasmid_fasta, wt_orf_fasta, start, end)


@pytest.fixture(scope="session")
def truth_df(data_dir):
    import pandas as pd
    return pd.read_csv(data_dir / "test_barcodes_truth.csv", index_col="barcode")


@pytest.fixture(scope="session")
def plate_info(data_dir):
    import pandas as pd
    return pd.read_csv(data_dir / "test_plate_info.csv")
