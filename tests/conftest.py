"""Provide shared testing utilities."""
import json
import os
from pathlib import Path

import pytest

from dcd_mapping.schemas import AlignmentResult, ScoresetMetadata, TxSelectResult

FIXTURE_DATA_DIR = Path(__file__).parents[0].resolve() / "fixtures"


def pytest_sessionstart(session) -> None:
    """Initialize testing environment."""
    os.environ["MAVEDB_STORAGE_DIR"] = str(FIXTURE_DATA_DIR.absolute())


@pytest.fixture(scope="session")
def fixture_data_dir():
    """Provide test data directory."""
    return FIXTURE_DATA_DIR


@pytest.fixture(scope="module")
def scoreset_metadata_fixture(fixture_data_dir: Path):
    """Provide scoreset metadata fixtures."""
    fixture_file = fixture_data_dir / "scoreset_metadata.json"
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for d in data["scoreset_metadata"]:
        formatted_data = ScoresetMetadata(**d)
        results[formatted_data.urn] = formatted_data
    return results


@pytest.fixture(scope="session")
def align_result_fixture(fixture_data_dir: Path):
    """Provide fixtures for alignment results."""
    fixture_file = fixture_data_dir / "align_result.json"
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = AlignmentResult(**result)
        results[urn] = formatted_result
    return results


@pytest.fixture(scope="session")
def transcript_results_fixture(fixture_data_dir: Path):
    """Provide fixtures for transcript selection results."""
    fixture_file = fixture_data_dir / "transcript_result.json"
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for urn, result in data.items():
        formatted_result = TxSelectResult(**result)
        results[urn] = formatted_result
    return results
