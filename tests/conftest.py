"""Provide shared testing utilities."""
import json
from pathlib import Path

import pytest

from mavemap.schemas import ScoresetMetadata


@pytest.fixture(scope="module")
def scoreset_metadata_fixture():
    """Provide scoreset metadata fixtures."""
    fixture_file = (
        Path(__file__).parents[0].resolve() / "fixtures" / "scoreset_metadata.json"
    )
    with open(fixture_file, "r") as f:
        data = json.load(f)
    results = {}
    for d in data["scoreset_metadata"]:
        formatted_data = ScoresetMetadata(**d)
        results[formatted_data.urn] = formatted_data
    return results
