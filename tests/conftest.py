"""Provide shared testing utilities."""
import asyncio
import json
from pathlib import Path

import pytest

from mavemap.schemas import ScoresetMetadata


@pytest.fixture(scope="session")
def event_loop(request):
    """Create an instance of the default event loop for each test case."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


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
