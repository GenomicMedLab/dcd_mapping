"""Test ``transcripts`` module."""
from typing import Dict

import pytest
from cool_seq_tool.schemas import TranscriptPriority

from mavemap.resources import get_scoreset_records
from mavemap.schemas import AlignmentResult, ScoresetMetadata
from mavemap.transcripts import select_transcript


@pytest.mark.asyncio(scope="module")
async def test_tx_src(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)  # TODO real fixture
    alignment_result = align_result_fixture[urn]

    actual = await select_transcript(metadata, records, alignment_result)

    assert actual
    assert actual.np == "NP_938033.1"
    assert actual.start == 269
    assert actual.is_full_match is True
    assert actual.nm == "NM_198291.3"
    assert actual.transcript_mode == TranscriptPriority.MANE_SELECT


@pytest.mark.asyncio(scope="module")
async def test_tx_scn5a(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_transcript(metadata, records, alignment_result)

    assert actual
    assert actual.np == "NP_000326.2"
    assert actual.start == 1619
    assert actual.is_full_match is True
    assert actual.nm == "NM_000335.5"
    assert actual.transcript_mode == TranscriptPriority.MANE_SELECT


@pytest.mark.asyncio(scope="module")
async def test_tx_hbb(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_transcript(metadata, records, alignment_result)
    assert actual is None


@pytest.mark.asyncio(scope="module")
async def test_tx_raf(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000061-h-1"""
    urn = "urn:mavedb:00000061-h-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_transcript(metadata, records, alignment_result)
    assert actual
    assert actual.np == "NP_002871.1"
    assert actual.start == 51
    assert actual.is_full_match is True
    assert actual.nm == "NM_002880.4"
    assert actual.transcript_mode == TranscriptPriority.MANE_SELECT


@pytest.mark.asyncio(scope="module")
async def test_tx_brca(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for"""
