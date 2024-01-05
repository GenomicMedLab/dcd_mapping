"""Test ``transcripts`` module."""
from typing import Dict

import pytest
from cool_seq_tool.schemas import Strand, TranscriptPriority

from mavemap.resources import get_scoreset_records
from mavemap.schemas import AlignmentResult, ScoresetMetadata, SequenceRange
from mavemap.transcripts import select_reference


@pytest.fixture(scope="session")
def align_result_fixture():
    """Provide fixtures for alignment results."""
    return {
        "urn:mavedb:00000041-a-1": AlignmentResult(
            chrom="chr20",
            strand=Strand.POSITIVE,
            coverage=100.0,
            ident_pct=0,  # TODO
            query_range=SequenceRange(start=0, end=750),
            query_subranges=[
                # TODO this looks wrong
                SequenceRange(start=0, end=52),
                SequenceRange(start=0, end=52),
                SequenceRange(start=0, end=52),
                SequenceRange(start=0, end=52),
                SequenceRange(start=0, end=52),
            ],
            hit_range=SequenceRange(start=37397802, end=37403325),
            hit_subranges=[
                # TODO this looks wrong
                SequenceRange(start=37397802, end=37397854),
                SequenceRange(start=37397802, end=37397854),
                SequenceRange(start=37397802, end=37397854),
                SequenceRange(start=37397802, end=37397854),
                SequenceRange(start=37397802, end=37397854),
            ],
        ),
        "urn:mavedb:00000098-a-1": AlignmentResult(
            chrom="chr3",
            strand=Strand.NEGATIVE,
            coverage=100.0,
            ident_pct=100.0,
            query_range=SequenceRange(start=0, end=12),
            query_subranges=[SequenceRange(start=0, end=12)],
            hit_range=SequenceRange(start=38551475, end=38551511),
            hit_subranges=[SequenceRange(start=38551475, end=38551511)],
        ),
        "urn:mavedb:00000018-a-1": AlignmentResult(
            chrom="chr11",
            strand=Strand.POSITIVE,
            coverage=100.0,
            ident_pct=100.0,
            query_range=SequenceRange(start=0, end=187),
            query_subranges=[SequenceRange(start=0, end=187)],
            hit_range=SequenceRange(start=5227021, end=5227208),
            hit_subranges=[SequenceRange(start=5227021, end=5227208)],
        ),
        "urn:mavedb:00000113-a-2": AlignmentResult(
            chrom="chr21",
            strand=Strand.NEGATIVE,
            coverage=100.0,
            ident_pct=100.0,
            query_range=SequenceRange(start=0, end=42),
            query_subranges=[
                SequenceRange(start=0, end=17),
                SequenceRange(start=17, end=42),
            ],
            hit_range=SequenceRange(start=25891793, end=25897623),
            hit_subranges=[
                SequenceRange(start=25897572, end=25897623),
                SequenceRange(start=25891793, end=25891868),
            ],
        ),
        "urn:mavedb:00000061-h-1": AlignmentResult(
            chrom="chr3",
            strand=Strand.NEGATIVE,
            coverage=100.0,
            ident_pct=100.0,
            query_range=SequenceRange(start=0, end=117),
            query_subranges=[
                SequenceRange(start=54, end=117),
                SequenceRange(start=0, end=54),
            ],
            hit_range=SequenceRange(start=12611999, end=12618568),
            hit_subranges=[
                SequenceRange(start=1261199, end=12612062),
                SequenceRange(start=12618514, end=12618568),
            ],
        ),
    }


@pytest.mark.asyncio
async def test_tx_src(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000041-a-1"""
    urn = "urn:mavedb:00000041-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)  # TODO real fixture
    alignment_result = align_result_fixture[urn]

    actual = await select_reference(metadata, records, alignment_result)

    assert actual
    assert actual.np == "NP_938033.1"
    assert actual.start == 269
    assert actual.is_full_match is True
    assert actual.nm == "NM_198291.3"
    assert actual.transcript_mode == TranscriptPriority.MANE_SELECT


@pytest.mark.asyncio
async def test_tx_scn5a(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000098-a-1"""
    urn = "urn:mavedb:00000098-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_reference(metadata, records, alignment_result)

    assert actual
    assert actual.np == "NP_000326.2"
    assert actual.start == 1619
    assert actual.is_full_match is True
    assert actual.nm == "NM_000335.5"
    assert actual.transcript_mode == TranscriptPriority.MANE_PLUS_CLINICAL


@pytest.mark.asyncio
async def test_tx_hbb(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000018-a-1"""
    urn = "urn:mavedb:00000018-a-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_reference(metadata, records, alignment_result)
    assert actual is None


@pytest.mark.asyncio
async def test_tx_app(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000113-a-2"""
    urn = "urn:mavedb:00000113-a-2"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_reference(metadata, records, alignment_result)
    assert actual is None


@pytest.mark.asyncio
async def test_tx_raf(
    scoreset_metadata_fixture: Dict[str, ScoresetMetadata],
    align_result_fixture: Dict[str, AlignmentResult],
):
    """Test transcript selection for urn:mavedb:00000061-h-1"""
    urn = "urn:mavedb:00000061-h-1"
    metadata = scoreset_metadata_fixture[urn]
    records = get_scoreset_records(urn)
    alignment_result = align_result_fixture[urn]

    actual = await select_reference(metadata, records, alignment_result)
    assert actual
    assert actual.np == "NP_002871.1"
    assert actual.start == 51
    assert actual.is_full_match is True
    assert actual.nm == "NM_002880.4"
    assert actual.transcript_mode == TranscriptPriority.MANE_SELECT
