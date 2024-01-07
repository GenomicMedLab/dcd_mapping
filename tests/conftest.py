"""Provide shared testing utilities."""
import json
from pathlib import Path

import pytest
from cool_seq_tool.schemas import Strand

from mavemap.schemas import AlignmentResult, ScoresetMetadata, SequenceRange


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


@pytest.fixture(scope="session")
def align_result_fixture():
    """Provide fixtures for alignment results."""
    return {
        "urn:mavedb:00000041-a-1": AlignmentResult(
            chrom="chr20",
            strand=Strand.POSITIVE,
            coverage=100.0,
            ident_pct=99.86666666666666,
            query_range=SequenceRange(start=0, end=750),
            query_subranges=[
                SequenceRange(start=0, end=52),
                SequenceRange(start=52, end=232),
                SequenceRange(start=232, end=309),
                SequenceRange(start=309, end=463),
                SequenceRange(start=463, end=595),
                SequenceRange(start=595, end=750),
            ],
            hit_range=SequenceRange(start=37397802, end=37403325),
            hit_subranges=[
                SequenceRange(start=37397802, end=37397854),
                SequenceRange(start=37400114, end=37400294),
                SequenceRange(start=37401601, end=37401678),
                SequenceRange(start=37402434, end=37402588),
                SequenceRange(start=37402748, end=37402880),
                SequenceRange(start=37403170, end=37403325),
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
        # "urn:mavedb:00000097-0-1": AlignmentResult(),
    }
