"""Provide core MaveDB mapping methods."""
import logging
from typing import List
from mavemap.schemas import ScoreRow, ScoresetMetadata
from mavemap.select import select_reference

from .align import AlignmentError, align
from .resources import get_scoreset_metadata, get_scoreset_records


_logger = logging.getLogger(__name__)


def map_scoreset(metadata: ScoresetMetadata, records: List[ScoreRow]) -> None:
    """Given information about a MAVE experiment, map to VRS.

    TODO actually do this.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :return: something (TODO)
    """
    try:
        alignment_result = align(metadata)
    except AlignmentError:
        _logger.error(f"Alignment failed for scoreset {metadata.urn}")
        return None

    select_reference(metadata, alignment_result)

    # TODO standardize to VRS



def map_scoreset_urn(scoreset_urn: str) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param scoreset_urn: identifier for a scoreset.
    :return: something (TODO)
    """
    metadata = get_scoreset_metadata(scoreset_urn)
    records = get_scoreset_records(scoreset_urn)

    return (map_scoreset(metadata, records))
