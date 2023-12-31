"""Provide core MaveDB mapping methods."""
import logging
from typing import List

from mavemap.align import AlignmentError, align
from mavemap.resources import get_scoreset_metadata, get_scoreset_records
from mavemap.schemas import ScoreRow, ScoresetMetadata
from mavemap.select import TxSelectError, select_reference
from mavemap.vrs_map import VrsMapError, vrs_map

_logger = logging.getLogger(__name__)


async def map_scoreset(
    metadata: ScoresetMetadata, records: List[ScoreRow], silent: bool = True
) -> None:
    """Given information about a MAVE experiment, map to VRS.

    :param metadata: salient data gathered from scoreset on MaveDB
    :param records: experiment scoring results
    :param silent:
    :return: something (TODO)
    """
    try:
        alignment_result = align(metadata, silent)
    except AlignmentError:
        _logger.error(f"Alignment failed for scoreset {metadata.urn}")
        return None

    try:
        transcript = await select_reference(metadata, alignment_result)
    except TxSelectError:
        _logger.error(f"Transcript selection failed for scoreset {metadata.urn}")
        return None

    try:
        _ = vrs_map(metadata, transcript, records)
    except VrsMapError:
        _logger.error(f"VRS mapping failed for scoreset {metadata.urn}")


async def map_scoreset_urn(scoreset_urn: str, silent: bool = True) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param scoreset_urn: identifier for a scoreset.
    :return: something (TODO)
    """
    metadata = get_scoreset_metadata(scoreset_urn)
    records = get_scoreset_records(scoreset_urn, silent)

    mapped = await map_scoreset(metadata, records, silent)
    return mapped  # TODO
