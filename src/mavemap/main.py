"""Provide core MaveDB mapping methods."""
import logging
from typing import List

import click

from mavemap.align import AlignmentError, align
from mavemap.resources import (
    ResourceAcquisitionError,
    get_scoreset_metadata,
    get_scoreset_records,
)
from mavemap.schemas import ScoreRow, ScoresetMetadata
from mavemap.transcripts import TxSelectError, select_transcript
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

    print("Alignment result:")  # TODO remove these print calls
    print(alignment_result)

    try:
        transcript = await select_transcript(
            metadata, records, alignment_result, silent
        )
    except TxSelectError:
        _logger.error(f"Transcript selection failed for scoreset {metadata.urn}")
        return None

    print("Transcript:")
    print(transcript)

    try:
        _ = vrs_map(metadata, transcript, records)
    except VrsMapError:
        _logger.error(f"VRS mapping failed for scoreset {metadata.urn}")


async def map_scoreset_urn(urn: str, silent: bool = True) -> None:
    """Perform end-to-end mapping for a scoreset.

    :param urn: identifier for a scoreset.
    :return: something (TODO)
    """
    try:
        metadata = get_scoreset_metadata(urn)
        records = get_scoreset_records(urn, silent)
    except ResourceAcquisitionError as e:
        msg = f"Unable to acquire resource from MaveDB: {e}"
        _logger.critical(msg)
        click.echo(f"Error: {msg}")
        return None
    mapped = await map_scoreset(metadata, records, silent)
    return mapped  # TODO not sure what this will be
