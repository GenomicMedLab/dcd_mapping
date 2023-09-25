"""Provide core MaveDB mapping methods."""
import logging
from mavemap.select import select_reference

from .align import AlignmentError, align
from .resources import get_scoreset_metadata, get_scoreset_records


_logger = logging.getLogger(__name__)


def map_scoreset(scoreset_urn: str) -> None:
    """Perform end-to-end mapping for a scoreset."""
    metadata = get_scoreset_metadata(scoreset_urn)
    _ = get_scoreset_records(scoreset_urn)

    try:
        alignment_result = align(metadata)
    except AlignmentError:
        _logger.error(f"Alignment failed for scoreset {metadata.urn}")
        return None

    select_reference(metadata, alignment_result)

    # TODO standardize to VRS
