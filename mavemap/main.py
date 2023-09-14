"""Provide core MaveDB mapping methods."""
from mavemap.select import select_reference

from .align import align
from .resources import get_scoreset_metadata, get_scoreset_records


def map_scoreset(scoreset_urn: str) -> None:
    """Perform end-to-end mapping for a scoreset."""
    metadata = get_scoreset_metadata(scoreset_urn)
    _ = get_scoreset_records(scoreset_urn)

    alignment_result = align(metadata)
    if not alignment_result:
        return  # TODO handle fail

    select_reference(metadata, alignment_result)

    # TODO standardize to VRS
