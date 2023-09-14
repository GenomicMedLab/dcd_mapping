"""Select best reference sequence."""
from typing import List

from mavemap.lookup import get_chromosome_identifier, get_gene_symbol, get_transcripts
from mavemap.schemas import AlignmentResult, ScoresetMetadata, TargetSequenceType


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


def _get_matching_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """TODO"""
    chromosome = get_chromosome_identifier(align_result.chrom)
    gene_symbol = get_gene_symbol(metadata)
    if not gene_symbol:
        raise TxSelectError
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = get_transcripts(
            gene_symbol, chromosome, hit_range.start, hit_range.end
        )
        transcript_matches.append(matches_list)
    return transcript_matches


def _select_protein_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> None:
    """TODO

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: TODO
    """
    _get_matching_transcripts(metadata, align_result)


def select_reference(metadata: ScoresetMetadata, align_result: AlignmentResult) -> None:
    """Select appropriate human reference sequence for scoreset.

    Fairly trivial for regulatory/other noncoding scoresets which report genomic
    variations.
    For protein scoresets, identify a matching RefSeq protein reference sequence.

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: TODO
    """
    if metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        _select_protein_reference(metadata, align_result)
    elif metadata.target_sequence_type == TargetSequenceType.DNA:
        pass
    else:
        raise ValueError  # TODO
