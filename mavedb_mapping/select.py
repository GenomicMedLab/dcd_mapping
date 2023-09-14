"""Select best reference sequence."""

from cool_seq_tool import CoolSeqTool
from cool_seq_tool.data_sources.seqrepo_access import SeqRepoAccess
from gene.query import QueryHandler

from mavedb_mapping.schemas import AlignmentResult, ScoresetMetadata, TargetSequenceType


def _get_chromosome_identifier(sr: SeqRepoAccess, chromosome: str) -> str:
    """Get latest NC_ identifier given a chromosome name.

    :param sr: SeqRepo interface
    :param chromosome: prefix-free chromosome name, e.g. ``"8"``, ``"X"``
    """
    result, _ = sr.chromosome_to_acs(chromosome)
    if not result:
        raise KeyError

    sorted_results = sorted(result)
    return sorted_results[-1]


def _get_gene_symbol(q: QueryHandler, metadata: ScoresetMetadata) -> str:
    """TODO"""
    if metadata.target_uniprot_id:
        uniprot_uri = f"uniprot:{metadata.target_uniprot_id.id}"
        result = q.normalize()
    return uniprot_uri, result


def _select_protein_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> None:
    """TODO

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: TODO
    """
    cst = CoolSeqTool()
    _ = _get_chromosome_identifier(cst.seqrepo_access, align_result.chrom)
    _ = _get_gene_symbol(cst.gene_query_handler, metadata)

    """
        select *
            from uta_20210129.tx_exon_aln_v
            where hgnc = '{gsymb}'
            and {locs[i][0]} between alt_start_i and alt_end_i
            or {locs[i][1]} between alt_start_i and alt_end_i
            and alt_ac = '{chrom}'
    """


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
