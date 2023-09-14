"""Handle API lookups to external (non-MaveDB) services.

This module should contain methods that we don't think are safely cache-able.
"""
import asyncio
from typing import List, Optional

import requests
from cool_seq_tool import CoolSeqTool

from mavemap.schemas import AlignmentResult, ScoresetMetadata


def get_clingen_id(hgvs: str) -> Optional[str]:
    """Fetch ClinGen ID.

    :param hgvs: HGVS ID todo ??
    :return: ClinGen ID if available
    :raise HTTPError: if request encounters an error
    """
    response = requests.get(f"https://reg.genome.network/allele?hgvs={hgvs}")
    try:
        response.raise_for_status()
    except requests.HTTPError:
        return None
    page = response.json()
    page = page["@id"]
    return page.split("/")[4]


class CoolSeqToolFactory:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Define construction of new instance."""
        if not hasattr(cls, "instance"):
            cls.instance = CoolSeqTool()
        return cls.instance


def _get_transcripts(
    gene_symbol: str, chromosome: str, start: int, end: int
) -> List[str]:
    """Get transcript accessions matching given parameters.

    TODO: need to set UTA schema as env var.
    TODO: may be able to successfully query with only one of gene symbol/chromosome ac.

    :param gene_symbol: HGNC-given gene symbol (usually, but not always, equivalent to
        symbols available in other nomenclatures.)
    :param chromosome: chromosome accession (e.g. ``"NC_000007.13"``)
    :param start: starting position
    :param end: ending position
    :return: TODO idk
    """
    uta = CoolSeqToolFactory().uta_db
    query = f"""
    SELECT tx_ac
    FROM {uta.schema}.tx_exon_aln_v
    WHERE hgnc = '{gene_symbol}'
    AND {start} BETWEEN alt_start_i AND alt_end_i
    OR {end} BETWEEN alt_start_i AND alt_end_i
    AND alt_ac = '{chromosome}';
    """
    result = asyncio.run(uta.execute_query(query))
    return [row["tx_ac"] for row in result.items()]


def _get_hgnc_symbol(term: str) -> Optional[str]:
    """Fetch HGNC symbol from gene term.

    :param term: gene referent
    :return: gene symbol if available
    """
    q = CoolSeqToolFactory().gene_query_handler
    result = q.normalize_unmerged(term)
    if "HGNC" in result.source_matches:
        hgnc = result.source_matches["HGNC"]  # type: ignore
        if len(hgnc.records) > 0:
            # probably fine to just use first match
            return hgnc.records[0].symbol


def get_gene_symbol(metadata: ScoresetMetadata) -> Optional[str]:
    """Acquire HGNC gene symbol given provided metadata from scoreset.

    :param ScoresetMetadata: data given by MaveDB API
    :return: gene symbol if available
    """
    if metadata.target_uniprot_ref:
        result = _get_hgnc_symbol(metadata.target_uniprot_ref.id)
        if result:
            return result

    # try taking the first word in the target name
    if metadata.target_gene_name:
        parsed_name = metadata.target_gene_name.split(" ")[0]
        return _get_hgnc_symbol(parsed_name)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


def get_matching_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """TODO"""
    gene_symbol = get_gene_symbol(metadata)
    if not gene_symbol:
        raise TxSelectError
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = _get_transcripts(
            gene_symbol, align_result.chrom, hit_range.start, hit_range.end
        )
        transcript_matches.append(matches_list)
    return transcript_matches
