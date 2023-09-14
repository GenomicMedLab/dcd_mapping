"""Handle API lookups to external (non-MaveDB) services.

This module should contain methods that we don't want to think about caching.
"""
import asyncio
from typing import List, Optional

import requests
from cool_seq_tool import CoolSeqTool

from mavemap.schemas import ManeData, ScoresetMetadata


def get_clingen_id(hgvs: str) -> Optional[str]:
    """Fetch ClinGen ID. TODO finish this.

    :param hgvs: HGVS ID todo ??
    :return: ClinGen ID if available
    :raise HTTPError: if request encounters an error
    """
    url = f"https://reg.genome.network/allele?hgvs={hgvs}"
    response = requests.get(url)
    response.raise_for_status()
    page = response.json()
    page = page["@id"]
    return page.split("/")[4]


def get_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """Get sequence directly from UniProt.

    :param uniprot_id: ID provided with target info
    :return: transcript accession if successful
    :raise HTTPError: if response comes with an HTTP error code
    """
    url = f"https://www.ebi.ac.uk/proteins/api/proteins?accession={uniprot_id.split(':')[1]}&format=json"
    response = requests.get(url)
    response.raise_for_status()
    json = response.json()
    return json[0]["sequence"]["sequence"]


class CoolSeqToolBuilder:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Provide ``CoolSeqTool`` instance. Construct it if unavailable."""
        if not hasattr(cls, "instance"):
            cls.instance = CoolSeqTool()
        return cls.instance


def get_protein_accession(transcript: str) -> Optional[str]:
    """Retrieve protein accession for a transcript.

    :param transcript: transcript accession, e.g. ``"NM_002529.3"``
    :return: protein accession if successful
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT pro_ac FROM {uta.schema}.associated_accessions
    WHERE tx_ac = '{transcript}'
    """
    result = asyncio.run(uta.execute_query(query))
    if result:
        return result[0]["pro_ac"]


def get_transcripts(
    gene_symbol: str, chromosome: str, start: int, end: int
) -> List[str]:
    """Get transcript accessions matching given parameters.

    TODO: may be able to successfully query with only one of gene symbol/chromosome ac.

    :param gene_symbol: HGNC-given gene symbol (usually, but not always, equivalent to
        symbols available in other nomenclatures.)
    :param chromosome: chromosome accession (e.g. ``"NC_000007.13"``)
    :param start: starting position
    :param end: ending position
    :return: TODO idk
    """
    uta = CoolSeqToolBuilder().uta_db
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
    q = CoolSeqToolBuilder().gene_query_handler
    result = q.normalize_unmerged(term)
    if "HGNC" in result.source_matches:
        hgnc = result.source_matches["HGNC"]  # type: ignore
        if len(hgnc.records) > 0:
            # probably fine to just use first match
            return hgnc.records[0].symbol


def get_gene_symbol(metadata: ScoresetMetadata) -> Optional[str]:
    """Acquire HGNC gene symbol given provided metadata from scoreset.

    Right now, we use two sources for normalizing:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

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


def get_mane_transcripts(transcripts: List[str]) -> List[ManeData]:
    """Get corresponding MANE data for transcripts.


    :param transcripts: candidate transcripts
    :return: complete MANE descriptions
    """
    mane = CoolSeqToolBuilder().mane_transcript_mappings
    mane_transcripts = mane.get_mane_from_transcripts(transcripts)
    return [ManeData(*list(r.values())) for r in mane_transcripts]


def get_chromosome_identifier(chromosome: str) -> str:
    """Get latest NC_ identifier given a chromosome name.

    :param chromosome: prefix-free chromosome name, e.g. ``"8"``, ``"X"``
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.chromosome_to_acs(chromosome)
    if not result:
        raise KeyError

    sorted_results = sorted(result)
    return sorted_results[-1]
