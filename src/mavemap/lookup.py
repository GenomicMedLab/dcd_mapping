"""Handle API lookups to external (non-MaveDB) services.

This module should contain methods that we don't want to think about caching.
"""
import logging
from typing import Dict, List, Optional

import requests
from cool_seq_tool.app import CoolSeqTool
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import TranscriptPriority
from ga4gh.core._internal.models import Extension, Gene
from ga4gh.vrs._internal.models import Allele, SequenceLocation
from ga4gh.vrs.extras.translator import Translator
from gene.database import create_db
from gene.query import QueryHandler
from gene.schemas import SourceName

from mavemap.schemas import GeneLocation, ManeData, ScoresetMetadata

_logger = logging.getLogger(__name__)

# ---------------------------------- Global ---------------------------------- #


class CoolSeqToolBuilder:
    """Singleton constructor for ``cool-seq-tool`` instance."""

    def __new__(cls) -> CoolSeqTool:
        """Provide ``CoolSeqTool`` instance. Construct it if unavailable.

        :return: singleton instance of CoolSeqTool
        """
        if not hasattr(cls, "instance"):
            cls.instance = CoolSeqTool()
        return cls.instance


def get_seqrepo() -> SeqRepoAccess:
    """Retrieve SeqRepo access instance."""
    cst = CoolSeqToolBuilder()
    return cst.seqrepo_access


class GeneNormalizerBuilder:
    """Singleton constructor for Gene Normalizer instance."""

    def __new__(cls) -> QueryHandler:
        """Provide Gene Normalizer instance. Construct it if unavailable.

        :return: singleton instance of ``QueryHandler`` for Gene Normalizer
        """
        if not hasattr(cls, "instance"):
            db = create_db()
            q = QueryHandler(db)
            cls.instance = q
        return cls.instance


class VrsTranslatorBuilder:
    """Singleton constructor for VRS-Python translator instance."""

    # TODO this looks.... very wrong?
    def __new__(cls) -> Translator:
        """Provide VRS-Python translator. Construct if unavailable.

        :return: singleton instances of Translator
        """
        if not hasattr(cls, "instance"):
            cst = CoolSeqToolBuilder()
            cls.instance = Translator(cst.seqrepo_access, normalize=False)
        return cls.instance


# ----------------------------------- UTA ----------------------------------- #


async def get_protein_accession(transcript: str) -> Optional[str]:
    """Retrieve protein accession for a transcript.

    :param transcript: transcript accession, e.g. ``"NM_002529.3"``
    :return: protein accession if successful
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT pro_ac FROM {uta.schema}.associated_accessions
    WHERE tx_ac = '{transcript}'
    """
    result = await uta.execute_query(query)
    if result:
        return result[0]["pro_ac"]


async def get_transcripts(
    gene_symbol: str, chromosome_ac: str, start: int, end: int
) -> List[str]:
    """Get transcript accessions matching given parameters (excluding non-coding RNA).

    TODO: may be able to successfully query with only one of gene symbol/chromosome ac.

    :param gene_symbol: HGNC-given gene symbol (usually, but not always, equivalent to
        symbols available in other nomenclatures.)
    :param chromosome: chromosome accession (e.g. ``"NC_000007.13"``)
    :param start: starting position
    :param end: ending position
    :return: candidate transcript accessions
    """
    uta = CoolSeqToolBuilder().uta_db
    query = f"""
    SELECT tx_ac
    FROM {uta.schema}.tx_exon_aln_v
    WHERE hgnc = '{gene_symbol}'
      AND ({start} BETWEEN alt_start_i AND alt_end_i OR {end} BETWEEN alt_start_i AND alt_end_i)
      AND alt_ac = '{chromosome_ac}'
      AND tx_ac NOT LIKE 'NR_%';
    """
    result = await uta.execute_query(query)
    return [row["tx_ac"] for row in result]


def get_mane_transcripts(transcripts: List[str]) -> List[ManeData]:
    """Get corresponding MANE data for transcripts.

    :param transcripts: candidate transcripts list
    :return: complete MANE descriptions
    """
    mane = CoolSeqToolBuilder().mane_transcript_mappings
    mane_transcripts = mane.get_mane_from_transcripts(transcripts)
    mane_data = []
    for result in mane_transcripts:
        mane_data.append(
            ManeData(
                ncbi_gene_id=result["#NCBI_GeneID"],
                ensembl_gene_id=result["Ensembl_Gene"],
                hgnc_gene_id=result["HGNC_ID"],
                symbol=result["symbol"],
                name=result["name"],
                refseq_nuc=result["RefSeq_nuc"],
                refseq_prot=result["RefSeq_prot"],
                ensembl_nuc=result["Ensembl_nuc"],
                ensembl_prot=result["Ensembl_prot"],
                transcript_priority=TranscriptPriority(
                    "_".join(result["MANE_status"].lower().split())
                ),
                grch38_chr=result["GRCh38_chr"],
                chr_start=result["chr_start"],
                chr_end=result["chr_end"],
                chr_strand=result["chr_strand"],
            )
        )
    return mane_data


# ------------------------------ Gene Normalizer ------------------------------ #


def _get_hgnc_symbol(term: str) -> Optional[str]:
    """Fetch HGNC symbol from gene term.

    :param term: gene referent
    :return: gene symbol if available
    """
    q = GeneNormalizerBuilder()
    result = q.normalize_unmerged(term)
    hgnc = result.source_matches.get(SourceName.HGNC)
    if hgnc and len(hgnc.records) > 0:
        # probably fine to just use first match
        return hgnc.records[0].symbol
    return None


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


def _normalize_gene(term: str) -> Optional[Gene]:
    """Fetch normalizer response for gene term.

    :param term: gene name or referent to normalize
    :return: normalized Gene if successful
    """
    q = GeneNormalizerBuilder()
    response = q.normalize(term)
    if response.match_type > 0:
        return response.gene
    else:
        return None


def _get_normalized_gene_response(
    metadata: ScoresetMetadata,
) -> Optional[Gene]:
    """Fetch best normalized concept given available scoreset metadata.

    :param metadata: salient scoreset metadata items
    :return: Normalized gene if available
    """
    if metadata.target_uniprot_ref:
        gene_descriptor = _normalize_gene(metadata.target_uniprot_ref.id)
        if gene_descriptor:
            return gene_descriptor

    # try taking the first word in the target name
    if metadata.target_gene_name:
        parsed_name = metadata.target_gene_name.split(" ")[0]
        gene_descriptor = _normalize_gene(parsed_name)
        if gene_descriptor:
            return gene_descriptor

    return None


def _get_genomic_interval(
    extensions: List[Extension], src_name: str
) -> Optional[GeneLocation]:
    """Extract start/end coords from extension list. Extensions in normalized genes
    can be of several different types, but we only want SequenceLocation data.

    :param extensions: extensions given in a descriptor
    :return: genomic interval if available
    """
    locations = [ext for ext in extensions if f"{src_name}_locations" in ext.name]
    if locations and len(locations[0].value) > 0:  # type: ignore
        location_values = [
            v
            for v in locations[0].value
            if v["type"] == "SequenceLocation"  # type: ignore
        ]
        if location_values:
            return GeneLocation(
                start=location_values[0]["start"],
                end=location_values[0]["end"],
                chromosome=get_chromosome_identifier_from_vrs_id(
                    f"ga4gh:{location_values[0]['sequenceReference']['refgetAccession']}"
                ),
            )
    return None


def get_gene_location(metadata: ScoresetMetadata) -> Optional[GeneLocation]:
    """Acquire gene location data from gene normalizer using metadata provided by
    scoreset.

    As with ``get_gene_symbol()``, we try to normalize from the following:
    1. UniProt ID, if available
    2. Target name: specifically, we try the first word in the name (this could
    cause some problems and we should double-check it)

    :param metadata: data given by MaveDB API
    :return: gene location data if available
    """
    gene_descriptor = _get_normalized_gene_response(metadata)
    if not gene_descriptor or not gene_descriptor.extensions:
        return None

    hgnc_locations: List[Extension] = [
        loc for loc in gene_descriptor.extensions if loc.name == "hgnc_locations"
    ]
    if hgnc_locations and len(hgnc_locations[0].value) > 0:  # type: ignore
        return GeneLocation(chromosome=hgnc_locations[0].value[0].chr)  # type: ignore

    for src_name in ("ensembl", "ncbi"):
        loc = _get_genomic_interval(gene_descriptor.extensions, src_name)
        if loc:
            return loc

    return None


# --------------------------------- SeqRepo --------------------------------- #
# TODO
# * some of these could be refactored into a single method
# * not clear if all of them are necessary
# * either way, they should all be renamed once we have a final idea of what's needed


def get_chromosome_identifier(chromosome: str) -> str:
    """Get latest NC_ identifier given a chromosome name.

    :param chromosome: prefix-free chromosome name, e.g. ``"8"``, ``"X"``
    :return: latest ID if available
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.chromosome_to_acs(chromosome)
    if not result:
        raise KeyError

    # make sure e.g. version .10 > version .9
    sorted_results = sorted(result, key=lambda i: int(i.split(".")[-1]))
    return sorted_results[-1]


def get_ucsc_chromosome_name(chromosome: str) -> str:
    """Get UCSC/GENCODE-style chromosome name, eg ``"chr1"`` instead of ``"1"`` or
    ``"NC_000001.11"``.

    :param chromosome: chromosome name/identifier
    :return: UCSC/GENCODE-style chromosome name
    :raise KeyError: if unable to find matching name
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(chromosome, "GRCh38")
    if not result:
        raise KeyError

    sorted_results = sorted([r for r in result if "chr" in r])
    try:
        return sorted_results[-1].split(":")[1]
    except IndexError:
        raise KeyError


def get_chromosome_identifier_from_vrs_id(sequence_id: str) -> Optional[str]:
    """Get NC_ identifier given a VRS sequence ID.

    :param sequence_id: identifier a la ``ga4gh:SQ.XXXXXX``
    :return: NC_ chromosome ID
    :raise KeyError: if unable to retrieve identifier
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    result, _ = sr.translate_identifier(sequence_id, "refseq")
    if not result:
        raise KeyError

    sorted_results = sorted(result)
    return sorted_results[-1]


def get_sequence(
    sequence_id: str, start: Optional[int] = None, end: Optional[int] = None
) -> str:
    """Get reference sequence given a sequence identifier.

    :param sequence_id: sequence identifier, e.g. ``"NP_938033.1"``
    :return: sequence
    :raise KeyError: if lookup fails
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    try:
        sequence = sr.get_sequence(sequence_id, start, end)
    except (KeyError, ValueError):
        _logger.error(f"Unable to acquire sequence for ID: {sequence_id}")
        raise KeyError
    if sequence is None:
        _logger.error(f"Unable to acquire sequence for ID: {sequence_id}")
        raise KeyError
    return sequence


def store_sequence(sequence: str, names: List[Dict]) -> None:
    """Store sequnce in SeqRepo.

    I'm a little queasy about this part -- it seems potentially dangerous to be
    modifying state outside of the mapper library itself, particularly if there
    are any needs for those changes to endure (and if there aren't, why are we
    modifying outside state in the first place?).

    Currently unused unless we really really need this functionality.

    :param sequence: raw sequence
    :param names: list of namespace/alias pairs,
        e.g. ``{"namespace": "GA4GH", "alias": "SQ.XXXXXX"}
    """
    sr = CoolSeqToolBuilder().seqrepo_access
    sr.sr.store(sequence, nsaliases=names)


# -------------------------------- VRS-Python -------------------------------- #


def hgvs_to_vrs(hgvs: str, alias_map: Dict) -> Allele:
    """Convert HGVS variation description to VRS object.

    # TODO incorporate alias map

    :param hgvs: MAVE-HGVS variation string
    :param alias_map: lookup for custom sequence IDs
    :return: Corresponding VRS allele as a Pydantic class
    """
    tr = VrsTranslatorBuilder()
    vrs_allele = tr.translate_from(hgvs, "hgvs")
    allele = Allele(**vrs_allele)

    if (
        not isinstance(allele.location, SequenceLocation)
        or not isinstance(allele.location.start, int)
        or not isinstance(allele.location.end, int)
    ):
        raise ValueError

    return allele


# ---------------------------------- Misc. ---------------------------------- #


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
