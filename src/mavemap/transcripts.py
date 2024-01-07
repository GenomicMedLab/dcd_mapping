"""Select best reference sequence."""
import logging
import re
from typing import List, Optional

from Bio.Data.CodonTable import IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from cool_seq_tool.schemas import TranscriptPriority
from gene.database.database import click

from mavemap.lookup import (
    get_chromosome_identifier,
    get_gene_symbol,
    get_mane_transcripts,
    get_sequence,
    get_transcripts,
    get_uniprot_sequence,
)
from mavemap.schemas import (
    AlignmentResult,
    ManeData,
    ScoreRow,
    ScoresetMetadata,
    TargetSequenceType,
    TargetType,
    TxSelectResult,
)

_logger = logging.getLogger(__name__)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


async def _get_compatible_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """Acquire matching transcripts"""
    if align_result.chrom.startswith("chr"):
        aligned_chrom = align_result.chrom[3:]
    else:
        aligned_chrom = align_result.chrom
    chromosome = get_chromosome_identifier(aligned_chrom)
    gene_symbol = get_gene_symbol(metadata)
    if not gene_symbol:
        raise TxSelectError
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = await get_transcripts(
            gene_symbol, chromosome, hit_range.start, hit_range.end
        )
        if matches_list:
            transcript_matches.append(matches_list)
    return transcript_matches


def _reduce_compatible_transcripts(matching_transcripts: List[List[str]]) -> List[str]:
    """Reduce list of list of transcripts to a list containing only entries present
    in each sublist

    :param matching_transcripts: list of list of transcript accession IDs
    :return: list of transcripts shared by all sublists
    """
    common_transcripts_set = set(matching_transcripts[0])
    for sublist in matching_transcripts[1:]:
        common_transcripts_set.intersection_update(sublist)
    common_transcripts = list(common_transcripts_set)
    return common_transcripts


def _choose_best_transcript(mane_transcripts: List[ManeData], urn: str) -> ManeData:
    """Choose best transcript (Select > Plus Clinical) given MANE status

    Todo:
    ----
     * Handle case where it's empty (ie implement longest compatible logic)

    :param mane_transcripts: list of MANE transcript descriptions
    :param urn: scoreset ID for error message
    :return: best transcript
    """
    if len(mane_transcripts) == 2:
        if mane_transcripts[0].transcript_priority == TranscriptPriority.MANE_SELECT:
            return mane_transcripts[0]
        else:
            _logger.warning(
                "Unexpected transcript priority return for %s: %s",
                urn,
                mane_transcripts,
            )
            breakpoint()  # TODO want to take a look at this when it comes up
            return mane_transcripts[1]
    elif len(mane_transcripts) == 1:
        return mane_transcripts[0]
    else:
        # TODO pretty sure here's where we need to handle the stuff at the end of this
        # code block fetching protein accessions
        """
        trans_lens = []
        for i in range(len(isect)):
            trans_lens.append(len(str(sr[isect[i]])))
        loc = trans_lens.index(max(trans_lens))
        nm = isect[loc]

        testquery = f"SELECT pro_ac FROM uta_20210129.associated_accessions WHERE tx_ac = '{nm}'"
        async def np():
            out = await utadb.execute_query(testquery)
            try:
                return out[0]['pro_ac']
            except:
                return out
        np = asyncio.run(np())

        if np != []:
            oseq = dat.at[j, 'target_sequence']

            if len(set(str(oseq))) > 4:
                stri = str(oseq)
            else:
                oseq = Seq(oseq)
                stri = str(oseq.translate(table=1)).replace('*', '')

            if str(sr[np]).find(stri) != -1:
                full_match = True
            else:
                full_match = False
            start = str(sr[np]).find(stri[:10])
            mappings_dict[dat.at[j,'urn']] = [np, start, dat.at[j, 'urn'], full_match, nm, 'Longest Compatible']
        """
        _logger.error(
            f"Unexpected number of MANE transcripts: {len(mane_transcripts)}, urn: {urn}"
        )
        raise NotImplementedError


def _get_protein_sequence(target_sequence: str) -> str:
    """Get protein sequence if necessary.

    It'd be nice if there was a more elegant way to check if the sequence was already a
    protein sequence.

    :param target_sequence: sequence set as baseline in MAVE experiment (might already
        be set to protein)
    :return: resulting protein sequence
    """
    if len(set(target_sequence)) > 4:
        protein_sequence = target_sequence
    else:
        protein_sequence = str(Seq(target_sequence).translate(table="1")).replace(
            "*", ""
        )
    return protein_sequence


async def _select_protein_reference(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> TxSelectResult:
    """Select preferred transcript for protein reference sequence

    :param metadata: Scoreset metadata from MaveDB
    :param align_result: alignment results
    :return: Best transcript and associated metadata
    :raise TxSelectError: if no matching MANE transcripts and unable to get UniProt ID/
    reference sequence
    """
    matching_transcripts = await _get_compatible_transcripts(metadata, align_result)
    if not matching_transcripts:
        common_transcripts = None
    else:
        common_transcripts = _reduce_compatible_transcripts(matching_transcripts)
    if not common_transcripts:
        if not metadata.target_uniprot_ref:
            raise TxSelectError(
                f"Unable to find matching transcripts for {metadata.urn}"
            )
        protein_sequence = get_uniprot_sequence(metadata.target_uniprot_ref.id)
        np_accession = metadata.target_uniprot_ref.id
        ref_sequence = get_uniprot_sequence(metadata.target_uniprot_ref.id)
        if not ref_sequence:
            raise ValueError(
                f"Unable to grab reference sequence from uniprot.org for {metadata.urn}"
            )
        nm_accession = None
        tx_mode = None
    else:
        mane_transcripts = get_mane_transcripts(common_transcripts)
        best_tx = _choose_best_transcript(mane_transcripts, metadata.urn)
        ref_sequence = get_sequence(best_tx.refseq_prot)
        nm_accession = best_tx.refseq_nuc
        np_accession = best_tx.refseq_prot
        tx_mode = best_tx.transcript_priority

    protein_sequence = _get_protein_sequence(metadata.target_sequence)
    # TODO -- look at these two lines
    is_full_match = ref_sequence.find(protein_sequence) != -1
    start = ref_sequence.find(protein_sequence[:10])

    protein_mapping_info = TxSelectResult(
        nm=nm_accession,
        np=np_accession,
        start=start,
        is_full_match=is_full_match,
        sequence=protein_sequence,
        transcript_mode=tx_mode,
    )
    return protein_mapping_info


def _offset_target_sequence(metadata: ScoresetMetadata, records: List[ScoreRow]) -> int:
    """Find start location in target sequence (it's not always 0)

    :param metadata:
    :param records:
    :return: starting index position (may be 0)
    """
    if not isinstance(records[0].hgvs_pro, str) or records[0].hgvs_pro.startswith("NP"):
        return 0  # TODO explain?
    protein_change_list = [rec.hgvs_pro.lstrip("p.") for rec in records]

    aa_dict = {}  # TODO name?
    for protein_change in protein_change_list:
        if protein_change == "_sy" or protein_change == "_wt":
            raise ValueError
        if ";" in protein_change:
            protein_changes = protein_change[1:-1].split(";")
        else:
            protein_changes = [protein_change]
        for change in protein_changes:
            aa = change[:3]
            if aa == "=" or change[-3:] not in IUPACData.protein_letters_3to1:
                continue
            if "=" in change:
                loc = change[3:-1]
            else:
                loc = change[3:-3]
            if loc not in aa_dict:
                loc = re.sub("[^0-9]", "", loc)
                if loc:
                    aa_dict[loc] = seq1(aa)

    err_locs = []
    protein_sequence = Seq(metadata.target_sequence).translate(table="1")
    for i in range(len(protein_sequence)):
        if str(i) in aa_dict and aa_dict[str(i)] != protein_sequence[i - 1]:
            err_locs.append(i)
    if len(err_locs) == 0:
        return 0
    aa_dict = {int(k): v for k, v in aa_dict.items()}
    aa_dict = sorted(aa_dict.items())
    aa_dict = dict(aa_dict)  # TODO what are we doing here

    """
    if len(err_locs) > 1:
        aa_dict = {int(k):v for k,v in aa_dict.items()}
        aa_dict = sorted(aa_dict.items())
        aa_dict = dict(aa_dict)
        locs = list(aa_dict.keys())[0:5]
        p0, p1, p2, p3, p4 = locs[0], locs[1], locs[2], locs[3], locs[4]
        offset = locs[0]

        seq = ''
        for key in aa_dict:
            seq = seq + aa_dict[key]

        for i in range(len(ts)):
            if ts[i] == aa_dict[p0] and ts[i + p1 - p0] == aa_dict[p1] and ts[i + p2 - p0] == aa_dict[p2] and ts[i + p3 - p0] == aa_dict[p3] and ts[i + p4 - p0] == aa_dict[p4]:
                if i + 1 == min(aa_dict.keys()) or i + 2 == min(aa_dict.keys()):
                    offset_within_ts[urn] = 0
                else:
                    offset_within_ts[urn] = i
                break
    """
    # this part probably goes after the result in ``select_reference()``
    """
    for key in offset_within_ts:
        mappings_dict[key][1] = offset_within_ts[key]
    """
    raise NotImplementedError


async def select_transcript(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
    silent: bool = False,
) -> Optional[TxSelectResult]:
    """Select appropriate human reference sequence for scoreset.

    * Unnecessary for regulatory/other noncoding scoresets which report genomic
    variations.
    * For protein scoresets, identify a matching RefSeq protein reference sequence.

    :param metadata: Scoreset metadata from MaveDB
    :param records:
    :param align_result: alignment results
    :return: Transcript description (accession ID, offset, selected sequence, etc)
    """
    msg = f"Selecting reference sequence for {metadata.urn}..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)

    if metadata.target_gene_category == TargetType.PROTEIN_CODING:
        output = await _select_protein_reference(metadata, align_result)
        if output and metadata.target_sequence_type == TargetSequenceType.DNA:
            if (
                metadata.urn == "urn:mavedb:00000053-a-1"
                or metadata.urn == "urn:mavedb:00000053-a-2"
            ):
                # target sequence is missing codon
                return None
            output.start = _offset_target_sequence(metadata, records)
    else:
        # can't provide transcripts for regulatory/noncoding scoresets
        return None

    msg = "Reference selection complete."
    if not silent:
        click.echo(msg)
    return output
