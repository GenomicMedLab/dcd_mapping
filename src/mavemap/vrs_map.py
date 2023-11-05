"""Map transcripts to VRS objects."""
from enum import Enum
from typing import Dict, List, Optional, Union
from Bio.Seq import Seq


from ga4gh.core import sha512t24u
from ga4gh.vrsatile.pydantic.vrs_models import (
    CURIE,
    Allele,
    Number,
    SequenceInterval,
    SequenceLocation,
    SimpleInterval,
    Text,
)

from mavemap.lookup import get_sequence, hgvs_to_vrs, store_sequence
from mavemap.schemas import AlignmentResult, HgvsTypePrefix, ScoreRow, ScoresetMetadata, TargetType


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _get_sequence_id(sequence: str) -> CURIE:
    """Make GA4GH sequence identifier

    :param sequence: raw string sequence
    :return: sequence digest (starting with ``"SQ."``)
    """
    sequence_id = f"ga4gh:SQ.{sha512t24u(sequence).encode('ascii')}"
    return sequence_id  # type: ignore


def _handle_hgvs_dup(sequence: str, allele: Allele) -> str:
    """Do the HGVS dup thing (TODO, not entirely clear to me what's going on)

    should reflect this line of code:
    # allele.state.sequence = 2*str(
        sr[str(allele.location.sequence_id)][
            allele.location.start.value:allele.location.end.value]
    )

    Gotta wonder if there are copy number/duplication semantics we're supposed to be
    employing here instead of this
    """
    start = allele.location.interval.start.value
    end = allele.location.interval.end.value
    return 2 * sequence[start:end]


def _get_haplotype_allele(
    hgvs_pro: str,
    np_acc: str,
    sequence_type: HgvsTypePrefix,
    sequence: str,
    offset: int = 0,
    mapped: bool = False,
    align_data: Optional[AlignmentResult] = None

) -> Union[Allele, Text]:
    """TODO

    * there's a seqrepo call in here referencing a sequence previously stored by
    its digest. I think we should just be passing this stuff directly.
    """
    base_variation = hgvs_pro.lstrip(sequence_type.value)
    if chr(91) in base_variation:  # TODO fix this, editor problem on my end
        base_variation = base_variation[1:][:-1]
        variation_list = list(set(base_variation.split(";")))
    else:
        variation_list = [base_variation]

    locations = {}
    alleles = []

    for variation in variation_list:
        hgvs_string = f"{np_acc}:{sequence_type.value}.{variation}"
        allele = hgvs_to_vrs(hgvs_string)

        if not mapped:
            sequence_id = _get_sequence_id(sequence)
            allele.location.sequence_id = sequence_id
            if "dup" in hgvs_string:
                allele.state.sequence = _handle_hgvs_dup(sequence, allele)
        else:
            if sequence_type != HgvsTypePrefix.LINEAR_GENOMIC:
                allele.location.interval.start.value += offset
                allele.location.interval.end.value += offset
                if "dup" in hgvs_string:
                    _handle_hgvs_dup(sequence, allele)
            else:
                start = allele.location.interval.start.value
                if not align_data:
                    # TODO not totally sure I'm in the right place here
                    # it seems like post-mapped protein HGVS variants might just need this?
                    raise ValueError
                else:
                    raise ValueError




def vrs_map(metadata: ScoresetMetadata, transcript: Dict, records: List[ScoreRow]):
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    Todo:
    ----
    * Store objects in SeqRepo

    :param metadata: salient MAVE scoreset metadata
    :param transcript: output of RefSeq transcript selection process
    :param records: scoreset records
    :return: TODO
    """
    if metadata.target_gene_category == TargetType.PROTEIN_CODING:
        np = transcript["np"]
        offset = transcript["start"]
        mappings_list = []
        scores_list = []
        accessions_list = []
        var_ids_pre_map = []
        var_ids_post_map = []
        spro = []  # TODO what are these
        accpro = []

        if len(set(metadata.target_sequence)) > 4:
            formatted_sequence = str(
                Seq(metadata.target_sequence).translate(table="1")
            ).replace("*", "")
        else:
            formatted_sequence = metadata.target_sequence

        target_sequence_id = (
            f"ga4gh:SQ.{sha512t24u(formatted_sequence.encode('ascii'))}"
        )

        for row in records:
            if row.hgvs_pro.startswith("NP"):
                var_ids_pre_map.append(hgvs_to_vrs(row.hgvs_pro))
                var_ids_post_map.append(hgvs_to_vrs(row.hgvs_pro))
                spro.append(row.score)
                accpro.append(row.accession)
            else:
                var_ids_pre_map.append(
                    _get_haplotype_allele(
                        row.hgvs_pro, np, HgvsTypePrefix.PROTEIN, formatted_sequence
                    )
                )
                spro.append(row.score)
                accpro.append(row.accession)
                if np.startswith("N"):
                    # TODO
                    """
                    # calling get haplotype with
                     - varm[j] = row.hgvs_pro
                     - np = np accession from transcript select
                     - 0 = ????    offset for post mapping
                     - 'p' = ????
                     - 'pre' = pre vs post mapping
                     - who knows what the empty strings are

                    var_ids_pre_map.append(get_haplotype_allele(varm[j], np, 0, 'p', tr, dp, ts, 'pre', '', '', '').as_dict())
                    var_ids_post_map.append(get_haplotype_allele(varm[j], np, offset, 'p', tr, dp, ts, 'post', '', '', '').as_dict())
                    spro.append(scores[j])
                    accpro.append(accessions[j])
                    """
                else:
                    # TODO
                    """
                    # calling get haplotype with
                    - varm[j] = row.hgvs_pro
                    var_ids_pre_map.append(get_haplotype_allele(varm[j], np, 0, 'p', tr, dp, ts, 'pre', '', '', ''))
                    var_ids_post_map.append(get_haplotype_allele(varm[j], np, offset, 'p', tr, dp, ts, 'post', ranges, hits, ''))
                    spro.append(scores[j])
                    accpro.append(accessions[j])
                    """

    else:
        raise ValueError  # TODO
