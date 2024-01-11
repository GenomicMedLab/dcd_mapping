"""Map transcripts to VRS objects."""
import logging
from typing import List, Optional, Tuple, Union

import click
from Bio.Seq import Seq
from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.core import sha512t24u
from ga4gh.vrs._internal.models import (
    Allele,
    Haplotype,
    LiteralSequenceExpression,
    SequenceLocation,
)
from ga4gh.vrs.normalize import normalize

from dcd_mapping.lookup import hgvs_to_vrs
from dcd_mapping.schemas import (
    AlignmentResult,
    ScoreRow,
    ScoresetMetadata,
    TargetType,
    TxSelectResult,
)

_logger = logging.getLogger(__name__)


class VrsMapError(Exception):
    """Raise in case of VRS mapping errors."""


def _get_sequence_id(sequence: str) -> str:
    """Make GA4GH sequence identifier

    :param sequence: raw string sequence
    :return: sequence identifier with digest (starting with ``"ga4gh:SQ."``)
    """
    sequence_id = f"ga4gh:SQ.{sha512t24u(sequence).encode('ascii')}"
    return sequence_id


def _make_dup_state() -> LiteralSequenceExpression:
    """TODO make the crazy dup sequence thing

    :return:
    """
    # 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
    raise NotImplementedError


def _update_allele_coords(allele: Allele, hits, ranges) -> Allele:  # noqa
    """TODO

    :param allele:
    :param hits:
    :return:
    """
    """
    start = allele.location.start.value
    if len(hits) == 1 and strand == 1:
        i = 0
        diff = start - hits[i][0]
        diff2 = allele.location.end.value - start
        allele.location.start.value = ranges[i][0] + diff
        allele.location.end.value = allele.location.start.value + diff2
    else:
        for i in range(len(hits)):
            if start >= hits[i][0] and start < hits[i][1]:
                break
        diff = start - hits[i][0]
        diff2 = allele.location.end.value - start
        if strand == 1: # positive orientation
            allele.location.start.value = ranges[i][0] + diff
            allele.location.end.value = allele.location.start.value + diff2
            if 'dup' in hgvs_string:
                allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
        else:
            allele.location.start.value = ranges[i][1] - diff - diff2
            allele.location.end.value = allele.location.start.value + diff2
            if 'dup' in hgvs_string:
                allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
            allele.state.sequence = str(Seq(str(allele.state.sequence)).reverse_complement())
    """
    raise NotImplementedError


def _get_haplotype_allele(
    hgvs_pro: str,
    ref_acc: str,
    offset: int,
    sequence_id: str,
    layer: AnnotationLayer,
    pre_map: bool = True,
) -> Union[Allele, Haplotype]:
    """TODO

    :param hgvs_pro: HGVS variation
    :param ref_acc: reference accession ID
    :param offset: offset from start of reference
    :param sequence_id: GA4GH sequence ID
    :param layer:
    :param pre_map: is pre-mapping (? TODO not really sure of this verbiage)
    """
    description = hgvs_pro.split(".", maxsplit=1)[-1]
    if chr(91) in description:
        description = description[1:-1]
        description_list = list(set(description.split(";")))
    else:
        description_list = [description]

    alleles = []
    for variation in description_list:
        hgvs_string = f"{ref_acc}:{layer.value}.{variation}"
        allele = hgvs_to_vrs(hgvs_string, {})  # TODO alias map

        if pre_map:
            allele.location.sequence_id = sequence_id
            if "dup" in hgvs_string:
                allele.state = _make_dup_state()
        else:
            if (
                not isinstance(allele.location, SequenceLocation)
                or not isinstance(allele.location.start, int)
                or not isinstance(allele.location.end, int)
            ):
                # shouldn't be possible anyway
                raise ValueError
            if layer == AnnotationLayer.GENOMIC:
                allele = _update_allele_coords(allele, hits, ranges)  # noqa
            else:
                allele.location.start += offset
                allele.location.end += offset
                if "dup" in hgvs_string:
                    allele.state = _make_dup_state()
        if allele.state.sequence == "N" and layer != AnnotationLayer.GENOMIC:
            # TODO
            """
            allele.state.sequence = str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])

            """
            raise NotImplementedError
        allele = normalize(allele)  # TODO inject custom data proxy?
        alleles.append(allele)

    if len(alleles) == 1:  # Not haplotype
        return alleles[0]
    else:
        return Haplotype(
            members=alleles,
            id=None,
            label=None,
            description=None,
            digest=None,
            type="Haplotype",
        )


def _map_protein_coding(
    metadata: ScoresetMetadata,
    transcript: TxSelectResult,
    records: List[ScoreRow],
) -> None:
    """TODO

    * there's a ``stri`` variable set if there are >4 unique chars in the target sequence??

    :param metadata: scoreset metadata
    :param transcript: selected transcript
    :param records: score results
    :return:
    """
    alias_map = {}
    mappings_list = []  # noqa
    var_ids_pre_map = []
    var_ids_post_map = []

    target_sequence = str(Seq(metadata.target_sequence).translate(table="1")).replace(
        "*", ""
    )
    sequence_id = _get_sequence_id(target_sequence)

    scores = [row.score for row in records]  # noqa
    accessions = [row.accession for row in records]  # noqa

    for row in records:
        if not isinstance(row.hgvs_pro, str):
            raise ValueError  # should be impossible -- TODO remove
        elif len(row.hgvs_pro) == 3 or row.hgvs_pro == "_wt" or row.hgvs_pro == "_sy":
            raise ValueError  # TODO should understand why this matters
        elif row.hgvs_pro.startswith("NP"):
            allele = hgvs_to_vrs(row.hgvs_pro, alias_map)
            var_ids_pre_map.append(allele)
            var_ids_post_map.append(allele)
        else:
            pre_allele = _get_haplotype_allele(
                row.hgvs_pro, transcript.np, 0, sequence_id, AnnotationLayer.PROTEIN
            )
            var_ids_pre_map.append(pre_allele)
            if transcript.np.startswith("N"):
                post_allele = _get_haplotype_allele(
                    row.hgvs_pro,
                    transcript.np,
                    transcript.start,
                    sequence_id,
                    AnnotationLayer.PROTEIN,
                    pre_map=False,
                )
            else:
                # TODO
                # I don't think this ever happens because it's referencing variables that
                # aren't set yet in the notebooks
                raise ValueError
            var_ids_post_map.append(post_allele)

        # TODO process NT column
    """
            tempdat = pd.DataFrame({'pre_mapping': var_ids_pre_map, 'mapped': var_ids_post_map})
            mappings_list.append(tempdat)
            scores_list.append(spro)
            accessions_list.append(accpro)

            # Process nt column if data present
            if vardat['hgvs_nt'].isnull().values.all() == False and '97' not in dat.at[i, 'urn']:
                var_ids_pre_map = []
                var_ids_post_map = []

                item = mave_blat_dict[dat.at[i, 'urn']]
                ranges = get_locs_list(item['hits'])
                hits = get_hits_list(item['hits'])
                ref = get_chr(dp, item['chrom'])
                ts = dat.at[i, 'target_sequence']
                strand = mave_blat_dict[dat.at[i, 'urn']]['strand']

                digest = 'SQ.' + sha512t24u(ts.encode('ascii'))
                alias_dict_list = [{'namespace': 'ga4gh', 'alias': digest}]
                sr.store(ts, nsaliases = alias_dict_list) # Add custom digest to SeqRepo

                ntlist = vardat['hgvs_nt']
                varm = vardat['hgvs_pro']
                sn = []
                accn = []

                for j in range(len(ntlist)):
                    if type(ntlist[j]) != str or ntlist[j] == '_wt' or ntlist[j] == '_sy':
                        continue
                    else:
                        try:
                            var_ids_pre_map.append(get_haplotype_allele(ntlist[j][2:], ref, 0, 'g', tr, dp, ts,'pre', ranges, hits, strand).as_dict())
                            var_ids_post_map.append(get_haplotype_allele(ntlist[j][2:], ref, 0, 'g', tr, dp, ts,'post', ranges, hits, strand).as_dict())
                            sn.append(scores[j])
                            accn.append(accessions[j])
                        except:
                            continue

                tempdat = pd.DataFrame({'pre_mapping': var_ids_pre_map, 'mapped': var_ids_post_map})
                mappings_list.append(tempdat)
                scores_list.append(sn)
                accessions_list.append(accn)

        vrs_mappings_dict[dat.at[i, 'urn']] = mappings_list
        scores_dict_coding[dat.at[i, 'urn']] = scores_list
        mavedb_ids_coding[dat.at[i, 'urn']] = accessions_list

    """


def _map_regulatory_noncoding(
    metadata: ScoresetMetadata,
    records: List[ScoreRow],
    align_result: AlignmentResult,
) -> List[Tuple[Allele, Allele]]:
    """Return VRS alleles representing pre- and post-mapped variation objects (?)

    :param metadata: metadata for URN
    :param records: list of MAVE experiment result rows
    :return: TODO
    """
    var_ids = []
    sequence_id = _get_sequence_id(metadata.target_sequence)

    for row in records:
        if row.hgvs_nt == "_wt" or row.hgvs_nt == "_sy":
            raise ValueError  # TODO unclear what's up here
        pre_map_allele = _get_haplotype_allele(
            row.hgvs_nt[2:], align_result.chrom, 0, sequence_id, AnnotationLayer.GENOMIC
        )  # TODO need query/hit ranges and strand for something
        if isinstance(pre_map_allele, Haplotype):
            breakpoint()  # TODO investigate
            raise NotImplementedError
        post_map_allele = _get_haplotype_allele(
            row.hgvs_nt[2:],
            align_result.chrom,
            0,
            sequence_id,
            AnnotationLayer.GENOMIC,
            False,
        )  # TODO need query/hit ranges and strand for something
        if isinstance(post_map_allele, Haplotype):
            breakpoint()  # TODO investigate
            raise NotImplementedError
        var_ids.append((pre_map_allele, post_map_allele))

    return var_ids


def vrs_map(
    metadata: ScoresetMetadata,
    align_result: AlignmentResult,
    transcript: Optional[TxSelectResult],
    records: List[ScoreRow],
    silent: bool = True,
) -> List[Tuple[Allele, Allele]]:
    """Given a description of a MAVE scoreset and an aligned transcript, generate
    the corresponding VRS objects.

    Todo:
    ----
    * Store objects in SeqRepo

    :param metadata: salient MAVE scoreset metadata
    :param transcript: output of transcript selection process
    :param records: scoreset records
    :param silent:
    :return: TODO
    """
    msg = f"Mapping {metadata.urn} to VRS..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)
    if metadata.target_gene_category == TargetType.PROTEIN_CODING and transcript:
        result = _map_protein_coding(metadata, transcript, records)
    else:
        result = _map_regulatory_noncoding(metadata, records, align_result)
    breakpoint()  # TODO tmp
    return result
