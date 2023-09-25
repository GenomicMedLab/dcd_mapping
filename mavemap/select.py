"""Select best reference sequence."""
from typing import List

from mavemap.lookup import get_chromosome_identifier, get_gene_symbol, get_transcripts
from mavemap.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    TargetSequenceType,
)


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
    matching_transcripts = _get_matching_transcripts(metadata, align_result)


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

"""WORKING HERE


        locs = get_locs_list(item['hits'])
        chrom = get_chr(dp, item['chrom'])


        async def mapq():
            transcript_lists = []
            for i in range(len(locs)):
                testquery = (f\"""select *
                            from uta_20210129.tx_exon_aln_v
                            where hgnc = '{gsymb}'
                            and {locs[i][0]} between alt_start_i and alt_end_i
                            or {locs[i][1]} between alt_start_i and alt_end_i
                            and alt_ac = '{chrom}'\""")

                out = await utadb.execute_query(testquery)
                tl = []
                for j in range(len(out)):
                    if out[j]['tx_ac'].startswith('NR_') == False:
                        tl.append(out[j]['tx_ac'])
                if tl != []:
                    transcript_lists.append(tl)
            return(transcript_lists)

        ts = asyncio.run(mapq())
        try:
            isect = list(set.intersection(*map(set,ts)))
        except:
            try: # Look for transcripts using uniprot id
                url = 'https://www.uniprot.org/uniprot/' + str(dat.at[j, 'uniprot_id']) + '.xml'
                page = requests.get(url)
                page = BeautifulSoup(page.text)
                page = page.find_all('sequence')
                up = page[1].get_text()

                stri = str(dat.at[j,'target_sequence'])
                if up.find(stri) != -1:
                    full_match = True
                else:
                    full_match = False
                start = up.find(stri[:10])
                mappings_dict[dat.at[j,'urn']] = [dat.at[j, 'uniprot_id'], start, dat.at[j, 'urn'], full_match]
                continue
            except:
                print([dat.at[j, 'urn'], 'no transcripts found'])
                mappings_dict[dat.at[j,'urn']] = []
                continue

        mane_trans = mane.get_mane_from_transcripts(isect)
        if mane_trans != []:
            if len(mane_trans) == 1:
                np = mane_trans[0]['RefSeq_prot']
                nm = mane_trans[0]['RefSeq_nuc']
                status = 'MANE Select'
            else:
                if mane_trans[0]['MANE_status'] == 'MANE Select':
                    np = mane_trans[0]['RefSeq_prot']
                    nm = mane_trans[0]['RefSeq_nuc']
                    status = 'MANE Select'
                else:
                    np = mane_trans[1]['RefSeq_prot']
                    nm = mane_trans[1]['RefSeq_nuc']
                    status = 'MANE Plus Clinical'

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
            mappings_dict[dat.at[j,'urn']] = [np, start, dat.at[j, 'urn'], full_match, nm, status]

        else:
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
mappings_dict
"""
