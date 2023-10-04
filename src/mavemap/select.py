"""Select best reference sequence."""
import logging
from typing import List

from Bio.Seq import Seq

from mavemap.lookup import (
    get_chromosome_identifier,
    get_gene_symbol,
    get_mane_transcript,
    get_reference_sequence,
    get_transcripts,
)
from mavemap.schemas import (
    AlignmentResult,
    ManeStatus,
    ScoresetMetadata,
    TargetSequenceType,
)

_logger = logging.getLogger(__name__)


class TxSelectError(Exception):
    """Raise for transcript selection failure."""


def _get_compatible_transcripts(
    metadata: ScoresetMetadata, align_result: AlignmentResult
) -> List[List[str]]:
    """Acquire matching transcripts"""
    chromosome = get_chromosome_identifier(align_result.chrom)
    gene_symbol = get_gene_symbol(metadata)
    if not gene_symbol:
        raise TxSelectError
    transcript_matches = []
    for hit_range in align_result.hit_subranges:
        matches_list = get_transcripts(
            gene_symbol, chromosome, hit_range.start, hit_range.end
        )
        if matches_list:
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
    matching_transcripts = _get_compatible_transcripts(metadata, align_result)

    # drop transcripts not associated with all hit subranges
    common_transcripts_set = set(matching_transcripts[0])
    for sublist in matching_transcripts[1:]:
        common_transcripts_set.intersection_update(sublist)
    common_transcripts = list(common_transcripts_set)

    # TODO if this fails, look up transcripts on uniprot?

    mane_transcripts = get_mane_transcript(common_transcripts)

    # TODO handle case where this is empty

    if len(mane_transcripts) == 2:
        if mane_transcripts[0].mane_status == ManeStatus.SELECT:
            mane_tx = mane_transcripts[0]
        else:
            mane_tx = mane_transcripts[1]
    elif len(mane_transcripts) == 1:
        mane_tx = mane_transcripts[0]
    else:
        _logger.error(
            f"Unexpected number of MANE transcripts: {len(mane_transcripts)}, urn: {metadata.urn}"
        )
        raise TxSelectError
    np, nm, status = mane_tx.refseq_prot, mane_tx.refseq_nuc, mane_tx.mane_status

    # TODO there's a thing here about taking the sequence as-is if it contains
    # more than four unique chars, that seems off
    # Check specific chars used instead?
    if len(set(metadata.target_sequence)) > 4:
        protein_sequence = metadata.target_sequence
    else:
        protein_sequence = str(
            Seq(metadata.target_sequence).translate(table="1")
        ).replace("*", "")

    ref_sequence = get_reference_sequence(np)
    is_full_match = ref_sequence.find(protein_sequence)
    start = ref_sequence.find(protein_sequence[:10])  # TODO seems potentially sus?


def select_reference(metadata: ScoresetMetadata, align_result: AlignmentResult) -> None:
    """Select appropriate human reference sequence for scoreset.

    Fairly trivial for regulatory/other noncoding scoresets which report genomic
    variations.
    For protein scoresets, identify a matching RefSeq protein reference sequence.
    More description here TODO.

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
## Helper functions

def get_start(string):
    return int(string.split(':')[0].strip('['))

def get_end(string):
    return int(string.split(':')[1].strip(']'))

def get_locs_list(hitsdat):
    locs_list = []
    for i in range(len(hitsdat.index)):
        start = get_start(hitsdat.at[i, 'hit_ranges'])
        end = get_end(hitsdat.at[i, 'hit_ranges'])
        locs_list.append([start,end])
    return locs_list

def get_hits_list(hitsdat):
    hits_list = []
    for i in range(len(hitsdat.index)):
        start = get_start(hitsdat.at[i, 'query_ranges'])
        end = get_end(hitsdat.at[i, 'query_ranges'])
        hits_list.append([start,end])
    return hits_list

def get_query_hits(dat):
    query_list = []
    hits_list = []
    for i in range(len(dat.index)):
        query_start = get_start(dat.at[i, 'query_ranges'])
        query_end = get_end(dat.at[i, 'query_ranges'])
        query_list.append([query_start, query_end])
        hit_start = get_start(dat.at[i, 'hit_ranges'])
        hit_end = get_end(dat.at[i, 'hit_ranges'])
        hits_list.append([hit_start, hit_end])
        return query_list, hits_list

def get_ga4gh(dp, ref):
    aliases = dp.get_metadata(ref)['aliases']
    f = filter(lambda x: 'ga4gh' in x, aliases)
    return 'ga4gh:' + list(f)[0].split(':')[1]

def get_chr(dp, chrom):
    aliases = dp.get_metadata('GRCh38:' + chrom)['aliases']
    f = filter(lambda x: 'refseq' in x, aliases)
    return list(f)[0].split(':')[1]

def modify_hgvs(var, ref, off, hp):
    if len(var) == 3 or var == '_wt' or var == '_sy' or '[' in var:
        return var
    var = ref + ':' + var
    var = hp.parse_hgvs_variant(var)
    var.posedit.pos.start.base = var.posedit.pos.start.base + off
    var.posedit.pos.end.base = var.posedit.pos.end.base + off
    return(str(var))

def blat_check(i):
    item = mave_blat_dict[dat.at[i, 'urn']]
    if item['uniprot'] == None:
        test = dat.at[i, 'target'].split(' ')
        for j in range(len(test)):
            try:
                out = qh.normalize(test[j]).gene_descriptor
                gene_dat = [out.label, out.extensions[2].value['chr']]
                if item['chrom'] != gene_dat[1]:
                    return False
                else:
                    return True
            except:
                continue

def get_haplotype_allele(var, ref, offset, l, tr, dp, ts, mapped, ranges, hits, strand):
    var = var.lstrip(f'{l}.')

    if '[' in var:
        var = var[1:][:-1]
        varlist = var.split(';')
        varlist = list(set(varlist))
    else:
        varlist = list()
        varlist.append(var)

    locs = {}
    alleles = []

    for i in range(len(varlist)):
        try:
            hgvs_string = ref + ':'+ l +'.' + varlist[i]
            allele = tr.translate_from(hgvs_string, 'hgvs')

            if mapped == 'pre':
                allele.location.sequence_id = 'ga4gh:SQ.' + sha512t24u(ts.encode('ascii'))
                if 'dup' in hgvs_string:
                    allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])

            else:
                if l != 'g':
                    allele.location.start.value = allele.location.start.value + offset
                    allele.location.end.value = allele.location.end.value + offset
                    if 'dup' in hgvs_string:
                        allele.state.sequence = 2*str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])

                else:
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

            if allele.state.sequence == 'N' and l != 'p':
                allele.state.sequence = str(sr[str(allele.location.sequence_id)][allele.location.start.value:allele.location.end.value])
            allele = normalize(allele, data_proxy = dp)
            allele.id = ga4gh_identify(allele)
            alleles.append(allele)
        except:
            vrstext = {'definition':ref + ':'+ l +'.' + varlist[i], 'type': 'Text'}
            return vrstext

    if len(alleles) == 1: # Not haplotype
        return alleles[0]
    else:
        return models.Haplotype(members = alleles)

def get_clingen_id(hgvs):
    url = 'https://reg.genome.network/allele?hgvs=' + hgvs
    page = requests.get(url).json()
    page = page['@id']
    try:
        return page.split('/')[4]
    except:
        return 'NA'

## UTA Transcript Selection
nest_asyncio.apply()
mane = MANETranscriptMappings()
utadb = UTADatabase(db_pwd = 'uta')
qh = QueryHandler()
dp = SeqRepoDataProxy(sr = sr)

mappings_dict = {}
mave_dat = pd.read_csv('results/mave_dat.csv')
dat = mave_dat
with open('results/mave_blat.pickle', 'rb') as fn:
    mave_blat_dict = pickle.load(fn)

for j in range(len(dat.index)):
    if dat.at[j, 'target_type'] == 'Protein coding' or dat.at[j, 'target_type'] == 'protein_coding':
        item = mave_blat_dict[dat.at[j,'urn']]
        #if blat_check(j) == False:
         #   mappings_dict[dat.at[j, 'urn']] = 'BLAT hit not found on correct chromosome'
          #  continue
        if item['chrom'] == 'NA':
            continue

        # hit subranges
        locs = get_locs_list(item['hits'])
        chrom = get_chr(dp, item['chrom'])  # probably an accession

        # get gene symbol from uniprot/target name data
        try:
            uniprot = dat.at[j, 'uniprot_id']
            gsymb = qh.normalize(str(f'uniprot:{uniprot}')).gene_descriptor.label
        except:
            temp = dat.at[j, 'target'].split(' ')
            if temp[0] == 'JAK':
                temp[0] = 'JAK1'
            gsymb = qh.normalize(temp[0]).gene_descriptor.label


        # for each hit subrange start/end pair, get compatible transcript accessions
        # drop any with NR_ prefixes
        # if non empty, append to overall list
        # return overall list of lists
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
