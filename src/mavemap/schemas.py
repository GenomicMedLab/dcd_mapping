"""Provide class definitions for commonly-used information objects."""
from enum import StrEnum
from typing import List, Optional

from cool_seq_tool.schemas import Strand, TranscriptPriority
from pydantic import BaseModel, StrictBool, StrictInt


class TargetSequenceType(StrEnum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class ReferenceGenome(StrEnum):
    """Define known reference genome names."""

    HG38 = "hg38"
    HG19 = "hg19"
    HG16 = "hg16"


class TargetType(StrEnum):
    """Define target gene types."""

    PROTEIN_CODING = "Protein coding"
    REGULATORY = "Regulatory"
    OTHER_NC = "Other noncoding"


class UniProtRef(BaseModel):
    """Store metadata associated with MaveDB UniProt reference"""

    id: str
    offset: int


class ScoresetMetadata(BaseModel):
    """Store all relevant metadata from metadata reported for scoreset by MaveDB"""

    urn: str
    target_gene_name: str
    target_gene_category: TargetType
    target_sequence: str
    target_sequence_type: TargetSequenceType
    target_reference_genome: ReferenceGenome
    target_uniprot_ref: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """Row from a MAVE score result"""

    hgvs_pro: str
    hgvs_nt: str
    score: str
    accession: str


class SequenceRange(BaseModel):
    """Define range over a sequence. Useful for expressing alignment query and hit results."""

    start: int
    end: int


class GeneLocation(BaseModel):
    """Gene location info, gathered from normalizer result. Likely to be incomplete."""

    chromosome: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None


class AlignmentResult(BaseModel):
    """Define BLAT alignment output."""

    chrom: str
    strand: Strand
    coverage: float
    ident_pct: float
    query_range: SequenceRange
    query_subranges: List[SequenceRange]
    hit_range: SequenceRange
    hit_subranges: List[SequenceRange]


class ManeData(BaseModel):
    """Structured MANE data retrieval result."""

    ncbi_gene_id: str
    ensembl_gene_id: str
    hgnc_gene_id: str
    symbol: str
    name: str
    refseq_nuc: str
    refseq_prot: str
    ensembl_nuc: str
    ensembl_prot: str
    transcript_priority: TranscriptPriority
    grch38_chr: str
    chr_start: int
    chr_end: int
    chr_strand: str


class TxSelectResult(BaseModel):
    """Define response object from transcript selection process."""

    nm: Optional[str] = None
    np: str
    start: StrictInt
    is_full_match: StrictBool
    transcript_mode: Optional[TranscriptPriority] = None
    sequence: str
