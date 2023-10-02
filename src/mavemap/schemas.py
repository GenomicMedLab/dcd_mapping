"""Provide class definitions for commonly-used information objects."""
from enum import Enum
from typing import Optional, List

from pydantic import BaseModel


class TargetGeneCategory(str, Enum):
    """Define target gene category options. Add more definitions as needed."""

    PROTEIN_CODING = "Protein coding"


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class ReferenceGenome(str, Enum):
    """Define known reference genome names."""



class UniProtRef(BaseModel):
    """Store metadata associated with MaveDB UniProt reference"""

    id: str
    offset: int



class ScoresetMetadata(BaseModel):
    """Store all relevant metadata from metadata reported for scoreset by MaveDB"""

    urn: str
    target_gene_name: Optional[str] = None
    target_gene_category: Optional[str] = None
    target_sequence: str
    target_sequence_type: TargetSequenceType
    target_reference_genome: Optional[ReferenceGenome] = None
    target_uniprot_ref: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """TODO"""

    hgvs_pro: Optional[str] = None
    hgvs_nt: Optional[str] = None
    # TODO ???
    score: str
    accession: str


class SequenceRange(BaseModel):
    """Define range over a sequence. Useful for expressing alignment query and hit results."""

    start: int
    end: int


class AlignmentResult(BaseModel):
    """Structured BLAT alignment output."""

    chrom: str
    strand: str
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
    mane_status: str
    grch38_chr: str
    chr_start: str
    chr_end: str
    chr_strand: str
