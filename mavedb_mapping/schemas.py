"""Provide class definitions for commonly-used information objects."""
from enum import Enum
from typing import List, Optional

from pydantic import BaseModel


class TargetGeneCategory(str, Enum):
    """Define target gene category options. Add more definitions as needed."""

    PROTEIN_CODING = "Protein coding"


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


class UniProtRef(BaseModel):
    """Record UniProt annotations from scoreset metadata"""

    id: str
    offset: int


class ScoresetMetadata(BaseModel):
    """Salient metadata for a scoreset."""

    urn: str
    target_gene_name: str
    target_gene_category: TargetGeneCategory
    target_sequence: str
    target_sequence_type: TargetSequenceType
    target_reference_genome: str
    target_uniprot_id: Optional[UniProtRef] = None


class ScoreRow(BaseModel):
    """Individual result in a scoreset."""

    hgvs_pro: str
    hgvs_nt: str
    score: str
    accession: str


class QueryRange(BaseModel):
    """Query start and end"""

    start: int
    end: int


class HitRange(BaseModel):
    """Hit start and end range."""

    start: int
    end: int


class AlignmentResult(BaseModel):
    """Store relevant data from alignment process."""

    chrom: str
    strand: str
    coverage: float
    ident_pct: float
    query_range: QueryRange
    query_subranges: List[QueryRange]
    hit_range: HitRange
    hit_subranges: List[HitRange]
