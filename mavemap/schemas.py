"""Provide class definitions for commonly-used information objects."""
from collections import namedtuple
from enum import Enum


class TargetGeneCategory(str, Enum):
    """Define target gene category options. Add more definitions as needed."""

    PROTEIN_CODING = "Protein coding"


class TargetSequenceType(str, Enum):
    """Define target sequence type. Add more definitions as needed."""

    PROTEIN = "protein"
    DNA = "dna"


UniProtRef = namedtuple("UniProtRef", ["id", "offset"])
ScoresetMetadata = namedtuple(
    "ScoresetMetadata",
    [
        "urn",
        "target_gene_name",
        "target_gene_category",
        "target_sequence",
        "target_sequence_type",
        "target_reference_genome",
        "target_uniprot_ref",
    ],
)
ScoreRow = namedtuple("ScoreRow", ["hgvs_pro", "hgvs_nt", "score", "accession"])
SequenceRange = namedtuple("SequenceRange", ["start", "end"])
AlignmentResult = namedtuple(
    "AlignmentResult",
    [
        "chrom",
        "strand",
        "coverage",
        "ident_pct",
        "query_range",
        "query_subranges",
        "hit_range",
        "hit_subranges",
    ],
)
