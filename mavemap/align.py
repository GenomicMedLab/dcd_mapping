"""Align MaveDB target sequences to a human reference genome."""
import logging
import subprocess
from pathlib import Path
from typing import Any, Dict, Generator, Optional

from Bio.SearchIO import read as read_blat
from Bio.SearchIO._model import QueryResult

from mavemap.resources import get_mapping_tmp_dir, get_ref_genome_file
from mavemap.schemas import (
    AlignmentResult,
    ScoresetMetadata,
    SequenceRange,
    TargetSequenceType,
)

class AlignmentError(Exception):
    """Raise when errors encountered during alignment."""


def _build_query_file(scoreset_metadata: ScoresetMetadata) -> Generator[Path, Any, Any]:
    """Construct BLAT query file.

    TODO double-check that yield behaves the way I think it does

    This function is broken out to enable mocking while testing.

    :param scoreset_metadata: MaveDB scoreset metadata object
    :return: Yielded Path to constructed file. Deletes file once complete.
    """
    query_file = get_mapping_tmp_dir() / "blat_query.fa"
    with open(query_file, "w") as f:
        f.write(">" + "query" + "\n")
        f.write(scoreset_metadata.target_sequence + "\n")
        f.close()
    yield query_file
    query_file.unlink()


def _run_blat_command(command: str, args: Dict) -> subprocess.CompletedProcess:
    """Execute BLAT binary with relevant params.

    Currently, we rely on a system-installed BLAT binary accessible in the containing
    environment's PATH. This is sort of awkward and it'd be nice to make use of some
    direct bindings or better packaging if that's possible.

    This function is broken out to enable mocking while testing.

    :param command: shell command to execute
    :param args: ``subprocess.run`` extra args (eg redirecting output for silent mode)
    :return: process result
    """
    return subprocess.run(command, shell=True, **args)


# TODO make output object an arg
def _get_blat_output(scoreset_metadata: ScoresetMetadata, query_file: Path, quiet: bool) -> QueryResult:
    """Run a BLAT query and returns a path to the output object.

    We create query and output files in the application's "temporary" folder, which
    should be deleted by the process once complete. This happens manually, but we could
    probably add a decorator or a context manager for a bit more elegance.

    :param scoreset_metadata: object containing scoreset attributes
    :param query_file: Path to BLAT query file
    :param quiet: suppress BLAT command output
    :return: BLAT query result
    :raise AlignmentError: if BLAT subprocess returns error code
    """
    reference_genome_file = get_ref_genome_file()
    # TODO is this min score value correct?
    # min_score = len(scoreset_metadata.target_sequence) // 2  # minimum match 50%
    min_score = 20
    out_file = get_mapping_tmp_dir() / "blat_out.psl"

    if scoreset_metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        command = f"blat {reference_genome_file} -q=prot -t=dnax -minScore={min_score} {query_file} {out_file}"
    elif scoreset_metadata.target_sequence_type == TargetSequenceType.DNA:
        command = f"blat {reference_genome_file} -q=dnax -t=dnax -minScore={min_score} {query_file} {out_file}"
    else:
        query_file.unlink()
        out_file.unlink()
        raise AlignmentError(
            f"Unknown target sequence type: {scoreset_metadata.target_sequence_type} for scoreset {scoreset_metadata.urn}"
        )
    if quiet:
        kwargs = {"stdout": subprocess.DEVNULL, "stderr": subprocess.STDOUT}
    else:
        kwargs = {}
    process = _run_blat_command(command, kwargs)
    if process.returncode != 0:
        query_file.unlink()
        out_file.unlink()
        raise AlignmentError(
            f"BLAT process returned error code {process.returncode}: {command}"
        )

    # the notebooks handle errors here by trying different BLAT arg configurations --
    # investigate, refer to older code if it comes up
    output = read_blat(out_file.absolute(), "blat-psl")

    # clean up
    query_file.unlink()
    out_file.unlink()

    return output


def _get_best_match(output: QueryResult) -> AlignmentResult:
    """Obtain best high-scoring pairs (HSP) object for query sequence.

    Initially, we do this naively (take first instance of best base score). In the
    future, we should try to match against items pulled from scoreset metadata, like
    UniProt accession info or locations of a named gene.

    :param output: BLAT result object
    """
    best_score = 0
    best_hsp = None
    for hit in output:
        for hsp in hit:
            if hsp.score > best_score:
                best_score = hsp.score
                best_hsp = hsp

    if best_hsp is None:
        raise AlignmentError("No parseable BLAT query results found")

    hsp = best_hsp
    chrom = hsp.hit_id.strip("chr")
    strand = hsp[0].query_strand
    coverage = 100 * (hsp.query_end - hsp.query_start) / output.seq_len  # type: ignore
    ident_pct = hsp.ident_pct

    query_subranges = []
    hit_subranges = []

    for hsp_fragment in hsp:
        query_subranges.append(
            SequenceRange(start=hsp_fragment.query_start, end=hsp_fragment.query_end)
        )
        hit_subranges.append(
            SequenceRange(start=hsp_fragment.hit_start, end=hsp_fragment.hit_end)
        )

    result = AlignmentResult(
        chrom=chrom,
        strand=strand,
        ident_pct=ident_pct,
        coverage=coverage,
        query_range=SequenceRange(start=hsp.query_start, end=hsp.query_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=hsp.hit_start, end=hsp.hit_end),
        hit_subranges=hit_subranges,
    )
    return result


def align(
    scoreset_metadata: ScoresetMetadata, quiet: bool = True
) -> AlignmentResult:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :param quiet: suppress BLAT process output if true
    :return: data wrapper containing alignment results
    """
    query_file = next(_build_query_file(scoreset_metadata))
    blat_output = _get_blat_output(scoreset_metadata, query_file, quiet)

    match = _get_best_match(blat_output)
    return match
