"""Align MaveDB target sequences to a human reference genome."""
import logging
import subprocess
import uuid
from pathlib import Path
from typing import Any, Dict, Generator, Optional

from Bio.SearchIO import HSP
from Bio.SearchIO import read as read_blat
from Bio.SearchIO._model import Hit, QueryResult
from cool_seq_tool.schemas import Strand
from gene.database.database import click

from mavemap.lookup import get_chromosome_identifier, get_gene_location
from mavemap.resources import get_mapping_tmp_dir, get_ref_genome_file
from mavemap.schemas import (
    AlignmentResult,
    GeneLocation,
    ScoresetMetadata,
    SequenceRange,
    TargetSequenceType,
)

_logger = logging.getLogger(__name__)


class AlignmentError(Exception):
    """Raise when errors encountered during alignment."""


def _build_query_file(scoreset_metadata: ScoresetMetadata) -> Generator[Path, Any, Any]:
    """Construct BLAT query file.

    TODO double-check that yield behaves the way I think it does

    This function is broken out to enable mocking while testing.

    :param scoreset_metadata: MaveDB scoreset metadata object
    :return: Yielded Path to constructed file. Deletes file once complete.
    """
    query_file = (
        get_mapping_tmp_dir() / f"blat_query_{scoreset_metadata.urn}_{uuid.uuid1()}.fa"
    )
    _logger.debug("Writing BLAT query to %", query_file)
    with open(query_file, "w") as f:
        f.write(">" + "query" + "\n")
        f.write(scoreset_metadata.target_sequence + "\n")
        f.close()
    yield query_file
    query_file.unlink()


def _run_blat_command(command: str, args: Dict) -> subprocess.CompletedProcess:
    """Execute BLAT binary with relevant params.

    This function is broken out to enable mocking while testing.

    Currently, we rely on a system-installed BLAT binary accessible in the containing
    environment's PATH. This is sort of awkward and it'd be nice to make use of some
    direct bindings or better packaging if that's possible.

    * Perhaps `gget`? https://pachterlab.github.io/gget/en/blat.html
    * ``PxBlat``? https://github.com/ylab-hi/pxblat

    :param command: shell command to execute
    :param args: ``subprocess.run`` extra args (eg redirecting output for silent mode)
    :return: process result
    """
    _logger.debug("Running BLAT command: %", command)
    return subprocess.run(command, shell=True, **args)


def _get_blat_output(
    scoreset_metadata: ScoresetMetadata, query_file: Path, silent: bool
) -> QueryResult:
    """Run a BLAT query and returns a path to the output object.

    We create query and output files in the application's "temporary" folder, which
    should be deleted by the process once complete. This happens manually, but we could
    probably add a decorator or a context manager for a bit more elegance.

    Ideally, we should see if we could pipe query output to STDOUT and then grab/process
    it that way instead of using a temporary intermediary file.

    :param scoreset_metadata: object containing scoreset attributes
    :param query_file: Path to BLAT query file
    :param silent: suppress BLAT command output
    :return: BLAT query result
    :raise AlignmentError: if BLAT subprocess returns error code
    """
    reference_genome_file = get_ref_genome_file(
        silent=silent
    )  # TODO hg38 by default--what about earlier builds?
    out_file = (
        get_mapping_tmp_dir() / f"blat_out_{scoreset_metadata.urn}_{uuid.uuid1()}.psl"
    )

    if scoreset_metadata.target_sequence_type == TargetSequenceType.PROTEIN:
        target_commands = "-q=prot -t=dnax"
    elif scoreset_metadata.target_sequence_type == TargetSequenceType.DNA:
        target_commands = "-q=dnax -t=dnax"
    else:
        query_file.unlink()
        out_file.unlink()
        raise AlignmentError(
            f"Unknown target sequence type: {scoreset_metadata.target_sequence_type} for scoreset {scoreset_metadata.urn}"
        )
    command = f"blat {reference_genome_file} {target_commands} -minScore=20 {query_file} {out_file}"
    if silent:
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
    # TODO
    # the notebooks handle errors here by trying different BLAT arg configurations --
    # investigate, refer to older code if it comes up
    output = read_blat(out_file.absolute(), "blat-psl")

    # clean up
    query_file.unlink()
    out_file.unlink()

    return output


def _get_best_hit(output: QueryResult, urn: str, chromosome: Optional[str]) -> Hit:
    """Get best hit from BLAT output.

    First, try to return hit corresponding to expected chromosome taken from scoreset
    metadata. If chromosome doesn't match any of the outputs or is unavailable, take
    the hit with the single highest-scoring HSP.

    :param output: BLAT output
    :param urn: scoreset URN to use in error messages
    :param chromosome: refseq chromosome ID, e.g. ``"NC_000001.11"``
    :return: best Hit
    :raise AlignmentError: if unable to get hits from output
    """
    if chromosome:
        if chromosome.startswith("refseq"):
            chromosome = chromosome[7:]

        for hit in output:
            hit_chr = hit.id
            if hit_chr.startswith("chr"):
                hit_chr = hit_chr[3:]
            hit_chr_ac = get_chromosome_identifier(hit_chr)
            if hit_chr_ac == chromosome:
                return hit
        else:
            if list(output):
                hit_chrs = [h.id for h in output]
                _logger.warning(
                    f"Failed to match hit chromosomes during alignment. URN: "
                    f"{urn}, expected chromosome: {chromosome}, hit chromosomes: {hit_chrs}"
                )

    best_score = 0
    best_score_hit = None
    for hit in output:
        best_local_score = max(hit, key=lambda i: i.score).score
        if best_local_score > best_score:
            best_score = best_local_score
            best_score_hit = hit

    if best_score_hit is None:
        _logger.error(f"Couldn't get hits from {urn} -- check BLAT output.")
        raise AlignmentError

    return best_score_hit


def _get_best_hsp(hit: Hit, urn: str, gene_location: Optional[GeneLocation]) -> HSP:
    """Retrieve preferred HSP from BLAT Hit object.

    If gene location data is available, prefer the HSP with the least distance
    between the start of the hit and the start coordinate of the gene. Otherwise,
    take the HSP with the highest score value.

    :param hit: hit object from BLAT result
    :param urn: scoreset identifier for use in error messages
    :param gene_location: location data acquired by normalizing scoreset metadata
    :return: Preferred HSP object
    :raise AlignmentError: if hit object appears to be empty (should be impossible)
    """
    best_hsp = None
    if gene_location and gene_location.start is not None:
        best_hsp = min(hit, key=lambda hsp: abs(hsp.hit_start - gene_location.start))
    else:
        best_hsp = max(hit, key=lambda hsp: hsp.score)
    if best_hsp is None:
        _logger.error(
            f"Unable to get best HSP from hit -- this should be impossible? urn: {urn}, hit: {str(hit)}"
        )
        raise AlignmentError
    return best_hsp


def _get_best_match(output: QueryResult, metadata: ScoresetMetadata) -> AlignmentResult:
    """Obtain best high-scoring pairs (HSP) object for query sequence.

    :param metadata: scoreset metadata
    :param output: BLAT result object
    :return: alignment result ??
    """
    location = get_gene_location(metadata)
    if location:
        chromosome = location.chromosome
    else:
        chromosome = None
    best_hit = _get_best_hit(output, metadata.urn, chromosome)
    best_hsp = _get_best_hsp(best_hit, metadata.urn, location)

    strand = Strand.POSITIVE if best_hsp[0].query_strand == 1 else Strand.NEGATIVE
    coverage = 100 * (best_hsp.query_end - best_hsp.query_start) / output.seq_len  # type: ignore
    identity = best_hsp.ident_pct  # type: ignore
    chrom = best_hsp.hit_id

    query_subranges = []
    hit_subranges = []
    for fragment in best_hsp:
        query_subranges.append(
            SequenceRange(start=fragment.query_start, end=fragment.query_end)
        )
        hit_subranges.append(
            SequenceRange(start=fragment.hit_start, end=fragment.hit_end)
        )

    result = AlignmentResult(
        chrom=chrom,
        strand=strand,
        ident_pct=identity,
        coverage=coverage,
        query_range=SequenceRange(start=best_hsp.query_start, end=best_hsp.query_end),
        query_subranges=query_subranges,
        hit_range=SequenceRange(start=best_hsp.hit_start, end=best_hsp.hit_end),
        hit_subranges=hit_subranges,
    )
    return result


def align(scoreset_metadata: ScoresetMetadata, silent: bool = True) -> AlignmentResult:
    """Align target sequence to a reference genome.

    :param scoreset_metadata: object containing scoreset metadata
    :param quiet: suppress BLAT process output if true
    :return: data wrapper containing alignment results
    """
    msg = f"Performing alignment for {scoreset_metadata.urn}..."
    if not silent:
        click.echo(msg)
    _logger.info(msg)
    query_file = next(_build_query_file(scoreset_metadata))
    blat_output = _get_blat_output(scoreset_metadata, query_file, silent)

    match = _get_best_match(blat_output, scoreset_metadata)
    msg = "Alignment complete."
    if not silent:
        click.echo(msg)
    _logger.info(msg)
    return match
