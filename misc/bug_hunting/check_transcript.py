import logging
import asyncio
import pickle

from cool_seq_tool.schemas import Strand

from dcd_mapping.align import align
from dcd_mapping.resources import get_scoreset_metadata, get_scoreset_records
from dcd_mapping.transcripts import select_transcript

logging.basicConfig(
    filename="dcd-check-transcript.log",
    format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
    level=logging.WARNING,
    force=True,
)
_logger = logging.getLogger(__name__)

with open("misc/bug_hunting/human_urns.txt", "r") as f:
    urns = [line.strip() for line in f.readlines()]

with open("notebooks/analysis/results/mave_blat.pickle", "rb") as f:
    mave_blat_dict = pickle.load(f)

with open("notebooks/analysis/results/mappings.pickle", "rb") as f:
    expected_mappings = pickle.load(f)

async def check_tx_results():
    for urn in urns:
        print(f"Checking {urn}...")
        if urn != "urn:mavedb:00000083-i-1":
            continue
        try:
            # prereqs
            metadata = get_scoreset_metadata(urn)
            records = get_scoreset_records(urn)
            alignment = align(metadata, use_cached=True)
        except Exception as e:
            _logger.error("%s error before transcript selection: %s", urn, e)
            continue
        try:
            tx_result = await select_transcript(metadata, records, alignment)
        except Exception as e:
            _logger.error("%s error during transcript selection: %s", urn, e)
            continue

        breakpoint()
        if urn not in mave_blat_dict:
            continue  # not in original experiment
        if urn not in expected_mappings and tx_result is not None:
            _logger.error("%s performed transcript selection when it shouldn't have", urn)
            continue
        if urn in expected_mappings and tx_result is None:
            _logger.error("%s skipped transcript selection when it shouldn't have", urn)
            continue

        expected = expected_mappings[urn]
        try:
            for (name, actual, expected) in [
                ("NP accession", tx_result.np, expected[0]),
                ("start position", tx_result.start, expected[1]),
                ("is full match", tx_result.is_full_match, expected[3]),
                ("NM accession", tx_result.nm, expected[4]),
                ("Tx priority", tx_result.transcript_mode, expected[5])
            ]:
                if actual != expected:
                    _logger.error(
                        "%s %s mismatch: %s (actual) vs %s (expected)",
                        urn,
                        name,
                        actual,
                        expected,
                    )
        except Exception as e:
            _logger.error("%s exception: %s", urn, e)

asyncio.run(check_tx_results())
