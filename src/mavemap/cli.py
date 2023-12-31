"""Provide command-line interface for accessing mapping functions."""
import asyncio
import logging

import click

from mavemap.main import map_scoreset_urn

_logger = logging.getLogger(__name__)


@click.command(no_args_is_help=True)
@click.argument("urn", nargs=1)
async def cli(urn: str) -> None:
    """Process user commands and call core `map_scoreset()` function.

    For example:

    % python3 -m mavemap.cl 'urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    """  # noqa: D301
    logging.basicConfig(filename="mavemap.log", level=logging.INFO, force=True)
    _logger.setLevel(logging.INFO)
    await map_scoreset_urn(urn, silent=False)


if __name__ == "__main__":
    asyncio.run(cli())
