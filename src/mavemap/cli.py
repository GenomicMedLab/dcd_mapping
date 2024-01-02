"""Provide command-line interface for accessing mapping functions."""
import asyncio
import logging

import click

from mavemap.main import map_scoreset_urn

_logger = logging.getLogger(__name__)


@click.command(no_args_is_help=True)
@click.argument("urn", nargs=1)
@click.option(
    "--debug",
    "-d",
    is_flag=True,
    show_default=True,
    default=False,
    help="Enable debug logging",
)
def cli(urn: str, debug: bool) -> None:
    """Process user commands and call core `map_scoreset()` function.

    For example:

    % python3 -m mavemap.cl 'urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    :param debug: if True, enable debug logging
    """  # noqa: D301
    logging.basicConfig(filename="mavemap.log", level=logging.INFO, force=True)
    if debug:
        _logger.setLevel(logging.DEBUG)
    else:
        _logger.setLevel(logging.INFO)
    asyncio.run(map_scoreset_urn(urn, silent=False))


if __name__ == "__main__":
    cli()
