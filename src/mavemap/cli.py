"""Provide command-line interface for accessing mapping functions."""
import logging

import click

from .main import map_scoreset_urn

_logger = logging.getLogger(__name__)


@click.command(no_args_is_help=True)
@click.argument("urn")
def cli(urn: str) -> None:
    """Process user commands and call core `map_scoreset()` function.

    For example:

    % dcd-map --urn='urn:mavedb:00000329-a-1'

    \f
    :param urn: scoreset URN
    """  # noqa: D301
    logging.basicConfig(filename="mavemap.log", level=logging.INFO, force=True)
    _logger.setLevel(logging.INFO)
    map_scoreset_urn(urn, silent=False)


if __name__ == "__main__":
    cli()
