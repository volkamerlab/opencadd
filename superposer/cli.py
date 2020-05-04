"""
Command Line Interface for superposer
"""
import os
import argparse
import logging

import atomium

from .api import align, METHODS
from ._version import get_versions as _get_versions

__version__ = _get_versions()["version"]
_logger = logging.getLogger(__name__)


def parse_cli(argv=None, greet=False):
    p = argparse.ArgumentParser(argv)
    p.add_argument(
        "structures",
        nargs="+",
        help="PDB IDs or paths to structure files to be aligned. At least two are needed. First one will be considered the target.",
    )
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Whether to print debugging info to stdout",
    )
    p.add_argument(
        "--method",
        default="theseus",
        help="Which alignment method to use",
        choices=[m.lower() for m in METHODS],
    )
    p.add_argument(
        "--method-options",
        default={},
        type=parse_method_options,
        help="Options to be passed to the chosen method. Syntax is `key: value; key: value;...`",
    )

    return p.parse_args()


def parse_method_options(string):
    options = {}
    if not string or not string.strip():
        return options
    fields = string.split(";")
    # use YAML (through ruamel_yaml) to parse each field
    for field in fields:
        # TODO: REPLACE WITH some_ruamel_yaml_function(field)  # -> {key: value}
        minidict = {"dummy": "value"}
        options.update(minidict)
    return options


def greeting():
    return (
        r"""
 ┌─┐┬ ┬┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┬─┐
 └─┐│ │├─┘├┤ ├┬┘├─┘│ │└─┐├┤ ├┬┘
 └─┘└─┘┴  └─┘┴└─┴  └─┘└─┘└─┘┴└─
 Brought to you by @volkamerlab
"""
    )[1:]


def main():
    args = parse_cli()
    debug_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=debug_level, format="%(message)s")

    _logger.log(100, "%s", greeting().format(version=__version__))

    # Delegate to the API method
    reference_id, *mobile_ids = args.structures

    opener = atomium.open if os.path.isfile(reference_id) else atomium.fetch
    reference_model = opener(reference_id).model

    for i, mobile_id in enumerate(mobile_ids, 1):
        opener = atomium.open if os.path.isfile(mobile_id) else atomium.fetch
        mobile_model = opener(mobile_id).model
        result, *_empty = align(
            [reference_model, mobile_model], method=METHODS[args.method], **args.method_options
        )
        _logger.log(
            100,
            "RMSD for alignment #%d between `%s` and `%s` is %.1fÅ",
            i,
            reference_id,
            mobile_id,
            result["scores"]["rmsd"],
        )
        for j, structure in enumerate(result["superposed"]):
            structure.write(f"superposed_{args.method}_{i}_{j}.pdb")
