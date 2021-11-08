"""
Command Line Interface for opencadd
"""
import argparse
import logging

from .api import align, METHODS
from ..core import Structure
from ...utils import PerLevelFormatter, EmojiPerLevelFormatter
from ..._version import get_versions as _get_versions

__version__ = _get_versions()["version"]
_logger = logging.getLogger(__name__)


def parse_cli(argv=None, greet=False):
    p = argparse.ArgumentParser(argv)
    p.add_argument(
        "structures",
        nargs="+",
        help="PDB IDs or paths to structure files to be aligned. At least two are needed. First one will be considered the target.",
    )
    p.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__)
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Whether to print debugging info to stdout",
    )
    p.add_argument(
        "--no-emoji", 
        action="store_true", 
        default=False,
        help="Disable emoji logging",
    )
    p.add_argument(
        "--method",
        default="theseus",
        help="Which alignment method to use",
        choices=[m.lower() for m in METHODS],
    )
    p.add_argument(
        "--select",
        default={},
        type=parse_selection,
        help="MDAnalysis like selection for the structures. Syntax is `Selection; Seletion`, Example:`backbone; backbone`",
    )
    p.add_argument(
        "--method-options",
        default={},
        type=parse_method_options,
        help="Options to be passed to the chosen method. Syntax is `key: value; key: value;...`",
    )

    return p.parse_args()


def parse_method_options(string):
    if not string or not string.strip():
        return {}
    else:
        options = dict(option.split(": ") for option in string.split("; "))
        return options


def parse_selection(string):
    if not string:
        return []
    else:
        selection = string.split("; ")
        return selection


def greeting():
    return (
        r"""
 ┌─┐┬ ┬┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┬─┐
 └─┐│ │├─┘├┤ ├┬┘├─┘│ │└─┐├┤ ├┬┘
 └─┘└─┘┴  └─┘┴└─┴  └─┘└─┘└─┘┴└─
 Brought to you by @volkamerlab
"""
    )[1:]


def configure_logger(level=logging.INFO, formatter=EmojiPerLevelFormatter):
    logger = logging.getLogger("opencadd")
    logger.setLevel(level)
    handler = logging.StreamHandler()
    formatter_instance = formatter()
    handler.setFormatter(formatter_instance)
    logger.addHandler(handler)


def main():
    args = parse_cli()
    # Configure logging
    formatter = PerLevelFormatter if args.no_emoji else EmojiPerLevelFormatter
    level = logging.DEBUG if args.verbose else logging.INFO
    configure_logger(level, formatter)
    _logger.log(101, greeting())

    # Delegate to the API method
    reference_id, *mobile_ids = args.structures

    _logger.debug(f"Fetching reference model `{reference_id}`")
    reference_model = Structure.from_string(reference_id)

    for i, mobile_id in enumerate(mobile_ids, 1):
        if mobile_id == reference_id and args.method == "theseus":
            _logger.error(f"Cannot align model `{mobile_id}` against itself using theseus!")
            continue
        _logger.debug(f"Fetching mobile model #{i} `{mobile_id}`")
        mobile_model = Structure.from_string(mobile_id)
        _logger.debug(
            f"Aligning reference `{reference_id}` and mobile `{mobile_id}` with method `{args.method}`",

        )
        result, *_empty = align(
            [reference_model, mobile_model],
            method=METHODS[args.method],
            user_select=args.select,
            **args.method_options,
        )

        # checks if the superposition is done, if not, there was no structural alignemnt found (for mmligner)
        if "superposed" in result:
            _logger.log(
                25,  # this the level id for results
                f"results for alignment #{i} between `{reference_id}` and `{mobile_id}`:\n"
                f"RMSD: {round(result['scores']['rmsd'], 3)} Å \n"
                f"coverage: {result['scores']['coverage']} \n"
                # size of selections; if no selection provided, size of structures
                f"lenght of selections: {reference_id} has {result['metadata']['reference_size']} residues; "
                f"{mobile_id} has {result['metadata']['mobile_size']} residues \n"
            )
            for j, structure in enumerate(result["superposed"], 1):
                structure.write(f"superposed_{args.method}_{i}_{j}.pdb")
                _logger.debug(
                    "Wrote superposed model #%d to %s", j, f"superposed_{args.method}_{i}_{j}.pdb"
                )
        else:
            _logger.log(
                25,
                f"`{args.method}` found no alignment for #{i} between `{reference_id}` and `{mobile_id}`."
            )
