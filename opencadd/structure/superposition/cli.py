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
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Whether to print debugging info to stdout",
    )
    p.add_argument("--no-emoji", action="store_true", default=False, help="Disable emoji logging")
    p.add_argument(
        "--method",
        default="theseus",
        help="Which alignment method to use",
        choices=[m.lower() for m in METHODS],
    )
    # p.add_argument(
    #     "--method-options",
    #     default={},
    #     type=parse_method_options,
    #     help="Options to be passed to the chosen method. Syntax is `key: value; key: value;...`",
    # )

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

    _logger.debug("Fetching reference model `%s`", reference_id)
    reference_model = Structure.from_string(reference_id)

    for i, mobile_id in enumerate(mobile_ids, 1):
        if mobile_id == reference_id and args.method == "theseus":
            _logger.error("Cannot align model `%s` against itself using theseus!", mobile_id)
            continue
        _logger.debug("Fetching mobile model #%d `%s`", i, mobile_id)
        mobile_model = Structure.from_string(mobile_id)
        _logger.debug(
            "Aligning reference `%s` and mobile `%s` with method `%s`",
            reference_id,
            mobile_id,
            args.method,
        )
        result, *_empty = align(
            [reference_model, mobile_model], method=METHODS[args.method], **args.method_options
        )

        # checks if the superposition is done, if not, there was no structural alignemnt found (for mmligner)
        if "superposed" in result:
            _logger.log(
                25,  # this the level id for results
                "RMSD for alignment #%d between `%s` and `%s` is %.1fÅ",
                i,
                reference_id,
                mobile_id,
                result["scores"]["rmsd"],
            )
            for j, structure in enumerate(result["superposed"], 1):
                structure.write(f"superposed_{args.method}_{i}_{j}.pdb")
                _logger.debug(
                    "Wrote superposed model #%d to %s", j, f"superposed_{args.method}_{i}_{j}.pdb"
                )
        else:
            _logger.log(
                25,
                "`%s` found no alignment for #%d between `%s` and `%s`.",
                args.method,
                i,
                reference_id,
                mobile_id,
            )
