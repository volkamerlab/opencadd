"""
CLI

Command Line Interface for structuralalignment
"""
import os
import argparse
import logging

import atomium

from .api import align, METHODS


def parse_cli():
    p = argparse.ArgumentParser()
    p.add_argument(
        "structures",
        nargs="+",
        help="PDB IDs or paths to structure files to be aligned. At least two are needed. First one will be considered the target.",
    )
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
        """
   _____ _                   _                   _
  / ____| |                 | |                 | |
 | (___ | |_ _ __ _   _  ___| |_ _   _ _ __ __ _| |
  \___ \| __| '__| | | |/ __| __| | | | '__/ _` | |
  ____) | |_| |  | |_| | (__| |_| |_| | | | (_| | |
 |_____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_|
     /\   | (_)                                | |
    /  \  | |_  __ _ _ __  _ __ ___   ___ _ __ | |_
   / /\ \ | | |/ _` | '_ \| '_ ` _ \ / _ \ '_ \| __|
  / ____ \| | | (_| | | | | | | | | |  __/ | | | |_
 /_/    \_\_|_|\__, |_| |_|_| |_| |_|\___|_| |_|\__|
                __/ |
               |___/

 v{version} · Brought to you by @volkamerlab
"""
    )[1:]


def main():
    args = parse_cli()
    debug_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=debug_level)

    from ._version import get_versions

    print(greeting().format(version=get_versions()["version"]))
    # Delegate to the API method
    reference_id, *mobile_ids = args.structures

    opener = atomium.open if os.path.isfile(reference_id) else atomium.fetch
    reference_model = opener(reference_id).model

    for mobile_id in mobile_ids:
        opener = atomium.open if os.path.isfile(mobile_id) else atomium.fetch
        mobile_model = opener(mobile_id).model
        result, *_empty = align(
            [reference_model, mobile_model], method=METHODS[args.method], **args.method_options
        )
        print(
            f"RMSD for alignment between `{reference_id}` and `{mobile_id}` is {result['scores']['rmsd']}Å"
        )
