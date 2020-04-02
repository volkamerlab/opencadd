"""
CLI

Command Line Interface for structuralalignment
"""
import os
import argparse

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
        minidict = {
            "dummy": "value"
        }  # REPLACE WITH some_ruamel_yaml_function(field)  # -> {key: value}
        options.update(minidict)
    return options


def main():
    args = parse_cli()

    # Build atomium models (IO)
    models = []
    for structure in args.structures:
        opener = atomium.open if os.path.isfile(structure) else atomium.fetch
        models.append(opener(structure).model)

    # Delegate to the API method
    return align(*models, method=METHODS[args.method], **args.method_options)
