"""
Read an mmCIF file from various sources.
"""

from typing import Literal, Optional, Union, Sequence
from pathlib import Path

import polars as pl

from opencadd import _typing
from opencadd.db.pdb import file
from . import _parser, _filestruct


__author__ = "Armin Ariamajd"


def from_pdb_id(
    pdb_id: str,
    biological_assembly_id: Optional[int] = None,
):
    mmcif_file = file.entry(
        pdb_id=pdb_id, file_format="cif", biological_assembly_id=biological_assembly_id
    )
    return from_file_content(content=mmcif_file)


def from_filepath(
        filepath: _typing.PathLike,
):
    path = Path(filepath)
    with open(path, "rb") as f:
        content = f.read()
    return from_file_content(content=content)


def from_file_content(
        content: Union[str, bytes],
) :
    if isinstance(content, bytes):
        content = content.decode()
    elif not isinstance(content, str):
        raise ValueError(
            "Parameter `content` expects either str or bytes, but the type of input argument "
            f"was: {type(content)}. Input was: {content}."
        )
    return _parser.CIFToDictParser().parse(content=content)


def from_parquet(filepath: _typing.PathLike):
    path = Path(filepath)
    return _filestruct.DDL2CIFFile(df=pl.read_parquet(path))
