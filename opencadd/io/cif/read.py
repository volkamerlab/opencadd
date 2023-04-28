"""
Read an mmCIF file from various sources.
"""

from typing import Literal, Optional, Union, Sequence
from pathlib import Path

from opencadd import _typing
from opencadd.db.pdb import file
from . import _parser as _parser


__author__ = "Armin Ariamajd"


def from_pdb_id(
    pdb_id: str,
    biological_assembly_id: Optional[int] = None,
    strictness: Literal[0, 1, 2, 3] = 0,
):
    mmcif_file = file.entry(
        pdb_id=pdb_id, file_format="cif", biological_assembly_id=biological_assembly_id
    )
    return from_file_content(content=mmcif_file, strictness=strictness)


def from_filepath(
        filepath: _typing.PathLike,
        strictness: Literal[0, 1, 2, 3] = 0,
):
    path = Path(filepath)
    with open(path, "rb") as f:
        content = f.read()
    return from_file_content(content=content, strictness=strictness)


def from_file_content(
        content: Union[str, bytes],
        strictness: Literal[0, 1, 2, 3] = 2,
) :
    if isinstance(content, bytes):
        content = content.decode()
    elif not isinstance(content, str):
        raise ValueError(
            "Parameter `content` expects either str or bytes, but the type of input argument "
            f"was: {type(content)}. Input was: {content}."
        )
    return _parser.CIFToDictParser().parse(content=content)
