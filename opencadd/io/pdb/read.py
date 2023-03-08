"""
Read a PDB file from various sources.
"""

from typing import Literal, Optional, Union, Sequence
from pathlib import Path

from opencadd import _typing
from opencadd.db.pdb import file
from . import struct, _parser, parts


__author__ = "Armin Ariamajd"


def from_pdb_id(
        pdb_id: str,
        biological_assembly_id: Optional[int] = None,
        parse_only: Optional[Sequence[Union[parts.Records, parts.Sections]]] = None,
        strictness: Literal[0, 1, 2, 3] = 0,
) -> struct.PDBStructure:
    """
    Parse a PDB file from a given PDB ID.

    This downloads the PDB file from the RCSB Protein Data Bank, and thus requires an internet connection.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    biological_assembly_id : int | str, optional, default: None
        Biological assembly ID of the entry. If not provided (i.e. when ``None``), the asymmetric unit will be
        downloaded. Notice that many records are only available in the PDB file of the asymmetric unit.
    strictness : {0, 1, 2, 3}, default: 0

    parse_only: bool, default: True
        Whether to return the fully parsed PDB file (True), or to return the parser (False) for when
        you want to only parse specific sections or records in the PDB file, or want to manually inspect
        the contents.

    Returns
    -------
    opencadd.io.pdb.struct.PDBStructure
    """
    pdb_file = file.entry(
        pdb_id=pdb_id, file_format="pdb", biological_assembly_id=biological_assembly_id
    )
    return from_file_content(content=pdb_file, strictness=strictness, parse_only=parse_only)


def from_filepath(
        filepath: _typing.PathLike,
        parse_only: Optional[Sequence[Union[parts.Records, parts.Sections]]] = None,
        strictness: Literal[0, 1, 2, 3] = 0,
):
    """
    Parse a PDB file from a local filepath.

    Parameters
    ----------
    filepath
    strictness
    parse_only

    Returns
    -------

    """
    path = Path(filepath)
    with open(path, "rb") as f:
        content = f.read()
    return from_file_content(content=content, strictness=strictness, parse_only=parse_only)


def from_url(
        url: str,
        strictness: Literal[0, 1, 2, 3] = 0,
        parse: bool = True,
) -> Union[_parser.PDBParser, struct.PDBStructure]:
    """
    Parse a PDB file from a URL.

    Parameters
    ----------
    url
    strictness
    parse

    Returns
    -------

    """
    pass


def from_file_content(
        content: Union[str, bytes],
        parse_only: Optional[Sequence[Union[parts.Records, parts.Sections]]] = None,
        strictness: Literal[0, 1, 2, 3] = 0,
) -> Union[_parser.PDBParser, struct.PDBStructure]:
    """
    Parse a PDB file from string or byte contents.

    Parameters
    ----------
    content
    strictness
    parse_only

    Returns
    -------

    """
    if isinstance(content, bytes):
        content = content.decode()
    elif not isinstance(content, str):
        raise ValueError(
            "Parameter `content` expects either str or bytes, but the type of input argument "
            f"was: {type(content)}. Input was: {content}."
        )
    if parse_only is None:
        records = (record.value for record in parts.Records)
    else:
        records = []
        for record_or_section in parse_only:
            if isinstance(record_or_section, parts.Records):
                records.append(record_or_section.value)
            elif isinstance(record_or_section, parts.Sections):
                records.extend(record_or_section.value)
            else:
                raise ValueError()
    return _parser.PDBParser(content=content, strictness=strictness).parse(records=records)
