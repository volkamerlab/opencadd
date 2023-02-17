"""
Read and write Protein Data Bank (PDB) files.
"""

from typing import Literal, Optional, Union, Sequence
from pathlib import Path
from enum import Enum

from opencadd.db.pdb import file
from opencadd.io.pdb import datastruct, _parser
from opencadd import _typing


class Records(Enum):
    """
    Enumeration of records in a PDB File.
    """
    HEADER = "header"
    OBSLTE = "obslte"
    TITLE = "title"
    SPLIT = "split"
    CAVEAT = "caveat"
    COMPND = "compnd"
    SOURCE = "source"
    KEYWDS = "keywds"
    EXPDTA = "expdta"
    NUMMDL = "nummdl"
    MDLTYP = "mdltyp"
    AUTHOR = "author"
    REVDAT = "revdat"
    SPRSDE = "sprsde"
    JRNL = "jrnl"
    REMARK = "remark"
    DBREF = "dbref"
    SEQADV = "seqadv"
    SEQRES = "seqres"
    MODRES = "modres"
    HET = "het"
    HETNAM = "hetnam"
    HELIX = "helix"
    SHEET = "sheet"
    SSBOND = "ssbond"
    LINK = "link"
    CISPEP = "cispep"
    SITE = "site"
    CRYST1 = "cryst1"
    ORIGX = "origx"
    SCALE = "scale"
    MTRIX = "mtrix"
    ATOM = "atom"
    ANISOU = "anisou"
    TER = "ter"
    CONECT = "conect"


class Sections(Enum):
    """
    Enumeration of Sections in a PDB file.
    """
    Title = (
        "header",
        "obslte",
        "title",
        "split",
        "caveat",
        "compnd",
        "source",
        "keywds",
        "expdta",
        "nummdl",
        "mdltyp",
        "author",
        "revdat",
        "sprsde",
        "jrnl",
    )

def from_pdb_id(
        pdb_id: str,
        biological_assembly_id: Optional[int] = None,
        parse_only: Optional[Sequence[Union[Records, Sections]]] = None,
        strictness: Literal[0, 1, 2, 3] = 0,
) -> datastruct.PDBFile:
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
    opencadd.io.pdb.datastruct.PDBFile
    """
    pdb_file = file.entry(
        pdb_id=pdb_id, file_format="pdb", biological_assembly_id=biological_assembly_id
    )
    return from_file_content(content=pdb_file, strictness=strictness, parse_only=parse_only)


def from_filepath(
        filepath: _typing.PathLike,
        parse_only: Optional[Sequence[Union[Records, Sections]]] = None,
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
) -> Union[_parser.PDBParser, datastruct.PDBFile]:
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
        parse_only: Optional[Sequence[Union[Records, Sections]]] = None,
        strictness: Literal[0, 1, 2, 3] = 0,
) -> Union[_parser.PDBParser, datastruct.PDBFile]:
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
        records = (record.value for record in Records)
    else:
        records = []
        for record_or_section in parse_only:
            if isinstance(record_or_section, Records):
                records.append(record_or_section.value)
            elif isinstance(record_or_section, Sections):
                records.extend(record_or_section.value)
            else:
                raise ValueError()
    return _parser.PDBParser(content=content, strictness=strictness).parse(records=records)

#from_pdb_id("3w32", biological_assembly_id=1)