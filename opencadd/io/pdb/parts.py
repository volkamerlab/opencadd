from enum import Enum


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
