from typing import NamedTuple, Union, Tuple, Callable, List
from enum import Enum

import numpy as np
import numpy.typing as npt


class Column(NamedTuple):
    slices: Union[slice, Tuple[slice]]
    pre_func: Callable = lambda x: x
    post_func: Callable = lambda x: x
    dtype: npt.DTypeLike = None

    @property
    def is_multi(self) -> bool:
        return isinstance(self.slices, tuple)

    @property
    def count_columns(self) -> int:
        return 1 if not self.is_multi else len(self.slices)


class Record(Enum):
    HEADER = "HEADER"
    OBSLTE = "OBSLTE"
    TITLE = "TITLE"
    SPLIT = "SPLIT"
    CAVEAT = "CAVEAT"
    COMPND = "COMPND"
    SOURCE = "SOURCE"
    KEYWDS = "KEYWDS"
    EXPDTA = "EXPDTA"
    NUMMDL = "NUMMDL"
    MDLTYP = "MDLTYP"
    AUTHOR = "AUTHOR"
    REVDAT = "REVDAT"
    SPRSDE = "SPRSDE"
    JRNL = "JRNL"
    REMARK = "REMARK"
    DBREF = {
        "columns_empty": np.concatenate(([6, 11, 13, 19, 25, 32, 41, 54, 61], np.arange(68, 80))),
        "columns": {
            "idCode": Column(slice(7, 11)),
            "chainID": Column(slice(12, 13)),
            "seqBegin": Column(slice(14, 18), dtype=int),
            "insertBegin": Column(slice(18, 19)),
            "seqEnd": Column(slice(20, 24), dtype=int),
            "insertEnd": Column(slice(24, 25)),
            "database": Column(slice(26, 32), post_func=np.char.strip),
            "dbAccession": Column(slice(33, 41), post_func=np.char.strip),
            "dbIdCode": Column(slice(42, 54), post_func=np.char.strip),
            "dbseqBegin": Column(slice(55, 60), dtype=int),
            "idbnsBeg": Column(slice(60, 61)),
            "dbseqEnd": Column(slice(62, 67), dtype=int),
            "dbinsEnd": Column(slice(67, 68)),
        },
    }
    DBREF1 = {
        "columns_empty": np.concatenate(
            ([6, 11, 13, 19, 25], np.arange(32, 47), np.arange(67, 80))
        ),
        "columns": {
            "idCode": Column(slice(7, 11)),
            "chainID": Column(slice(12, 13)),
            "seqBegin": Column(slice(14, 18), dtype=int),
            "insertBegin": Column(slice(18, 19)),
            "seqEnd": Column(slice(20, 24), dtype=int),
            "insertEnd": Column(slice(24, 25)),
            "database": Column(slice(26, 32), post_func=np.char.strip),
            "dbIdCode": Column(slice(47, 67), post_func=np.char.strip),
        },
    }

    SEQADV = "SEQADV"
    SEQRES = {
        "columns_empty": np.concatenate(
            ([6, 10, 12, 17, 18], np.arange(22, 69, 4), np.arange(70, 80))
        ),
        "columns": {
            "serNum": Column(slice(7, 10), dtype=int),
            "chainID": Column(slice(11, 12)),
            "numRes": Column(slice(13, 17), dtype=int),
            "resName": Column(
                tuple(slice(i, i+3) for i in range(19, 68, 4)),
                post_func=np.char.strip
            )
        },
    }
    MODRES = "MODRES"
    HET = "HET"
    HETNAM = "HETNAM"
    HETSYN = "HETSYN"
    FORMUL = "FORMUL"
    HELIX = "HELIX"
    SHEET = "SHEET"
    SSBOND = "SSBOND"
    LINK = "LINK"
    CISPEP = "CISPEP"
    SITE = "SITE"
    CRYST1 = "CRYST1"
    ORIGX1 = "ORIGX1"
    ORIGX2 = "ORIGX2"
    ORIGX3 = "ORIGX3"
    SCALE1 = "SCALE1"
    SCALE2 = "SCALE2"
    SCALE3 = "SCALE3"
    MTRIX1 = "MTRIX1"
    MTRIX2 = "MTRIX2"
    MTRIX3 = "MTRIX3"
    MODEL = "MODEL"
    ATOM = "ATOM"
    ANISOU = "ANISOU"
    TER = "TER"
    HETATM = "HETATM"
    ENDMDL = "ENDMDL"
    CONECT = "CONECT"
    MASTER = "MASTER"
    END = "END"

    @property
    def empty_columns(self) -> np.ndarray:
        return self.value["columns_empty"]

    @property
    def columns(self) -> dict:
        return self.value["columns"]

    @classmethod
    @property
    def names(cls) -> List:
        return cls.__dict__['_member_names_']

