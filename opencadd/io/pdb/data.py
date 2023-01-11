from typing import List
from enum import Enum
import numpy as np
import pandas as pd

from opencadd.io.parsing import Column, MultiColumn
from opencadd.io.pdb import fields


class Record(Enum):

    HEADER = (
        True,
        True,
        True,
        {
            "classification": Column((10, 50), fields.header_classification),
            "deposition_date": Column((50, 59), fields.Date.from_pdb),
            "pdb_id": Column((62, 66)),
        },
        np.concatenate((np.arange(6, 10), np.arange(59, 62), np.arange(66, 80))),
    )

    OBSLTE = (
        False,
        True,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "replacement_date": Column((11, 20), fields.Date.from_pdb),
            "pdb_id": Column((21, 25)),
            "replaced_pdb_ids": MultiColumn(
                ((31, 35), (36, 40), (41, 45), (46, 50), (51, 55), (56, 60), (61, 65), (66, 70),
                 (71, 75))
            )
        },
        np.concatenate(
            ([6, 7, 10, 20], np.arange(25, 31), np.arange(35, 71, 5), np.arange(74, 80))
        )
    )

    TITLE = (
        True,
        True,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "title": Column((10, 80), fields.String.from_pdb)
        },
        np.array([6, 7])
    )

    SPLIT = (
        False,
        True,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "pdb_ids": MultiColumn(
                ((11, 15), (16, 20), (21, 25), (26, 30), (31, 35), (36, 40), (41, 45), (46, 50),
                 (51, 55), (56, 60), (61, 65), (66, 70), (71, 75), (76, 80))
            )
        }
    )

    CAVEAT = (
        False,
        True,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "pdb_id": Column((11, 15)),
            "comment": Column((19, 79), fields.String.from_pdb),
        }
    )

    COMPND = (
        True,
        True,
        False,
        {
            "continuation": Column((7, 10)),
            "compound": Column(
                (10, 80),
                lambda arr: fields.compound_compound(fields.SpecificationList.from_pdb(arr))
            )
        }
    )

    SOURCE = (
        True,
        True,
        False,
        {
            "continuation": Column((7, 10)),
            "source_name": Column(
                (10, 79),
                lambda arr: fields.source_source_name(fields.SpecificationList.from_pdb(arr)))
        }
    )

    KEYWDS = (
        True,
        True,
        False,
        {
            "continuation": Column((8, 10)),
            "keywords": Column((10, 79), fields.List.from_pdb)
        }
    )

    EXPDTA = (
        True,
        True,
        False,
    )

    NUMMDL = (
        False,
        True,
        True,

    )

    MDLTYP = (
        False,
        True,
        False,
    )

    AUTHOR = (
        True,
        True,
        False,

    )

    REVDAT = (
        True,
        False,
        True,
    )

    SPRSDE = (
        False,
        True,
        False,
    )

    JRNL = (
        False,
    )

    REMARK = "REMARK"

    DBREF = (
        False,
        False,
        True,
        {
            # "pdb_id": Column((7, 11)),
            "chain_id": Column((12, 13)),
            "residue_num_begin": Column((14, 18), fields.Integer.from_pdb),
            "residue_icode_begin": Column((18, 19)),
            "residue_num_end": Column((20, 24), fields.Integer.from_pdb),
            "residue_icode_end": Column((24, 25)),
            "db": Column((26, 32), np.char.strip),
            "db_chain_accession": Column((33, 41), np.char.strip),
            "db_chain_id": Column((42, 54), np.char.strip),
            "db_residue_num_begin": Column((55, 60), fields.Integer.from_pdb),
            "db_residue_icode_begin": Column((60, 61)),
            "db_residue_num_end": Column((62, 67), fields.Integer.from_pdb),
            "db_residue_icode_end": Column((67, 68)),
        },
        np.concatenate(([6, 11, 13, 19, 25, 32, 41, 54, 61], np.arange(68, 80))),
    )

    DBREF1 = (
        False,
        False,
        True,
        {
            # "pdb_id": Column((7, 11)),
            "chain_id": Column((12, 13)),
            "residue_num_begin": Column((14, 18), fields.Integer.from_pdb),
            "residue_icode_begin": Column((18, 19)),
            "residue_num_end": Column((20, 24), fields.Integer.from_pdb),
            "residue_icode_end": Column((24, 25)),
            "db": Column((26, 32), np.char.strip),
            "db_chain_id": Column((47, 67), np.char.strip),
        },
    )

    DBREF2 = (
        False,
        False,
        True,
        {
            # "pdb_id": Column((7, 11)),
            "chain_id": Column((12, 13)),
            "db_chain_accession": Column((18, 40), np.char.strip),
            "db_residue_num_begin": Column((45, 55), fields.Integer.from_pdb),
            "db_residue_num_end": Column((57, 67), fields.Integer.from_pdb),
        }
    )

    SEQADV = (
        False,
        False,
        True,
        {
            # "pdb_id": Column((7, 11)),
            "residue_name": Column((12, 15), np.char.strip),
            "chain_id": Column((16, 17)),
            "residue_num": Column((18, 22), fields.Integer.from_pdb),
            "residue_icode": Column((22, 23)),
            "database": Column((24, 28), np.char.strip),
            "dbAccession": Column((29, 38), np.char.strip),
            "dbRes": Column((39, 42), np.char.strip),
            "dbSeq": Column((43, 48), fields.Integer.from_pdb),
            "description": Column((49, 70), np.char.strip)
        },
        np.concatenate(([6, 11, 15, 17, 23, 28, 38, 42, 48], np.arange(70, 80)))
    )

    SEQRES = (
        False,
        False,
        False,
        {
            "serial": Column((7, 10), fields.Integer.from_pdb),
            "chain_id": Column((11, 12)),
            "residue_count": Column((13, 17), fields.Integer.from_pdb),
            "residue_name": MultiColumn([(c, c+3) for c in range(19, 68, 4)], np.char.strip)
        },
        np.concatenate(([6, 10, 12, 17, 18], np.arange(22, 69, 4), np.arange(70, 80))),
    )

    MODRES = (
        False,
        False,
        True,
        {
            # "pdb_id": Column((7, 11)),
            "residue_name": Column((12, 15), np.char.strip),
            "chain_id": Column((16, 17)),
            "residue_num": Column((18, 22), fields.Integer.from_pdb),
            "residue_icode": Column((22, 23)),
            "residue_name_std": Column((24, 27), np.char.strip),
            "description": Column((29, 70), np.char.strip)
        },
        np.concatenate(([6, 11, 15, 17, 23, 27, 28], np.arange(70, 80)))
    )

    HET = (
        False,
        False,
        True,
        {
            "heterogen_id": Column((7, 10), np.char.strip),
            "chain_id": Column((12, 13)),
            "residue_num": Column((13, 17), fields.Integer.from_pdb),
            "residue_icode": Column((17, 18)),
            "hetatm_count": Column((20, 25), fields.Integer.from_pdb),
            "description": Column((30, 70), np.char.strip)
        }
    )

    HETNAM = (
        False,
        False,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "id": Column((11, 14), np.char.strip),
            "name": Column((15, 70))
        }
    )

    HETSYN = (
        False,
        False,
        False,
        {
            "continuation": Column((8, 10), fields.Continuation.from_pdb),
            "id": Column((11, 14), np.char.strip),
            "synonyms": Column((15, 17))
        }
    )

    FORMUL = (
        False,
        False,
        False,
        {
            "component_num": Column((8, 10), fields.Integer.from_pdb),
            "id": Column((12, 15), np.char.strip),
            "continuation": Column((16, 18), fields.Continuation.from_pdb),
            "is_water": Column((18, 19)),
            "formula": Column((19, 70))
        }
    )

    HELIX = (
        False,
        False,
        True,
    )

    SHEET = (
        False,
        False,
        True,
    )

    SSBOND = (
        False,
        False,
        True,
    )

    LINK = (
        False,
        False,
        True,
    )

    CISPEP = (
        False,
        False,
        True,
    )

    SITE = (
        False,
        False,
        False,
    )

    CRYST1 = (
        True,
        True,
        True,
    )

    ORIGX1 = (
        True,
        True,
        True,
    )

    ORIGX2 = (
        True,
        True,
        True,
    )

    ORIGX3 = (
        True,
        True,
        True,
    )

    SCALE1 = (
        True,
        True,
        True,
    )

    SCALE2 = (
        True,
        True,
        True,
    )

    SCALE3 = (
        True,
        True,
        True,
    )

    MTRIX1 = (
        False,
        False,
        True,
    )

    MTRIX2 = (
        False,
        False,
        True,
    )

    MTRIX3 = (
        False,
        False,
        True,
    )

    MODEL = (
        False,
        False,
        True,
    )

    ATOM = (
        False,
        False,
        True,
        {
            "serial": Column((6, 11), fields.Integer.from_pdb),
            "atom_name": Column((12, 16), np.char.strip),
            "alt_location": Column((16, 17)),
            "residue_name": Column((17, 20), np.char.strip),
            "chain_id": Column((21, 22)),
            "residue_num": Column((22, 26), fields.Integer.from_pdb),
            "residue_icode": Column((26, 27)),
            "x": Column((30, 38), fields.Real.from_pdb),
            "y": Column((38, 46), fields.Real.from_pdb),
            "z": Column((46, 54), fields.Real.from_pdb),
            "occupancy": Column((54, 60), fields.Real.from_pdb),
            "temp_factor": Column((60, 66), fields.Real.from_pdb),
            "element": Column((76, 78), np.char.strip),
            "charge": Column((78, 80), fields.atom_charge)
        },
        np.concatenate(([11, 20, 27, 28, 29], np.arange(66, 77)))
    )

    ANISOU = (
        False,
        False,
        True,
        {
            "serial": Column((6, 11), fields.Integer.from_pdb),
            "atom_name": Column((12, 16), np.char.strip),
            "alt_location": Column((16, 17)),
            "residue_name": Column((17, 20), np.char.strip),
            "chain_id": Column((21, 22)),
            "residue_num": Column((22, 26), fields.Integer.from_pdb),
            "residue_icode": Column((26, 27)),
            "u11": Column((28, 35), fields.Integer.from_pdb),
            "u22": Column((35, 42), fields.Integer.from_pdb),
            "u33": Column((42, 49), fields.Integer.from_pdb),
            "u12": Column((49, 56), fields.Integer.from_pdb),
            "u13": Column((56, 63), fields.Integer.from_pdb),
            "u23": Column((63, 70), fields.Integer.from_pdb),
            "element": Column((76, 78), np.char.strip),
            "charge": Column((78, 80), fields.atom_charge)
        }
    )

    TER = (
        False,
        False,
        True,
        {
            "serial": Column((6, 11), fields.Integer.from_pdb),
            "residue_name": Column((17, 20), np.char.strip),
            "chain_id": Column((21, 22)),
            "residue_num": Column((22, 26), fields.Integer.from_pdb),
            "residue_icode": Column((26, 27)),
        }
    )

    HETATM = (
        False,
        False,
        True,
    )

    ENDMDL = (
        False,
        False,
        True,
    )

    CONECT = (
        False,
        False,
        True,
        {
            "serial_ref": ((7, 11), int),
            "serial_b1": ((12, 16), int),
            "serial_b2": ((17, 21), int),
            "serial_b3": ((22, 26), int),
            "serial_b4": ((27, 31), int),
        }
    )

    MASTER = (
        True,
        True,
        True,
        {
            "numRemark": Column((10, 15), fields.Integer.from_pdb),
            "numHet": Column((20, 25), fields.Integer.from_pdb),
            "numHelix": Column((25, 30), fields.Integer.from_pdb),
            "numSheet": Column((30, 35), fields.Integer.from_pdb),
            "numSite": Column((40, 45), fields.Integer.from_pdb),
            "numXform": Column((45, 50), fields.Integer.from_pdb),
            "numCoord": Column((50, 55), fields.Integer.from_pdb),
            "numTer": Column((55, 60), fields.Integer.from_pdb),
            "numConect": Column((60, 65), fields.Integer.from_pdb),
            "numSeq": Column((65, 70), fields.Integer.from_pdb),
        },
        np.concatenate((np.arange(6, 10), np.arange(35, 40), np.arange(70, 80)))
    )

    END = (
        True,
        True,
        True,
        dict(),
        np.arange(6, 80)
    )

    @property
    def is_mandatory(self) -> bool:
        """Whether the record is a mandatory record that must exist in every PDB file."""
        return self.value[0]

    @property
    def is_one_time(self) -> bool:
        """Whether the record only appears once in a PDB file."""
        return self.value[1]

    @property
    def is_single_line(self) -> bool:
        """Whether the record only takes a single line per appearance."""
        return self.value[2]

    @property
    def columns(self) -> dict:
        return self.value[3]

    @property
    def empty_columns(self) -> np.ndarray:
        return self.value[4]

    @classmethod
    @property
    def names(cls) -> List:
        return cls.__dict__['_member_names_']


