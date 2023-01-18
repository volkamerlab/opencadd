

from typing import Tuple, Sequence, Union, List, NamedTuple, Dict, Callable
import datetime
import functools


import numpy as np
import pandas as pd

from opencadd.io.parsing import extract_column
from opencadd.io.pdb import fields
from opencadd.io.parsing import Column, MultiColumn


class Record:

    def __init__(
            self,
            name: str,
            is_mandatory: bool,
            is_one_time: bool,
            is_single_line: bool,
            column_data: dict,
            inner_structure: dict = None,
    ):
        self._name = name
        self._is_mandatory = is_mandatory
        self._is_one_time = is_one_time
        self._is_single_line = is_single_line
        self._inner_structure = inner_structure

        self._empty_columns: np.ndarray
        pass

    @property
    def is_mandatory(self) -> bool:
        """Whether the record is a mandatory record that must exist in every PDB file."""
        return self._is_mandatory

    @property
    def is_one_time(self) -> bool:
        """Whether the record only appears once in a PDB file."""
        return self._is_one_time

    @property
    def is_single_line(self) -> bool:
        """Whether the record only takes a single line per appearance."""
        return self._is_single_line

    @property
    def columns(self) -> dict:
        return self.value[3]

    @property
    def empty_columns(self) -> np.ndarray:
        return self._empty_columns


class Records:
    def __init__(self):

        field_continuation = ("continuation", (8, 10), fields.Continuation)

        self.HEADER = Record(
            name="HEADER",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "classification": ("classification", (10, 50), fields.String),
                "deposition_date": ("depDate", (50, 59), fields.Date),
                "pdb_id": ("idCode", (62, 66), fields.IDcode),
            },
        )
        self.OBSLTE = Record(
            name="OBSLTE",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "replacement_date": ("repDate", (11, 20), fields.Date),
                "pdb_id": ("idCode", (21, 25), fields.IDcode),
                "replaced_pdb_ids": ("rIdCode", [(i, i+4) for i in range(31, 72, 5)], fields.IDcode)
            },
        )
        self.TITLE = Record(
            name="TITLE",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "title": ("title", (10, 80), fields.String)
            },
        )
        self.SPLIT = Record(
            name="SPLIT",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "pdb_ids": ("idCode", [(i, i+4) for i in range(11, 77, 5)], fields.IDcode)
            },
        )
        self.CAVEAT = Record(
            name="CAVEAT",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "pdb_id": ("idCode", (11, 15), fields.IDcode),
                "description": ("comment", (19, 79), fields.String),
            },
        )
        self.COMPND = Record(
            name="COMPND",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": ("continuation", (7, 10), fields.Continuation),
                "compound": ("compound", (10, 80), fields.SpecificationList)
            },
            inner_structure={
                "mol_id": ("MOL_ID", fields.Integer),
                "name": ("MOLECULE", fields.LString),
                "chain_ids": ("CHAIN", fields.List),
                "fragment": ("FRAGMENT", fields.List),
                "synonyms": ("SYNONYM", fields.List),
                "enzyme_commission_num": ("EC", fields.List),
                "engineered": ("ENGINEERED", fields.LString),
                "mutation": ("MUTATION", fields.LString),
                "description": ("OTHER_DETAILS", fields.LString)
            },
        )
        self.SOURCE = Record(
            name="SOURCE",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": ("continuation", (7, 10), fields.Continuation),
                "source_name": ("srcName", (10, 79), fields.SpecificationList),
            },
            inner_structure={
                "mol_id": ("MOL_ID", fields.Integer),
                "synthetic": ("SYNTHETIC", fields.LString),
                "fragment": ("FRAGMENT", fields.LString),
                "organism": ("ORGANISM_COMMON", fields.LString),
                "organism_sci": ("ORGANISM_SCIENTIFIC", fields.LString),
                "organism_tax_id": ("ORGANISM_TAXID", fields.LString),
                "strain": ("STRAIN", fields.LString),
                "variant": ("VARIANT", fields.LString),
                "cell_line": ("CELL_LINE", fields.LString),
                "atcc_id": ("ATCC", fields.LString),
                "organ": ("ORGAN", fields.LString),
                "tissue": ("TISSUE", fields.LString),
                "cell": ("CELL", fields.LString),
                "organelle": ("ORGANELLE", fields.LString),
                "secretion": ("SECRETION", fields.LString),
                "cell_loc": ("CELLULAR_LOCATION", fields.LString),
                "plasmid": ("PLASMID", fields.LString),
                "gene": ("GENE", fields.LString),
                "expsys": ("EXPRESSION_SYSTEM_COMMON", fields.LString),
                "expsys_sci": ("EXPRESSION_SYSTEM", fields.LString),
                "expsys_tax_id": ("EXPRESSION_SYSTEM_TAXID", fields.LString),
                "expsys_strain": ("EXPRESSION_SYSTEM_STRAIN", fields.LString),
                "expsys_variant": ("EXPRESSION_SYSTEM_VARIANT", fields.LString),
                "expsys_cell_line": ("EXPRESSION_SYSTEM_CELL_LINE", fields.LString),
                "expsys_atcc_id": ("EXPRESSION_SYSTEM_ATCC_NUMBER", fields.LString),
                "expsys_organ": ("EXPRESSION_SYSTEM_ORGAN", fields.LString),
                "expsys_tissue": ("EXPRESSION_SYSTEM_TISSUE", fields.LString),
                "expsys_cell": ("EXPRESSION_SYSTEM_CELL", fields.LString),
                "expsys_organelle": ("EXPRESSION_SYSTEM_ORGANELLE", fields.LString),
                "expsys_cell_loc": ("EXPRESSION_SYSTEM_CELLULAR_LOCATION", fields.LString),
                "expsys_vector_type": ("EXPRESSION_SYSTEM_VECTOR_TYPE", fields.LString),
                "expsys_vector": ("EXPRESSION_SYSTEM_VECTOR", fields.LString),
                "expsys_plasmid": ("EXPRESSION_SYSTEM_PLASMID", fields.LString),
                "expsys_gene": ("EXPRESSION_SYSTEM_GENE", fields.LString),
                "details": ("OTHER_DETAILS", fields.LString),
            },
        ),
        self.KEYWDS = Record(
            name="KEYWDS",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "keywords": ("keywds", (10, 79), fields.List)
            },
        )
        self.EXPDTA = Record(
            name="EXPDTA",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "technique": ("technique", (10, 79), fields.SList)
            },
        )
        self.NUMMDL = Record(
            name="NUMMDL",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "count_models": ("modelNumber", (10, 14), fields.Integer)
            }
        )
        self.MDLTYP = Record(
            name="MDLTYP",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "description": ("comment", (10, 80), fields.SList)
            },
        )
        self.AUTHOR = Record(
            name="AUTHOR",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "authors": ("authorList", (10, 79), fields.List)
            },
        )
        self.REVDAT = Record(
            name="REVDAT",
            is_mandatory=True,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "modification_num": ("modNum", (7, 10), fields.Integer),
                "continuation": ("continuation", (10, 12), fields.Continuation),
                "modification_date": ("modDate", (13, 22), fields.Date),
                "pdb_id": ("modId", (23, 27), fields.IDcode),
                "modification_type": ("modType", (31, 32), fields.Integer),
                "modification_details": ("record", [(i, i + 5) for i in range(39, 61, 7)], fields.LString)
            }
        )
        self.SPRSDE = Record(
            name="SPRSDE",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "superseded_date": ("sprsdeDate", (11, 20), fields.Date),
                "pdb_id": ("idCode", (21, 25), fields.IDcode),
                "superseded_pdb_ids": ("sIdCode", [(i, i + 4) for i in range(31, 72, 5)], fields.IDcode)
            },
        )
        self.JRNL = Record(
            name="JRNL",
            is_mandatory=False,
            is_one_time=True,
            is_single_line=False,
            column_data={
                "tag": ("tag", (12, 16), fields.LString),
                "continuation": ("continuation", (16, 18), fields.Continuation),
                "data": ("data", (19, 80), fields.LString),
            },
            inner_structure={
                "AUTH": {
                    "authors": ("authorList", (19, 79), fields.List)
                },
                "TITL": {
                    "title": ("title", (19, 79), fields.LString)
                },
            }
        )
        self.REMARK = Record(
            name="REMARK",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "remark_num": ("remarkNum", (7, 10), fields.Integer),
                "content": ((11, 79), fields.LString),
            }
        )
        self.DBREF = Record(
            name="DBREF",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "pdb_id": ("idCode", (7, 11), fields.IDcode),
                "chain_id": ("chainID", (12, 13), fields.Character),
                "residue_num_begin": ("seqBegin", (14, 18), fields.Integer),
                "residue_icode_begin": ("insertBegin", (18, 19), fields.AChar),
                "residue_num_end": ("seqEnd", (20, 24), fields.Integer),
                "residue_icode_end": ("insertEnd", (24, 25), fields.AChar),
                "db": ("database", (26, 32), fields.LString),
                "db_chain_accession": ("dbAccession", (33, 41), fields.LString),
                "db_chain_id": ("dbIdCode", (42, 54), fields.LString),
                "db_residue_num_begin": ("dbseqBegin", (55, 60), fields.Integer),
                "db_residue_icode_begin": ("idbnsBeg", (60, 61), fields.AChar),
                "db_residue_num_end": ("dbseqEnd", (62, 67), fields.Integer),
                "db_residue_icode_end": ("dbinsEnd", (67, 68), fields.AChar),
            },
        )
        self.DBREF1 = Record(
            name="DBREF1",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "pdb_id": ("idCode", (7, 11), fields.IDcode),
                "chain_id": ("chainID", (12, 13), fields.Character),
                "residue_num_begin": ("seqBegin", (14, 18), fields.Integer),
                "residue_icode_begin": ("insertBegin", (18, 19), fields.AChar),
                "residue_num_end": ("seqEnd", (20, 24), fields.Integer),
                "residue_icode_end": ("insertEnd", (24, 25), fields.AChar),
                "db": ("database", (26, 32), fields.LString),
                "db_chain_id": ("dbIdCode", (47, 67), fields.LString),
            },
        )
        self.DBREF2 = Record(
            name="DBREF2",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "pdb_id": ("idCode", (7, 11), fields.IDcode),
                "chain_id": ("chainID", (12, 13), fields.Character),
                "db_chain_accession": ("dbAccession", (18, 40), fields.LString),
                "db_residue_num_begin": ("seqBegin", (45, 55), fields.Integer),
                "db_residue_num_end": ("seqEnd", (57, 67), fields.Integer),
            }
        )
        self.SEQADV = Record(
            name="SEQADV",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "pdb_id": ("idCode", (7, 11), fields.IDcode),
                "residue_name": ("resName", (12, 15), fields.ResidueName),
                "chain_id": ("chainID", (16, 17), fields.Character),
                "residue_num": ("seqNum", (18, 22), fields.Integer),
                "residue_icode": ("iCode", (22, 23), fields.AChar),
                "db": ("database", (24, 28), fields.LString),
                "db_chain_accession": ("dbAccession", (29, 38), fields.LString),
                "db_residue_name": ("dbRes", (39, 42), fields.ResidueName),
                "db_residue_num": ("dbSeq", (43, 48), fields.Integer),
                "description": ("conflict", (49, 70), fields.LString)
            },
        )
        self.SEQRES = Record(
            name="SEQRES",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "serial": ("serNum", (7, 10), fields.Integer),
                "chain_id": ("chainID", (11, 12), fields.Character),
                "residue_count": ("numRes", (13, 17), fields.Integer),
                "residue_name": ("resName", [(c, c + 3) for c in range(19, 68, 4)], fields.ResidueName)
            },
        )
        self.MODRES = Record(
            name="MODRES",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "pdb_id": ("idCode", (7, 11), fields.IDcode),
                "residue_name": ("resName", (12, 15), fields.ResidueName),
                "chain_id": ("chainID", (16, 17), fields.Character),
                "residue_num": ("seqNum", (18, 22), fields.Integer),
                "residue_icode": ("iCode", (22, 23), fields.AChar),
                "residue_name_std": ("stdRes", (24, 27), fields.ResidueName),
                "description": ("comment", (29, 70), fields.LString)
            },
        )
        self.HET = Record(
            name="HET",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "het_id": ("hetID", (7, 10), fields.LString),
                "chain_id": ("chainID", (12, 13), fields.Character),
                "residue_num": ("seqNum", (13, 17), fields.Integer),
                "residue_icode": ("iCode", (17, 18), fields.AChar),
                "hetatm_count": ("numHetAtoms", (20, 25), fields.Integer),
                "description": ("text", (30, 70), fields.LString)
            }
        )
        self.HETNAM = Record(
            name="HETNAM",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "het_id": ("hetID", (11, 14), fields.LString),
                "name": ("text", (15, 70), fields.String)
            }
        )
        self.HETSYN = Record(
            name="HETSYN",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "continuation": field_continuation,
                "he_id": ("hetID", (11, 14), fields.LString),
                "synonyms": ("hetSynonyms", (15, 17), fields.SList)
            }
        )
        self.FORMUL = Record(
            name="FORMUL",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "component_num": ("compNum", (8, 10), fields.Integer),
                "het_id": ("hetID", (12, 15), fields.LString),
                "continuation": ("continuation", (16, 18), fields.Continuation),
                "is_water": ("asterisk", (18, 19), fields.Character),
                "formula": ("text", (19, 70), fields.String)
            }
        )
        self.HELIX = Record(
            name="HELIX",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serNum", (7, 10), fields.Integer),
                "helix_id": ("helixID", (11, 14), fields.LString),
                "residue_name_begin": ("initResName", (15, 18), fields.ResidueName),
                "chain_id_begin": ("initChainID", (19, 20), fields.Character),
                "residue_num_begin": ("initSeqNum", (21, 25), fields.Integer),
                "residue_icode_begin": ("initICode", (25, 26), fields.AChar),
                "residue_name_end": ("endResName", (27, 30), fields.ResidueName),
                "chain_id_end": ("endChainID", (31, 32), fields.Character),
                "residue_num_end": ("endSeqNum", (33, 37), fields.Integer),
                "residue_icode_end": ("endICode", (37, 38), fields.AChar),
                "class": ("helixClass", (38, 40), fields.Integer),
                "description": ("comment", (40, 70), fields.LString),
                "length": ("length", (71, 76), fields.Integer),
            }
        )
        self.SHEET = Record(
            name="SHEET",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "strand": ("strand", (7, 10), fields.Integer),
                "sheet_id": ("sheetID", (11, 14), fields.LString),
                "count_strands": ("numStrands", (14, 16), fields.Integer),
                "residue_name_begin": ("initResName", (17, 20), fields.ResidueName),
                "chain_id_begin": ("initChainID", (21, 22), fields.Character),
                "residue_num_begin": ("initSeqNum", (22, 26), fields.Integer),
                "residue_icode_begin": ("initICode", (26, 27), fields.AChar),
                "residue_name_end": ("endResName", (28, 31), fields.ResidueName),
                "chain_id_end": ("endChainID", (32, 33), fields.Character),
                "residue_num_end": ("endSeqNum", (33, 37), fields.Integer),
                "residue_icode_end": ("endICode", (37, 38), fields.AChar),
                "sense": ("sense", (38, 40), fields.Integer),
                "atom_curr": ("curAtom", (41, 45), fields.Atom),
                "residue_name_curr": ("curResName", (45, 48), fields.ResidueName),
                "chain_id_curr": ("curChainID", (49, 50), fields.Character),
                "residue_num_curr": ("curResSeq", (50, 54), fields.Integer),
                "residue_icode_curr": ("curICode", (54, 55), fields.AChar),
                "atom_prev": ("prevAtom", (56, 60), fields.Atom),
                "residue_name_prev": ("prevResName", (60, 63), fields.ResidueName),
                "chain_id_end_prev": ("prevChainID", (64, 65), fields.Character),
                "residue_num_end_prev": ("prevResSeq", (65, 69), fields.Integer),
                "residue_icode_end_prev": ("prevICode", (69, 70), fields.AChar)
            }
        )
        self.SSBOND = Record(
            name="SSBOND",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serNum", (7, 10), fields.Integer),
                "residue_name_1": ("CYS1", (11, 14), fields.LString),
                "chain_id_1": ("ChainID1", (15, 16), fields.Character),
                "residue_num_1": ("seqNum1", (17, 21), fields.Integer),
                "residue_icode_1": ("iCode1", (21, 22), fields.Atom),
                "residue_name_2": ("CYS2", (25, 28), fields.LString),
                "chain_id_2": ("chainID2", (29, 30), fields.Character),
                "residue_num_2": ("seqNum2", (31, 35), fields.Integer),
                "residue_icode_2": ("iCode2", (35, 36), fields.AChar),
                "symmetry_op_1": ("sym1", (59, 65), fields.SymOP),
                "symmetry_op_2": ("sym2", (66, 72), fields.SymOP),
                "bond_length": ("length", (73, 78), fields.Real),
            }
        )
        self.LINK = Record(
            name="LINK",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "atom_name_1": ("name1", (12, 16), fields.Atom),
                "alt_location_1": ("altLoc1", (16, 17), fields.Character),
                "residue_name_1": ("resName1", (17, 20), fields.ResidueName),
                "chain_id_1": ("chainID1", (21, 22), fields.Character),
                "residue_num_1": ("resSeq1", (22, 26), fields.Integer),
                "residue_icode_1": ("iCode1", (26, 27), fields.AChar),
                "atom_name_2": ("name2", (42, 46), fields.Atom),
                "alt_location_2": ("altLoc2", (46, 47), fields.Character),
                "residue_name_2": ("resName2", (47, 50), fields.ResidueName),
                "chain_id_2": ("chainID2", (51, 52), fields.Character),
                "residue_num_2": ("resSeq2", (52, 56), fields.Integer),
                "residue_icode_2": ("iCode2", (56, 57), fields.AChar),
                "symmetry_op_1": ("sym1", (59, 65), fields.SymOP),
                "symmetry_op_2": ("sym2", (66, 72), fields.SymOP),
                "distance": ("length", (73, 78), fields.Real)
            }
        )
        self.CISPEP = Record(
            name="CISPEP",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serNum", (7, 10), fields.Integer),
                "residue_name_1": ("pep1", (11, 14), fields.LString),
                "chain_id_1": ("chainID1", (15, 16), fields.Character),
                "residue_num_1": ("seqNum1", (17, 21), fields.Integer),
                "residue_icode_1": ("iCode1", (21, 22), fields.AChar),
                "residue_name_2": ("pep2", (25, 28), fields.LString),
                "chain_id_2": ("chainID2", (29, 30), fields.Character),
                "residue_num_2": ("seqNum2", (31, 35), fields.Integer),
                "residue_icode_2": ("iCode2", (35, 36), fields.AChar),
                "model_num": ("modNum", (43, 46), fields.Integer),
                "angle": ("measure", (53, 59), fields.Real)
            }
        )
        self.SITE = Record(
            name="SITE",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=False,
            column_data={
                "serial": ("seqNum", (7, 10), fields.Integer),
                "site_id": ("siteID", (11, 14), fields.LString),
                "residue_count": ("numRes", (15, 17), fields.Integer),
                "residue_name": ("resName", ((18, 21), (29, 32), (40, 43), (51, 54)), fields.ResidueName),
                "chain_id": ("chainID", ((22, 23), (33, 34), (44, 45), (55, 56)), fields.Character),
                "residue_num": ("seq", ((23, 27), (34, 38), (45, 49), (56, 60)), fields.Integer),
                "residue_icode": ("iCode", ((27, 28), (38, 39), (49, 50), (60, 61)), fields.AChar),
            }
        )
        self.CRYST1 = Record(
            name="CRYST1",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "a": ("a", (6, 15), fields.Real),
                "b": ("b", (15, 24), fields.Real),
                "c": ("c", (24, 33), fields.Real),
                "alpha": ("alpha", (33, 40), fields.Real),
                "beta": ("beta", (40, 47), fields.Real),
                "gamma": ("gamma", (47, 54), fields.Real),
                "space_group": ("space_group", (55, 66), fields.LString),
                "z": ("z", (66, 70), fields.Integer),
            }
        )
        self.ORIGX1 = Record(
            name="ORIGX1",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "o_1_1": ("o[1][1]", (10, 20), fields.Real),
                "o_1_2": ("o[1][2]", (20, 30), fields.Real),
                "o_1_3": ("o[1][3]", (30, 40), fields.Real),
                "t_1": ("t[1]", (45, 55), fields.Real),
            }
        )
        self.ORIGX2 = Record(
            name="ORIGX2",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "o_2_1": ("o[2][1]", (10, 20), fields.Real),
                "o_2_2": ("o[2][2]", (20, 30), fields.Real),
                "o_2_3": ("o[2][3]", (30, 40), fields.Real),
                "t_2": ("t[2]", (45, 55), fields.Real),
            }
        )
        self.ORIGX3 = Record(
            name="ORIGX3",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "o_3_1": ("o[3][1]", (10, 20), fields.Real),
                "o_3_2": ("o[3][2]", (20, 30), fields.Real),
                "o_3_3": ("o[3][3]", (30, 40), fields.Real),
                "t_3": ("t[3]", (45, 55), fields.Real),
            }
        )
        self.SCALE1 = Record(
            name="SCALE1",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "s_1_1": ("s[1][1]", (10, 20), fields.Real),
                "s_1_2": ("s[1][2]", (20, 30), fields.Real),
                "s_1_3": ("s[1][3]", (30, 40), fields.Real),
                "u_1": ("u[1]", (45, 55), fields.Real),
            }
        )
        self.SCALE2 = Record(
            name="SCALE2",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "s_2_1": ("s[2][1]", (10, 20), fields.Real),
                "s_2_2": ("s[2][2]", (20, 30), fields.Real),
                "s_2_3": ("s[2][3]", (30, 40), fields.Real),
                "u_2": ("u[2]", (45, 55), fields.Real),
            }
        )
        self.SCALE3 = Record(
            name="SCALE3",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "s_3_1": ("s[3][1]", (10, 20), fields.Real),
                "s_3_2": ("s[3][2]", (20, 30), fields.Real),
                "s_3_3": ("s[3][3]", (30, 40), fields.Real),
                "u_3": ("u[3]", (45, 55), fields.Real),
            }
        )
        self.MTRIX1 = Record(
            name="MTRIX1",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (7, 10), fields.Integer),
                "m_1_1": ("m[1][1]", (10, 20), fields.Real),
                "m_1_2": ("m[1][2]", (20, 30), fields.Real),
                "m_1_3": ("m[1][3]", (30, 40), fields.Real),
                "v_1": ("v[1]", (45, 55), fields.Real),
                "coordinates_are_contained": ("iGiven", (59, 60), fields.Integer)
            }
        )
        self.MTRIX2 = Record(
            name="MTRIX2",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (7, 10), fields.Integer),
                "m_2_1": ("m[2][1]", (10, 20), fields.Real),
                "m_2_2": ("m[2][2]", (20, 30), fields.Real),
                "m_2_3": ("m[2][3]", (30, 40), fields.Real),
                "v_2": ("v[2]", (45, 55), fields.Real),
                "coordinates_are_contained": ("iGiven", (59, 60), fields.Integer)
            }
        )
        self.MTRIX3 = Record(
            name="MTRIX3",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (7, 10), fields.Integer),
                "m_3_1": ("m[3][1]", (10, 20), fields.Real),
                "m_3_2": ("m[3][2]", (20, 30), fields.Real),
                "m_3_3": ("m[3][3]", (30, 40), fields.Real),
                "v_3": ("v[3]", (45, 55), fields.Real),
                "coordinates_are_contained": ("iGiven", (59, 60), fields.Integer)
            }
        )
        self.MODEL = Record(
            name="MODEL",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (10, 14), fields.Integer)
            }
        )
        self.ATOM = Record(
            name="ATOM",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (6, 11), fields.Integer),
                "atom_name": ("name", (12, 16), fields.Atom),
                "alt_location": ("altLoc", (16, 17), fields.Character),
                "residue_name": ("resName", (17, 20), fields.ResidueName),
                "chain_id": ("chainID", (21, 22), fields.Character),
                "residue_num": ("resSeq", (22, 26), fields.Integer),
                "residue_icode": ("iCode", (26, 27), fields.Character),
                "x": ("x", (30, 38), fields.Real),
                "y": ("y", (38, 46), fields.Real),
                "z": ("z", (46, 54), fields.Real),
                "occupancy": ("occupancy", (54, 60), fields.Real),
                "temp_factor": ("tempFactor", (60, 66), fields.Real),
                "element": ("element", (76, 78), fields.LString),
                "charge": ("charge", (78, 80), fields.LString)
            },
        )
        self.ANISOU = Record(
            name="ANISOU",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (6, 11), fields.Integer),
                "atom_name": ("name", (12, 16), fields.Atom),
                "alt_location": ("altLoc", (16, 17), fields.Character),
                "residue_name": ("resName", (17, 20), fields.ResidueName),
                "chain_id": ("chainID", (21, 22), fields.Character),
                "residue_num": ("resSeq", (22, 26), fields.Integer),
                "residue_icode": ("iCode", (26, 27), fields.AChar),
                "u_1_1": ("u[1][1]", (28, 35), fields.Integer),
                "u_2_2": ("u[2][2]", (35, 42), fields.Integer),
                "u_3_3": ("u[3][3]", (42, 49), fields.Integer),
                "u_1_2": ("u[1][2]", (49, 56), fields.Integer),
                "u_1_3": ("u[1][3]", (56, 63), fields.Integer),
                "u_2_3": ("u[2][3]", (63, 70), fields.Integer),
                "element": ("element", (76, 78), fields.LString),
                "charge": ("charge", (78, 80), fields.LString)
            }
        )
        self.TER = Record(
            name="TER",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (6, 11), fields.Integer),
                "residue_name": ("resName", (17, 20), fields.ResidueName),
                "chain_id": ("chainID", (21, 22), fields.Character),
                "residue_num": ("resSeq", (22, 26), fields.Integer),
                "residue_icode": ("iCode", (26, 27), fields.AChar),
            }
        )
        self.HETATM = Record(
            name="HETATM",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial": ("serial", (6, 11), fields.Integer),
                "atom_name": ("name", (12, 16), fields.Atom),
                "alt_location": ("altLoc", (16, 17), fields.Character),
                "residue_name": ("resName", (17, 20), fields.ResidueName),
                "chain_id": ("chainID", (21, 22), fields.Character),
                "residue_num": ("resSeq", (22, 26), fields.Integer),
                "residue_icode": ("iCode", (26, 27), fields.Character),
                "x": ("x", (30, 38), fields.Real),
                "y": ("y", (38, 46), fields.Real),
                "z": ("z", (46, 54), fields.Real),
                "occupancy": ("occupancy", (54, 60), fields.Real),
                "temp_factor": ("tempFactor", (60, 66), fields.Real),
                "element": ("element", (76, 78), fields.LString),
                "charge": ("charge", (78, 80), fields.LString)
            },
        )
        self.ENDMDL = Record(
            name="ENDMDL",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data=dict()
        )
        self.CONECT = Record(
            name="CONECT",
            is_mandatory=False,
            is_one_time=False,
            is_single_line=True,
            column_data={
                "serial_ref": ("serial1", (7, 11), fields.Integer),
                "serial_b1": ("serial2", (12, 16), fields.Integer),
                "serial_b2": ("serial3", (17, 21), fields.Integer),
                "serial_b3": ("serial4", (22, 26), fields.Integer),
                "serial_b4": ("serial5", (27, 31), fields.Integer),
            }
        )
        self.MASTER = Record(
            name="MASTER",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data={
                "remark": ("numRemark", (10, 15), fields.Integer),
                "0": ("0", (15, 20), fields.Integer),
                "het": ("numHet", (20, 25), fields.Integer),
                "helix": ("numHelix", (25, 30), fields.Integer),
                "sheet": ("numSheet", (30, 35), fields.Integer),
                "turn": ("numTurn", (35, 40), fields.Integer),
                "site": ("numSite", (40, 45), fields.Integer),
                "xform": ("numXform", (45, 50), fields.Integer),
                "coord": ("numCoord", (50, 55), fields.Integer),
                "ter": ("numTer", (55, 60), fields.Integer),
                "conect": ("numConect", (60, 65), fields.Integer),
                "seqres": ("numSeq", (65, 70), fields.Integer),
            },
        )
        self.END = Record(
            name="END",
            is_mandatory=True,
            is_one_time=True,
            is_single_line=True,
            column_data=dict(),
        )
        return

    @property
    def names(cls) -> List:
        return cls.__dict__['_member_names_']


def chain(*funcs: Callable):
    """
    Chain functions.

    Parameters
    ----------
    funcs : Callable
        Functions to chain.

    Returns
    -------
    Callable
        For n functions f_n, return F(x) = f_n(f_{n-1}(...(f_2(f_1(x))))
    """
    return functools.partial(lambda val: functools.reduce(lambda x, y: y(x), funcs, val))




"""Creates 1-letter amino acid codes from DataFrame

        Non-canonical amino-acids are converted as follows:
        ASH (protonated ASP) => D
        CYX (disulfide-bonded CYS) => C
        GLH (protonated GLU) => E
        HID/HIE/HIP (different protonation states of HIS) = H
        HYP (hydroxyproline) => P
        MSE (selenomethionine) => M
        """
amino3to1dict = {
    "ASH": "A",
    "ALA": "A",
    "CYX": "C",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PYL": "O",
    "HYP": "P",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "SEL": "U",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}









class Obsolete:
    def __init__(
            self,
            replacement_date: datetime.date,
            replaced_pdb_ids: Sequence[str],
    ):
        self.replacement_date: datetime.date = replacement_date
        self.replaced_pdb_ids: np.ndarray = np.asarray(replaced_pdb_ids)
        return

    def __repr__(self):
        return f"Obsolete({self.replacement_date}, {self.replaced_pdb_ids})"

    def __str__(self):
        return f"Replaced on {self.replacement_date} by {self.replaced_pdb_ids}."


class Split:
    pass


class Caveat:
    pass



class Source:
    pass


class Keywords:
    pass


class ExperimentalData:
    pass


class NumModels:
    pass


class ModelType:
    pass


class Author:
    pass


class RevisionData:
    pass


class SPRSDE:
    pass


class Journal:
    pass


class Remark2:
    pass


class Remark3:
    pass


class Cryst1:
    """
    CRYST1 record of the PDB file, containing the unit cell parameters, space group and Z value.

    Notes
    -----
    * If the structure was not determined by crystallographic means, CRYST1 contains the unitary values, i.e.:
        * a = b = c = 1.0
        * α = β = γ = 90 degrees
        * space group = P 1
        * Z = 1
    """
    def __init__(
            self,
            unit_cell_lengths: Union[Tuple[float, float, float], np.ndarray],
            unit_cell_angles: Union[Tuple[float, float, float], np.ndarray],
            space_group: str,
            z: int,
    ):
        self._unit_cell_lengths = np.asarray(unit_cell_lengths)
        self._unit_cell_angles = np.asarray(unit_cell_angles)
        self._space_group = space_group
        self._z = z
        return

    @property
    def length(self) -> np.ndarray:
        """
        Lattice parameters a, b, and c, i.e. the lengths of the unit cell, in Ångstrom.

        Returns
        -------
        numpy.ndarray, shape: (3,), dtype: float64
            array(a, b, c)
        """
        return self._unit_cell_lengths

    @property
    def angle(self) -> np.ndarray:
        """
        Lattice parameters α, β, and γ, i.e. the angles between the edges of the unit cell, in degrees.

        Returns
        -------
        numpy.ndarray, shape: (3,), dtype: float64
            array(α, β, γ)
        """
        return self._unit_cell_angles

    @property
    def space_group(self) -> str:
        """
        Hermnn-Mauguin space group symbol.

        Returns
        -------
        str
            Full international Table's Hermann-Mauguin symbol, without parenthesis, and the screw axis
            described as a two-digit number.
            Examples: 'P 1 21 1' (instead of 'P 21'), 'P 43 21 2'.

        Notes
        -----
        * For a rhombohedral space group in the hexagonal setting, the lattice type symbol used is H.
        """
        return self._space_group

    @property
    def z(self) -> int:
        """
        Z-value of the unit cell, i.e. the number of polymeric chains in a unit cell.
        In the case of heteropolymers, Z is the number of occurrences of the most populous chain.

        Returns
        -------
        int

        Notes
        -----
        As an example, given two chains A and B, each with a different sequence, and the space group 'P 2'
        that has two equipoints in the standard unit cell, the following table gives the correct Z value:
            Asymmetric Unit Content     Z value
            -----------------------     -------
                                  A     2
                                 AA     4
                                 AB     2
                                AAB     4
                               AABB     4
        """
        return self._z

    def __repr__(self):
        return f"Cryst1({tuple(self.length)}, {tuple(self.angle)}, {self.space_group}, {self.z})"

    def __str__(self):
        return (
            f"Space group: {self.space_group}\n"
            f"Z-values: {self.z}\n"
            f"Unit cell parameters:\n"
            f"\tLengths (a, b, c): {tuple(self.length)}\n"
            f"\tAngles (α, β, γ): {tuple(self.angle)}"
        )
