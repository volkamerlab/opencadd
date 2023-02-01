

from typing import Tuple, Sequence, Union, List, NamedTuple, Dict, Callable
import datetime
import functools


import numpy as np
import pandas as pd

#from opencadd.io.parsing import extract_column
from . import _fields as fields



class Record:

    def __init__(
            self,
            field_data: dict,
            index: int = None,
            name: str = None,
            is_mandatory: bool = False,
            is_one_time: bool = False,
            is_single_line: bool = False,
            inner_structure: dict = None,
    ):
        self._idx: int = index
        self._name = name
        self._is_mandatory = is_mandatory
        self._is_one_time = is_one_time
        self._is_single_line = is_single_line
        if inner_structure is not None:
            self._inner_structure_names = []
            self._inner_structure_mapping = dict()
            for name, (key, dtype) in inner_structure.items():
                self._inner_structure_names.append(name)
                self._inner_structure_mapping[key] = (name, dtype)
        else:
            self._inner_structure_names = None
            self._inner_structure_mapping = None

        self._fields = []
        self._field_names = []
        self._field_names_pdb = []
        self._fields_dict = dict()
        for field_name, field_config in field_data.items():
            self._field_names.append(field_name)
            field_name_orig, field_columns, field_dtype = field_config
            field_is_multicolumn = isinstance(field_columns[0], Sequence)
            field_class = fields.MultiColumn if field_is_multicolumn else fields.Column
            self._field_names_pdb.append(field_name_orig)
            field = field_class(field_columns, field_dtype)
            self._fields.append(field)
            self._fields_dict[field_name] = field
        self._empty_columns: np.ndarray
        pass

    def __eq__(self, other):
        return self.name == other.name

    @property
    def index(self) -> int:
        return self._idx

    @property
    def name(self) -> str:
        return self._name

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
    def is_one_time_single_line(self) -> bool:
        return self._is_one_time and self._is_single_line

    @property
    def is_one_time_multiple_lines(self):
        return self._is_one_time and (not self._is_single_line)

    @property
    def is_multiple_times_single_line(self) -> bool:
        return (not self._is_one_time) and self._is_single_line

    @property
    def fields(self) -> List:
        return self._fields

    @property
    def fields_dict(self) -> dict:
        return self._fields_dict

    @property
    def field_names(self):
        return self._field_names

    @property
    def field_names_pdb(self):
        return self._field_names_pdb

    @property
    def empty_columns(self) -> np.ndarray:
        return self._empty_columns

    def extract(self, char_table: np.ndarray):
        pass

    @property
    def inner_structure_names(self):
        return self._inner_structure_names

    @property
    def inner_structure_mapping(self):
        return self._inner_structure_mapping


HEADER = Record(
    index=0,
    name="HEADER",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "classification": ("classification", (10, 50), fields.LString),
        "deposition_date": ("depDate", (50, 59), fields.Date),
        "pdb_id": ("idCode", (62, 66), fields.IDcode),
    },
)


OBSLTE = Record(
    index=1,
    name="OBSLTE",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "replacement_date": ("repDate", (11, 20), fields.Date),
        "pdb_id": ("idCode", (21, 25), fields.IDcode),
        "replacement_pdb_ids": ("rIdCode", [(i, i+4) for i in range(31, 72, 5)], fields.IDcode)
    },
)


TITLE = Record(
    index=2,
    name="TITLE",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "title": ("title", (10, 80), fields.String)
    },
)


SPLIT = Record(
    index=3,
    name="SPLIT",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "pdb_ids": ("idCode", [(i, i+4) for i in range(11, 77, 5)], fields.IDcode)
    },
)


CAVEAT = Record(
    index=4,
    name="CAVEAT",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "pdb_id": ("idCode", (11, 15), fields.IDcode),
        "description": ("comment", (19, 79), fields.String),
    },
)


COMPND = Record(
    index=5,
    name="COMPND",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
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


SOURCE = Record(
    index=6,
    name="SOURCE",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
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
)


KEYWDS = Record(
    index=7,
    name="KEYWDS",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "keywords": ("keywds", (10, 79), fields.List)
    },
)


EXPDTA = Record(
    index=8,
    name="EXPDTA",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "technique": ("technique", (10, 79), fields.SList)
    },
)


NUMMDL = Record(
    index=9,
    name="NUMMDL",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "count_models": ("modelNumber", (10, 14), fields.Integer)
    }
)


MDLTYP = Record(
    index=10,
    name="MDLTYP",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "description": ("comment", (10, 80), fields.SList)
    },
)


AUTHOR = Record(
    index=11,
    name="AUTHOR",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "authors": ("authorList", (10, 79), fields.List)
    },
)


REVDAT = Record(
    index=12,
    name="REVDAT",
    is_mandatory=True,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "mod_num": ("modNum", (7, 10), fields.Integer),
        "continuation": ("continuation", (10, 12), fields.Continuation),
        "date": ("modDate", (13, 22), fields.Date),
        "pdb_id": ("modId", (23, 27), fields.IDcode),
        "is_initial": ("modType", (31, 32), fields.Integer),
        "details": ("record", [(i, i + 6) for i in range(39, 61, 7)], fields.LString)
    }
)


SPRSDE = Record(
    index=13,
    name="SPRSDE",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "date": ("sprsdeDate", (11, 20), fields.Date),
        "pdb_id": ("idCode", (21, 25), fields.IDcode),
        "superseded_pdb_ids": ("sIdCode", [(i, i + 4) for i in range(31, 72, 5)], fields.IDcode)
    },
)


JRNL = Record(
    index=14,
    name="JRNL",
    is_mandatory=False,
    is_one_time=True,
    is_single_line=False,
    field_data={
        "tag": ("tag", (12, 16), fields.LString),
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "data": ("data", (19, 80), fields.LString),
    },
    # inner_structure={
    #     "AUTH": {"authors": ("authorList", (19, 79), fields.List)},
    #     "TITL": {"title": ("title", (19, 79), fields.LString)},
    # }
)


JRNL_AUTH = Record(
    name="AUTH",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "authors": ("authorList", (19, 79), fields.List)
    }
)

JRNL_TITL = Record(
    name="TITL",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "title": ("title", (19, 79), fields.String)
    }
)

JRNL_EDIT = Record(
    name="EDIT",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "editors": ("editorList", (19, 79), fields.List)
    }
)

JRNL_REF = Record(
    name="REF",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "publication_name": ("pubName", (19, 47), fields.String),
        "volume": ("volume", (51, 55), fields.LString),
        "page": ("page", (56, 61), fields.LString),
        "year": ("year", (62, 66), fields.Integer),
    }
)

JRNL_PUBL = Record(
    name="PUBL",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "publisher": ("pub", (19, 70), fields.String)
    }
)

JRNL_REFN = Record(
    name="REFN",
    is_one_time=True,
    is_single_line=True,
    field_data={
        "issn_essn": ("ISSN/ESSN", (35, 39), fields.List),
        "issn": ("issn", (40, 65), fields.LString)
    }
)

JRNL_PMID = Record(
    name="PMID",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "pubmed_id": ("PMID", (19, 79), fields.List)
    }
)

JRNL_DOI = Record(
    name="DOI",
    is_one_time=True,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "doi": ("DOI", (19, 79), fields.List)
    }
)


REMARK = Record(
    index=15,
    name="REMARK",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "remark_num": ("remarkNum", (7, 10), fields.Integer),
        "content": ("content", (11, 80), fields.LString),
    }
)


DBREF = Record(
    index=16,
    name="DBREF",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


DBREF1 = Record(
    index=17,
    name="DBREF1",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


DBREF2 = Record(
    index=18,
    name="DBREF2",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "pdb_id": ("idCode", (7, 11), fields.IDcode),
        "chain_id": ("chainID", (12, 13), fields.Character),
        "db_chain_accession": ("dbAccession", (18, 40), fields.LString),
        "db_residue_num_begin": ("dbseqBegin", (45, 55), fields.Integer),
        "db_residue_num_end": ("dbseqEnd", (57, 67), fields.Integer),
    }
)


SEQADV = Record(
    index=19,
    name="SEQADV",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


SEQRES = Record(
    index=20,
    name="SEQRES",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "serial": ("serNum", (7, 10), fields.Integer),
        "chain_id": ("chainID", (11, 12), fields.Character),
        "residue_count": ("numRes", (13, 17), fields.Integer),
        "residue_name": ("resName", [(c, c + 3) for c in range(19, 68, 4)], fields.ResidueName)
    },
)


MODRES = Record(
    index=21,
    name="MODRES",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "pdb_id": ("idCode", (7, 11), fields.IDcode),
        "residue_name": ("resName", (12, 15), fields.ResidueName),
        "chain_id": ("chainID", (16, 17), fields.Character),
        "residue_num": ("seqNum", (18, 22), fields.Integer),
        "residue_icode": ("iCode", (22, 23), fields.AChar),
        "residue_name_std": ("stdRes", (24, 27), fields.ResidueName),
        "description": ("comment", (29, 70), fields.LString)
    },
)


HET = Record(
    index=22,
    name="HET",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "het_id": ("hetID", (7, 10), fields.LString),
        "chain_id": ("chainID", (12, 13), fields.Character),
        "residue_num": ("seqNum", (13, 17), fields.Integer),
        "residue_icode": ("iCode", (17, 18), fields.AChar),
        "hetatm_count": ("numHetAtoms", (20, 25), fields.Integer),
        "description": ("text", (30, 70), fields.LString)
    }
)


HETNAM = Record(
    index=23,
    name="HETNAM",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "het_id": ("hetID", (11, 14), fields.LString),
        "name": ("text", (15, 70), fields.String)
    }
)


HETSYN = Record(
    index=24,
    name="HETSYN",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "continuation": ("continuation", (8, 10), fields.Continuation),
        "het_id": ("hetID", (11, 14), fields.LString),
        "synonyms": ("hetSynonyms", (15, 70), fields.SList)
    }
)


FORMUL = Record(
    index=25,
    name="FORMUL",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "component_num": ("compNum", (8, 10), fields.Integer),
        "het_id": ("hetID", (12, 15), fields.LString),
        "continuation": ("continuation", (16, 18), fields.Continuation),
        "is_water": ("asterisk", (18, 19), fields.Character),
        "formula": ("text", (19, 70), fields.String)
    }
)


HELIX = Record(
    index=26,
    name="HELIX",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


SHEET = Record(
    index=27,
    name="SHEET",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


SSBOND = Record(
    index=28,
    name="SSBOND",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serNum", (7, 10), fields.Integer),
        "residue_name_1": ("CYS1", (11, 14), fields.LString),
        "chain_id_1": ("ChainID1", (15, 16), fields.Character),
        "residue_num_1": ("seqNum1", (17, 21), fields.Integer),
        "residue_icode_1": ("iCode1", (21, 22), fields.Atom),
        "residue_name_2": ("CYS2", (25, 28), fields.LString),
        "chain_id_2": ("chainID2", (29, 30), fields.Character),
        "residue_num_2": ("seqNum2", (31, 35), fields.Integer),
        "residue_icode_2": ("iCode2", (35, 36), fields.AChar),
        "symmetry_op_num_1": ("sym1", (59, 62), fields.Integer),
        "translation_vector_1": ("sym1", [(i, i+1) for i in range(62, 65)], fields.SymOP),
        "symmetry_op_num_2": ("sym2", (66, 69), fields.Integer),
        "translation_vector_2": ("sym2", [(i, i+1) for i in range(69, 72)], fields.SymOP),
        "bond_length": ("length", (73, 78), fields.Real),
    }
)


LINK = Record(
    index=29,
    name="LINK",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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
        "symmetry_op_num_1": ("sym1", (59, 62), fields.Integer),
        "translation_vector_1": ("sym1", [(i, i+1) for i in range(62, 65)], fields.SymOP),
        "symmetry_op_num_2": ("sym2", (66, 69), fields.Integer),
        "translation_vector_2": ("sym2", [(i, i+1) for i in range(69, 72)], fields.SymOP),
        "distance": ("length", (73, 78), fields.Real)
    }
)


CISPEP = Record(
    index=30,
    name="CISPEP",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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


SITE = Record(
    index=31,
    name="SITE",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=False,
    field_data={
        "serial": ("seqNum", (7, 10), fields.Integer),
        "site_id": ("siteID", (11, 14), fields.LString),
        "residue_count": ("numRes", (15, 17), fields.Integer),
        "residue_name": ("resName", ((18, 21), (29, 32), (40, 43), (51, 54)), fields.ResidueName),
        "chain_id": ("chainID", ((22, 23), (33, 34), (44, 45), (55, 56)), fields.Character),
        "residue_num": ("seq", ((23, 27), (34, 38), (45, 49), (56, 60)), fields.Integer),
        "residue_icode": ("iCode", ((27, 28), (38, 39), (49, 50), (60, 61)), fields.AChar),
    }
)


CRYST1 = Record(
    index=32,
    name="CRYST1",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "a": ("a", (6, 15), fields.Real),
        "b": ("b", (15, 24), fields.Real),
        "c": ("c", (24, 33), fields.Real),
        "alpha": ("alpha", (33, 40), fields.Real),
        "beta": ("beta", (40, 47), fields.Real),
        "gamma": ("gamma", (47, 54), fields.Real),
        "space_group": ("sGroup", (55, 66), fields.LString),
        "z": ("z", (66, 70), fields.Integer),
    }
)


ORIGX1 = Record(
    index=33,
    name="ORIGX1",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "o_1_1": ("o[1][1]", (10, 20), fields.Real),
        "o_1_2": ("o[1][2]", (20, 30), fields.Real),
        "o_1_3": ("o[1][3]", (30, 40), fields.Real),
        "t_1": ("t[1]", (45, 55), fields.Real),
    }
)


ORIGX2 = Record(
    index=34,
    name="ORIGX2",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "o_2_1": ("o[2][1]", (10, 20), fields.Real),
        "o_2_2": ("o[2][2]", (20, 30), fields.Real),
        "o_2_3": ("o[2][3]", (30, 40), fields.Real),
        "t_2": ("t[2]", (45, 55), fields.Real),
    }
)


ORIGX3 = Record(
    index=35,
    name="ORIGX3",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "o_3_1": ("o[3][1]", (10, 20), fields.Real),
        "o_3_2": ("o[3][2]", (20, 30), fields.Real),
        "o_3_3": ("o[3][3]", (30, 40), fields.Real),
        "t_3": ("t[3]", (45, 55), fields.Real),
    }
)


SCALE1 = Record(
    index=36,
    name="SCALE1",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "s_1_1": ("s[1][1]", (10, 20), fields.Real),
        "s_1_2": ("s[1][2]", (20, 30), fields.Real),
        "s_1_3": ("s[1][3]", (30, 40), fields.Real),
        "u_1": ("u[1]", (45, 55), fields.Real),
    }
)


SCALE2 = Record(
    index=37,
    name="SCALE2",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "s_2_1": ("s[2][1]", (10, 20), fields.Real),
        "s_2_2": ("s[2][2]", (20, 30), fields.Real),
        "s_2_3": ("s[2][3]", (30, 40), fields.Real),
        "u_2": ("u[2]", (45, 55), fields.Real),
    }
)


SCALE3 = Record(
    index=38,
    name="SCALE3",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
        "s_3_1": ("s[3][1]", (10, 20), fields.Real),
        "s_3_2": ("s[3][2]", (20, 30), fields.Real),
        "s_3_3": ("s[3][3]", (30, 40), fields.Real),
        "u_3": ("u[3]", (45, 55), fields.Real),
    }
)


MTRIX1 = Record(
    index=39,
    name="MTRIX1",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (7, 10), fields.Integer),
        "m_1_1": ("m[1][1]", (10, 20), fields.Real),
        "m_1_2": ("m[1][2]", (20, 30), fields.Real),
        "m_1_3": ("m[1][3]", (30, 40), fields.Real),
        "v_1": ("v[1]", (45, 55), fields.Real),
        "is_given": ("iGiven", (59, 60), fields.Character)
    }
)


MTRIX2 = Record(
    index=40,
    name="MTRIX2",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (7, 10), fields.Integer),
        "m_2_1": ("m[2][1]", (10, 20), fields.Real),
        "m_2_2": ("m[2][2]", (20, 30), fields.Real),
        "m_2_3": ("m[2][3]", (30, 40), fields.Real),
        "v_2": ("v[2]", (45, 55), fields.Real),
        "is_given": ("iGiven", (59, 60), fields.Character)
    }
)


MTRIX3 = Record(
    index=41,
    name="MTRIX3",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (7, 10), fields.Integer),
        "m_3_1": ("m[3][1]", (10, 20), fields.Real),
        "m_3_2": ("m[3][2]", (20, 30), fields.Real),
        "m_3_3": ("m[3][3]", (30, 40), fields.Real),
        "v_3": ("v[3]", (45, 55), fields.Real),
        "is_given": ("iGiven", (59, 60), fields.Character)
    }
)


MODEL = Record(
    index=42,
    name="MODEL",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (10, 14), fields.Integer)
    }
)


ATOM = Record(
    index=43,
    name="ATOM",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (6, 11), fields.Integer),
        "atom_name": ("name", (12, 16), fields.Atom),
        "alt_location": ("altLoc", (16, 17), fields.Character),
        "residue_name": ("resName", (17, 20), fields.ResidueName),
        "chain_id": ("chainID", (21, 22), fields.Character),
        "residue_num": ("resSeq", (22, 26), fields.Integer),
        "residue_icode": ("iCode", (26, 27), fields.AChar),
        "x": ("x", (30, 38), fields.Real),
        "y": ("y", (38, 46), fields.Real),
        "z": ("z", (46, 54), fields.Real),
        "occupancy": ("occupancy", (54, 60), fields.Real),
        "temp_factor": ("tempFactor", (60, 66), fields.Real),
        "element": ("element", (76, 78), fields.Element),
        "charge": ("charge", ((79, 80), (78, 79)), fields.LString)
    },
)


ANISOU = Record(
    index=44,
    name="ANISOU",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
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
        "element": ("element", (76, 78), fields.Element),
        "charge": ("charge", ((79, 80), (78, 79)), fields.LString)
    }
)


TER = Record(
    index=45,
    name="TER",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (6, 11), fields.Integer),
        "residue_name": ("resName", (17, 20), fields.ResidueName),
        "chain_id": ("chainID", (21, 22), fields.Character),
        "residue_num": ("resSeq", (22, 26), fields.Integer),
        "residue_icode": ("iCode", (26, 27), fields.AChar),
    }
)


HETATM = Record(
    index=46,
    name="HETATM",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial": ("serial", (6, 11), fields.Integer),
        "atom_name": ("name", (12, 16), fields.Atom),
        "alt_location": ("altLoc", (16, 17), fields.Character),
        "residue_name": ("resName", (17, 20), fields.ResidueName),
        "chain_id": ("chainID", (21, 22), fields.Character),
        "residue_num": ("resSeq", (22, 26), fields.Integer),
        "residue_icode": ("iCode", (26, 27), fields.AChar),
        "x": ("x", (30, 38), fields.Real),
        "y": ("y", (38, 46), fields.Real),
        "z": ("z", (46, 54), fields.Real),
        "occupancy": ("occupancy", (54, 60), fields.Real),
        "temp_factor": ("tempFactor", (60, 66), fields.Real),
        "element": ("element", (76, 78), fields.Element),
        "charge": ("charge", (78, 80), fields.LString)
    },
)


ENDMDL = Record(
    index=47,
    name="ENDMDL",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data=dict()
)


CONECT = Record(
    index=48,
    name="CONECT",
    is_mandatory=False,
    is_one_time=False,
    is_single_line=True,
    field_data={
        "serial_ref": ("serial", (7, 11), fields.Integer),
        "serial_bonded": ("serial", [(i, i+4) for i in range(12, 28, 5)], fields.Integer),
    }
)


MASTER = Record(
    index=49,
    name="MASTER",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data={
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


END = Record(
    index=50,
    name="END",
    is_mandatory=True,
    is_one_time=True,
    is_single_line=True,
    field_data=dict(),
)


records: Sequence[Record] = (
    HEADER,
    OBSLTE,
    TITLE,
    SPLIT,
    CAVEAT,
    COMPND,
    SOURCE,
    KEYWDS,
    EXPDTA,
    NUMMDL,
    MDLTYP,
    AUTHOR,
    REVDAT,
    SPRSDE,
    JRNL,
    REMARK,
    DBREF,
    DBREF1,
    DBREF2,
    SEQADV,
    SEQRES,
    MODRES,
    HET,
    HETNAM,
    HETSYN,
    FORMUL,
    HELIX,
    SHEET,
    SSBOND,
    LINK,
    CISPEP,
    SITE,
    CRYST1,
    ORIGX1,
    ORIGX2,
    ORIGX3,
    SCALE1,
    SCALE2,
    SCALE3,
    MTRIX1,
    MTRIX2,
    MTRIX3,
    MODEL,
    ATOM,
    ANISOU,
    TER,
    HETATM,
    ENDMDL,
    CONECT,
    MASTER,
    END,
)


names: np.ndarray = np.array([rec.name for rec in records])

count = names.size

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
