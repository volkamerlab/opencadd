"""

"""

from typing import Tuple, Sequence, Union, List, NamedTuple, Dict, Callable
import datetime
import functools

import numpy as np
import pandas as pd

from . import _fields as fields


class Record:

    def __init__(
            self,
            field_data: dict,
            index: int = None,
            name: str = None,
            is_mandatory: bool = False,
            is_one_time: bool = False,
            is_one_line: bool = False,
            key_unique: str = None,
            key_continuation: str = "continuation",
            keys_repeat: Sequence[str] = None,
            inner_structure: dict = None,
    ):
        self._idx: int = index
        self._name = name
        self._is_mandatory = is_mandatory
        self._is_one_time = is_one_time
        self._is_one_line = is_one_line
        self._key_unique = key_unique
        self._key_continuation = key_continuation
        self._keys_repeat = keys_repeat
        if inner_structure is not None:
            self._inner_structure_names = []
            self._inner_structure_mapping = dict()
            for name, (key, dtype) in inner_structure.items():
                self._inner_structure_names.append(name)
                self._inner_structure_mapping[key] = (name, dtype)
        else:
            self._inner_structure_names = None
            self._inner_structure_mapping = None

        if self._is_one_time:
            if self._is_one_line:
                self._frequency = "1t1l"
                self._extract = self._extract_record_one_time_one_line
            else:
                self._frequency = "1tml"
                self._extract = self._extract_record_one_time_multiple_lines
        else:
            if self._is_one_line:
                self._frequency = "mt1l"
                self._extract = self._extract_record_multiple_times_one_line
            else:
                self._frequency = "mtml"
                self._extract = self._extract_record_multiple_times_multiple_lines

        self._fields_dict = field_data
        all_column_indices = np.array([])
        for field in field_data.values():
            all_column_indices = np.concatenate((all_column_indices, field.indices))
        self._empty_columns: np.ndarray = np.setdiff1d(np.arange(80), all_column_indices)
        self._to_dataframe = self._name not in (
            "REMARK", "ATOM", "HETATM", "ANISOU", "MTRIX1", "MTRIX2", "MTRIX3", "CONECT"
        )
        return

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
        return self._is_one_line

    @property
    def key_unique(self):
        return self._key_unique

    @property
    def key_continuation(self):
        return self._key_continuation

    @property
    def keys_repeat(self):
        return self._keys_repeat

    @property
    def fields_dict(self) -> dict:
        return self._fields_dict

    @property
    def empty_columns(self) -> np.ndarray:
        return self._empty_columns

    def extract(self, record_char_table: np.ndarray):
        return self._extract(record_char_table)

    @property
    def inner_structure_names(self):
        return self._inner_structure_names

    @property
    def inner_structure_mapping(self):
        return self._inner_structure_mapping

    def _extract_record_one_time_one_line(self, record_lines: np.ndarray):
        data = {
            field_name: field.extract(char_table=record_lines[:1])[0]  # Take only the first occurrence
            for field_name, field in self.fields_dict.items()
        }
        keys = list(data.keys())
        if len(keys) == 1:
            return data[keys[0]]
        return data

    def _extract_record_one_time_multiple_lines(self, record_lines: np.ndarray):
        continuation = self.fields_dict["continuation"].extract(char_table=record_lines)
        line_idx_in_order = np.argsort(continuation, kind="stable")
        record_lines_ordered = record_lines[line_idx_in_order]
        data = {
            field_name: field.extract(char_table=record_lines_ordered)
            for field_name, field in self.fields_dict.items() if field_name != "continuation"
        }
        keys = list(data.keys())
        if len(keys) == 1:
            return data[keys[0]]
        return data

    def _extract_record_multiple_times_one_line(self, record_lines: np.ndarray):
        data = {
            field_name: field.extract(char_table=record_lines)
            for field_name, field in self.fields_dict.items()
        }
        return pd.DataFrame(data) if self._to_dataframe else data

    def _extract_record_multiple_times_multiple_lines(self, record_lines: np.ndarray):
        data = {
            field_name: field.extract(char_table=record_lines)
            for field_name, field in self.fields_dict.items()
        }
        field_uniques = data.pop(self.key_unique)
        field_continuation = data.pop(self.key_continuation)
        _, idx_unique_fields = np.unique(field_uniques, return_index=True)
        idx_unique_fields.sort()
        unique_keys = field_uniques[idx_unique_fields]
        rows = []
        keys_repeat = tuple() if self.keys_repeat is None else self.keys_repeat
        for unique_key in unique_keys:
            # correlated lines are merged to create a multiple-times/single-line dataframe format:
            mask = field_uniques == unique_key
            idx_lines = np.argwhere(mask).reshape(-1)
            continuation = field_continuation[mask]
            idx_mapping = np.argsort(continuation)
            idx_lines_sorted = idx_lines[idx_mapping]
            row = {self.key_unique: unique_key}
            for field_name, field_vals in data.items():
                fields_sorted = field_vals[idx_lines_sorted]
                non_empty_fields = fields_sorted[fields_sorted != ""]
                if non_empty_fields.size == 0:
                    row[field_name] = ""
                else:
                    field_val_casted = self.fields_dict[field_name].cast_to_dtype(non_empty_fields)
                    row[field_name] = field_val_casted[0 if field_name in keys_repeat else slice(None)]
            rows.append(row)
        return pd.DataFrame(rows).set_index(self.key_unique, drop=False)


HEADER = Record(
    index=0,
    name="HEADER",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data={
        "classification": fields.Column((10, 50), fields.LString),
        "dep_date": fields.Column((50, 59), fields.Date),
        "pdb_id": fields.Column((62, 66), fields.IDcode),
    },
)


OBSLTE = Record(
    index=1,
    name="OBSLTE",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "rep_date": fields.Column((11, 20), fields.Date, only_first=True),
        "pdb_id": fields.Column((21, 25), fields.IDcode, only_first=True),
        "rep_pdb_id": fields.Column([(i, i+4) for i in range(31, 72, 5)], fields.IDcode, only_non_empty=True)
    },
)


TITLE = Record(
    index=2,
    name="TITLE",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "title": fields.Column((10, 80), fields.String, strip=False)
    },
)


SPLIT = Record(
    index=3,
    name="SPLIT",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "split_pdb_id": fields.Column([(i, i+4) for i in range(11, 77, 5)], fields.IDcode, only_non_empty=True)
    },
)


CAVEAT = Record(
    index=4,
    name="CAVEAT",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "pdb_id": fields.Column((11, 15), fields.IDcode, only_first=True),
        "description": fields.Column((19, 79), fields.String, strip=False),
    },
)


COMPND = Record(
    index=5,
    name="COMPND",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((7, 10), fields.Continuation),
        "compound": fields.Column((10, 80), fields.SpecificationList, strip=False)
    },
    inner_structure={
        "mol_id": ("MOL_ID", fields.Integer),
        "molecule": ("MOLECULE", fields.LString),
        "chain_id": ("CHAIN", fields.List),
        "fragment": ("FRAGMENT", fields.List),
        "synonym": ("SYNONYM", fields.List),
        "ec": ("EC", fields.List),
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
    is_one_line=False,
    field_data={
        "continuation": fields.Column((7, 10), fields.Continuation),
        "src_name": fields.Column((10, 79), fields.SpecificationList, strip=False),
    },
    inner_structure={
        "mol_id": ("MOL_ID", fields.Integer),
        "synthetic": ("SYNTHETIC", fields.LString),
        "fragment": ("FRAGMENT", fields.LString),
        "organism_com": ("ORGANISM_COMMON", fields.LString),
        "organism_sci": ("ORGANISM_SCIENTIFIC", fields.LString),
        "organism_tax_id": ("ORGANISM_TAXID", fields.LString),
        "strain": ("STRAIN", fields.LString),
        "variant": ("VARIANT", fields.LString),
        "cell_line": ("CELL_LINE", fields.LString),
        "atcc": ("ATCC", fields.LString),
        "organ": ("ORGAN", fields.LString),
        "tissue": ("TISSUE", fields.LString),
        "cell": ("CELL", fields.LString),
        "organelle": ("ORGANELLE", fields.LString),
        "secretion": ("SECRETION", fields.LString),
        "cell_loc": ("CELLULAR_LOCATION", fields.LString),
        "plasmid": ("PLASMID", fields.LString),
        "gene": ("GENE", fields.LString),
        "exp_sys_com": ("EXPRESSION_SYSTEM_COMMON", fields.LString),
        "exp_sys_sci": ("EXPRESSION_SYSTEM", fields.LString),
        "exp_sys_tax_id": ("EXPRESSION_SYSTEM_TAXID", fields.LString),
        "exp_sys_strain": ("EXPRESSION_SYSTEM_STRAIN", fields.LString),
        "exp_sys_variant": ("EXPRESSION_SYSTEM_VARIANT", fields.LString),
        "exp_sys_cell_line": ("EXPRESSION_SYSTEM_CELL_LINE", fields.LString),
        "exp_sys_atcc": ("EXPRESSION_SYSTEM_ATCC_NUMBER", fields.LString),
        "exp_sys_organ": ("EXPRESSION_SYSTEM_ORGAN", fields.LString),
        "exp_sys_tissue": ("EXPRESSION_SYSTEM_TISSUE", fields.LString),
        "exp_sys_cell": ("EXPRESSION_SYSTEM_CELL", fields.LString),
        "exp_sys_organelle": ("EXPRESSION_SYSTEM_ORGANELLE", fields.LString),
        "exp_sys_cell_loc": ("EXPRESSION_SYSTEM_CELLULAR_LOCATION", fields.LString),
        "exp_sys_vector_type": ("EXPRESSION_SYSTEM_VECTOR_TYPE", fields.LString),
        "exp_sys_vector": ("EXPRESSION_SYSTEM_VECTOR", fields.LString),
        "exp_sys_plasmid": ("EXPRESSION_SYSTEM_PLASMID", fields.LString),
        "exp_sys_gene": ("EXPRESSION_SYSTEM_GENE", fields.LString),
        "description": ("OTHER_DETAILS", fields.LString),
    },
)


KEYWDS = Record(
    index=7,
    name="KEYWDS",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "keywds": fields.Column((10, 79), fields.List, strip=False)
    },
)


EXPDTA = Record(
    index=8,
    name="EXPDTA",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "technique": fields.Column((10, 79), fields.SList)
    },
)


NUMMDL = Record(
    index=9,
    name="NUMMDL",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=True,
    field_data={
        "num_model": fields.Column((10, 14), fields.Integer),
    }
)


MDLTYP = Record(
    index=10,
    name="MDLTYP",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "description": fields.Column((10, 80), fields.SList)
    },
)


AUTHOR = Record(
    index=11,
    name="AUTHOR",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "author_list": fields.Column((10, 79), fields.List, strip=False)
    },
)


REVDAT = Record(
    index=12,
    name="REVDAT",
    is_mandatory=True,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "mod_num": fields.Column((7, 10), fields.Integer),
        "continuation": fields.Column((10, 12), fields.Continuation),
        "mod_date": fields.Column((13, 22), fields.Date, cast=False),
        "mod_pdb_id": fields.Column((23, 27), fields.IDcode, cast=False),
        "is_init": fields.Column((31, 32), fields.Integer, cast=False),
        "record": fields.Column([(i, i + 6) for i in range(39, 61, 7)], fields.LString, cast=False)
    },
    key_unique="mod_num",
    keys_repeat=("mod_date", "mod_pdb_id", "is_init"),
)


SPRSDE = Record(
    index=13,
    name="SPRSDE",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "sprsde_date": fields.Column((11, 20), fields.Date, only_first=True),
        "pdb_id": fields.Column((21, 25), fields.IDcode, only_first=True),
        "sprsde_pdb_id": fields.Column([(i, i + 4) for i in range(31, 72, 5)], fields.IDcode, only_non_empty=True)
    },
)


JRNL = Record(
    index=14,
    name="JRNL",
    is_mandatory=False,
    is_one_time=True,
    is_one_line=False,
    field_data={
        "tag": fields.Column((12, 16), fields.LString),
        "continuation": fields.Column((16, 18), fields.Continuation),
        "data": fields.Column((19, 80), fields.LString),
    },
)

JRNL_AUTH = Record(
    name="AUTH",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "author": fields.Column((19, 79), fields.List, strip=False)
    }
)

JRNL_TITL = Record(
    name="TITL",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "title": fields.Column((19, 80), fields.String, strip=False)
    }
)

JRNL_EDIT = Record(
    name="EDIT",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "editor": fields.Column((19, 79), fields.List, strip=False)
    }
)

JRNL_REF = Record(
    name="REF",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "pub_name": fields.Column((19, 47), fields.String, strip=False),
        "vol": fields.Column((51, 55), fields.LString, only_first=True),
        "page": fields.Column((56, 61), fields.LString, only_first=True),
        "year": fields.Column((62, 66), fields.Integer, only_first=True),
    },
    keys_repeat=("vol", "page", "year"),
)

JRNL_PUBL = Record(
    name="PUBL",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "pub": fields.Column((19, 70), fields.String, strip=False)
    }
)

JRNL_REFN = Record(
    name="REFN",
    is_one_time=True,
    is_one_line=True,
    field_data={
        "issn_essn": fields.Column((35, 39), fields.LString),
        "issn": fields.Column((40, 65), fields.LString)
    }
)

JRNL_PMID = Record(
    name="PMID",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "pm_id": fields.Column((19, 79), fields.String)
    }
)

JRNL_DOI = Record(
    name="DOI",
    is_one_time=True,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((16, 18), fields.Continuation),
        "doi": fields.Column((19, 79), fields.String)
    }
)


REMARK = Record(
    index=15,
    name="REMARK",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,  # it is actually multiple-lines (False), but its extraction pattern fits one-line
    field_data={
        "remark_num": fields.Column((7, 10), fields.Integer),
        "content": fields.Column((11, 80), fields.LString, cast=False, strip=False),
    }
)


DBREF = Record(
    index=16,
    name="DBREF",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "pdb_id": fields.Column((7, 11), fields.IDcode),
        "chain_id": fields.Column((12, 13), fields.Character),
        "res_num_begin": fields.Column((14, 18), fields.Integer),
        "res_icode_begin": fields.Column((18, 19), fields.AChar),
        "res_num_end": fields.Column((20, 24), fields.Integer),
        "res_icode_end": fields.Column((24, 25), fields.AChar),
        "db": fields.Column((26, 32), fields.LString),
        "db_chain_accession": fields.Column((33, 41), fields.LString),
        "db_chain_id": fields.Column((42, 54), fields.LString),
        "db_res_num_begin": fields.Column((55, 60), fields.Integer),
        "db_res_icode_begin": fields.Column((60, 61), fields.AChar),
        "db_res_num_end": fields.Column((62, 67), fields.Integer),
        "db_res_icode_end": fields.Column((67, 68), fields.AChar),
    },
)


DBREF1 = Record(
    index=17,
    name="DBREF1",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "pdb_id": fields.Column((7, 11), fields.IDcode),
        "chain_id": fields.Column((12, 13), fields.Character),
        "res_num_begin": fields.Column((14, 18), fields.Integer),
        "res_icode_begin": fields.Column((18, 19), fields.AChar),
        "res_num_end": fields.Column((20, 24), fields.Integer),
        "res_icode_end": fields.Column((24, 25), fields.AChar),
        "db": fields.Column((26, 32), fields.LString),
        "db_chain_id": fields.Column((47, 67), fields.LString),
    },
)


DBREF2 = Record(
    index=18,
    name="DBREF2",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "pdb_id": fields.Column((7, 11), fields.IDcode),
        "chain_id": fields.Column((12, 13), fields.Character),
        "db_chain_accession": fields.Column((18, 40), fields.LString),
        "db_res_num_begin": fields.Column((45, 55), fields.Integer),
        "db_res_num_end": fields.Column((57, 67), fields.Integer),
    }
)


SEQADV = Record(
    index=19,
    name="SEQADV",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "pdb_id": fields.Column((7, 11), fields.IDcode),
        "res_name": fields.Column((12, 15), fields.ResidueName),
        "chain_id": fields.Column((16, 17), fields.Character),
        "res_num": fields.Column((18, 22), fields.Integer),
        "res_icode": fields.Column((22, 23), fields.AChar),
        "db": fields.Column((24, 28), fields.LString),
        "db_chain_accession": fields.Column((29, 38), fields.LString),
        "db_res_name": fields.Column((39, 42), fields.ResidueName),
        "db_res_num": fields.Column((43, 48), fields.LString),
        "description": fields.Column((49, 70), fields.LString)
    },
)


SEQRES = Record(
    index=20,
    name="SEQRES",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "serial": fields.Column((7, 10), fields.Integer),
        "chain_id": fields.Column((11, 12), fields.Character),
        "num_res": fields.Column((13, 17), fields.Integer, cast=False),
        "res_name": fields.Column([(c, c + 3) for c in range(19, 68, 4)], fields.ResidueName, cast=False)
    },
    key_unique="chain_id",
    key_continuation="serial",
    keys_repeat=("num_res",)
)


MODRES = Record(
    index=21,
    name="MODRES",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "pdb_id": fields.Column((7, 11), fields.IDcode),
        "res_name": fields.Column((12, 15), fields.ResidueName),
        "chain_id": fields.Column((16, 17), fields.Character),
        "res_num": fields.Column((18, 22), fields.Integer),
        "res_icode": fields.Column((22, 23), fields.AChar),
        "std_res_name": fields.Column((24, 27), fields.ResidueName),
        "description": fields.Column((29, 70), fields.LString)
    },
)


HET = Record(
    index=22,
    name="HET",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "het_id": fields.Column((7, 10), fields.LString),
        "chain_id": fields.Column((12, 13), fields.Character),
        "res_num": fields.Column((13, 17), fields.Integer),
        "res_icode": fields.Column((17, 18), fields.AChar),
        "num_hetatm": fields.Column((20, 25), fields.Integer),
        "description": fields.Column((30, 70), fields.LString)
    }
)


HETNAM = Record(
    index=23,
    name="HETNAM",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "het_id": fields.Column((11, 14), fields.LString),
        "name": fields.Column((15, 70), fields.String, strip=False, cast=False)
    },
    key_unique="het_id"
)


HETSYN = Record(
    index=24,
    name="HETSYN",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "continuation": fields.Column((8, 10), fields.Continuation),
        "het_id": fields.Column((11, 14), fields.LString),
        "synonym": fields.Column((15, 70), fields.SList, cast=False)
    },
    key_unique="het_id"
)


FORMUL = Record(
    index=25,
    name="FORMUL",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "comp_num": fields.Column((8, 10), fields.Integer, cast=False),
        "het_id": fields.Column((12, 15), fields.LString),
        "continuation": fields.Column((16, 18), fields.Continuation),
        "is_water": fields.Column((18, 19), fields.Character, cast=False),
        "formula": fields.Column((19, 70), fields.String, strip=False, cast=False)
    },
    key_unique="het_id",
    keys_repeat=("comp_num", "is_water")
)


HELIX = Record(
    index=26,
    name="HELIX",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((7, 10), fields.Integer),
        "helix_id": fields.Column((11, 14), fields.LString),
        "res_name_begin": fields.Column((15, 18), fields.ResidueName),
        "chain_id_begin": fields.Column((19, 20), fields.Character),
        "res_num_begin": fields.Column((21, 25), fields.Integer),
        "res_icode_begin": fields.Column((25, 26), fields.AChar),
        "res_name_end": fields.Column((27, 30), fields.ResidueName),
        "chain_id_end": fields.Column((31, 32), fields.Character),
        "res_num_end": fields.Column((33, 37), fields.Integer),
        "res_icode_end": fields.Column((37, 38), fields.AChar),
        "helix_class": fields.Column((38, 40), fields.Integer),
        "description": fields.Column((40, 70), fields.LString),
        "length": fields.Column((71, 76), fields.Integer),
    }
)


SHEET = Record(
    index=27,
    name="SHEET",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "strand": fields.Column((7, 10), fields.Integer),
        "sheet_id": fields.Column((11, 14), fields.LString),
        "num_strand": fields.Column((14, 16), fields.Integer),
        "res_name_begin": fields.Column((17, 20), fields.ResidueName),
        "chain_id_begin": fields.Column((21, 22), fields.Character),
        "res_num_begin": fields.Column((22, 26), fields.Integer),
        "res_icode_begin": fields.Column((26, 27), fields.AChar),
        "res_name_end": fields.Column((28, 31), fields.ResidueName),
        "chain_id_end": fields.Column((32, 33), fields.Character),
        "res_num_end": fields.Column((33, 37), fields.Integer),
        "res_icode_end": fields.Column((37, 38), fields.AChar),
        "sense": fields.Column((38, 40), fields.Integer),
        "atom_name_cur": fields.Column((41, 45), fields.Atom),
        "res_name_cur": fields.Column((45, 48), fields.ResidueName),
        "chain_id_cur": fields.Column((49, 50), fields.Character),
        "res_num_cur": fields.Column((50, 54), fields.Integer),
        "res_icode_cur": fields.Column((54, 55), fields.AChar),
        "atom_name_pre": fields.Column((56, 60), fields.Atom),
        "res_name_pre": fields.Column((60, 63), fields.ResidueName),
        "chain_id_pre": fields.Column((64, 65), fields.Character),
        "res_num_pre": fields.Column((65, 69), fields.Integer),
        "res_icode_pre": fields.Column((69, 70), fields.AChar)
    }
)


SSBOND = Record(
    index=28,
    name="SSBOND",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((7, 10), fields.Integer),
        "res_name_1": fields.Column((11, 14), fields.LString),
        "chain_id_1": fields.Column((15, 16), fields.Character),
        "res_num_1": fields.Column((17, 21), fields.Integer),
        "res_icode_1": fields.Column((21, 22), fields.Atom),
        "res_name_2": fields.Column((25, 28), fields.LString),
        "chain_id_2": fields.Column((29, 30), fields.Character),
        "res_num_2": fields.Column((31, 35), fields.Integer),
        "res_icode_2": fields.Column((35, 36), fields.AChar),
        "sym_op_1": fields.Column((59, 62), fields.Integer),
        "trans_vec_1": fields.Column([(i, i+1) for i in range(62, 65)], fields.SymOP, strip=False),
        "sym_op_2": fields.Column((66, 69), fields.Integer),
        "trans_vec_2": fields.Column([(i, i+1) for i in range(69, 72)], fields.SymOP, strip=False),
        "length": fields.Column((73, 78), fields.Real),
    }
)


LINK = Record(
    index=29,
    name="LINK",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "atom_name_1": fields.Column((12, 16), fields.Atom),
        "alt_loc_1": fields.Column((16, 17), fields.Character),
        "res_name_1": fields.Column((17, 20), fields.ResidueName),
        "chain_id_1": fields.Column((21, 22), fields.Character),
        "res_num_1": fields.Column((22, 26), fields.Integer),
        "res_icode_1": fields.Column((26, 27), fields.AChar),
        "atom_name_2": fields.Column((42, 46), fields.Atom),
        "alt_loc_2": fields.Column((46, 47), fields.Character),
        "res_name_2": fields.Column((47, 50), fields.ResidueName),
        "chain_id_2": fields.Column((51, 52), fields.Character),
        "res_num_2": fields.Column((52, 56), fields.Integer),
        "res_icode_2": fields.Column((56, 57), fields.AChar),
        "sym_op_1": fields.Column((59, 62), fields.Integer),
        "trans_vec_1": fields.Column([(i, i+1) for i in range(62, 65)], fields.SymOP, strip=False),
        "sym_op_2": fields.Column((66, 69), fields.Integer),
        "trans_vec_2": fields.Column([(i, i+1) for i in range(69, 72)], fields.SymOP, strip=False),
        "length": fields.Column((73, 78), fields.Real)
    }
)


CISPEP = Record(
    index=30,
    name="CISPEP",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((7, 10), fields.Integer),
        "res_name_1": fields.Column((11, 14), fields.LString),
        "chain_id_1": fields.Column((15, 16), fields.Character),
        "res_num_1": fields.Column((17, 21), fields.Integer),
        "res_icode_1": fields.Column((21, 22), fields.AChar),
        "res_name_2": fields.Column((25, 28), fields.LString),
        "chain_id_2": fields.Column((29, 30), fields.Character),
        "res_num_2": fields.Column((31, 35), fields.Integer),
        "res_icode_2": fields.Column((35, 36), fields.AChar),
        "model_num": fields.Column((43, 46), fields.Integer),
        "angle": fields.Column((53, 59), fields.Real)
    }
)


SITE = Record(
    index=31,
    name="SITE",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=False,
    field_data={
        "serial": fields.Column((7, 10), fields.Integer),
        "site_id": fields.Column((11, 14), fields.LString),
        "num_res": fields.Column((15, 17), fields.Integer, cast=False),
        "res_name": fields.Column(((18, 21), (29, 32), (40, 43), (51, 54)), fields.ResidueName, cast=False),
        "chain_id": fields.Column(((22, 23), (33, 34), (44, 45), (55, 56)), fields.Character, cast=False),
        "res_num": fields.Column(((23, 27), (34, 38), (45, 49), (56, 60)), fields.Integer, cast=False),
        "res_icode": fields.Column(((27, 28), (38, 39), (49, 50), (60, 61)), fields.AChar, cast=False),
    },
    key_unique="site_id",
    key_continuation="serial",
    keys_repeat=("num_res",)
)


CRYST1 = Record(
    index=32,
    name="CRYST1",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data={
        "a": fields.Column((6, 15), fields.Real),
        "b": fields.Column((15, 24), fields.Real),
        "c": fields.Column((24, 33), fields.Real),
        "alpha": fields.Column((33, 40), fields.Real),
        "beta": fields.Column((40, 47), fields.Real),
        "gamma": fields.Column((47, 54), fields.Real),
        "space_group": fields.Column((55, 66), fields.LString),
        "z": fields.Column((66, 70), fields.Integer),
    }
)

_field_data_xform = {
    "m": fields.Column([(i, i+10) for i in range(10, 31, 10)], fields.Real),
    "v": fields.Column((45, 55), fields.Real)
}

ORIGX1 = Record(
    index=33,
    name="ORIGX1",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


ORIGX2 = Record(
    index=34,
    name="ORIGX2",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


ORIGX3 = Record(
    index=35,
    name="ORIGX3",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


SCALE1 = Record(
    index=36,
    name="SCALE1",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


SCALE2 = Record(
    index=37,
    name="SCALE2",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


SCALE3 = Record(
    index=38,
    name="SCALE3",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data=_field_data_xform
)


_field_data_mtrix = _field_data_xform | {
    "serial": fields.Column((7, 10), fields.Integer),
    "is_given": fields.Column((59, 60), fields.Character),
}

MTRIX1 = Record(
    index=39,
    name="MTRIX1",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data=_field_data_mtrix
)


MTRIX2 = Record(
    index=40,
    name="MTRIX2",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data=_field_data_mtrix
)


MTRIX3 = Record(
    index=41,
    name="MTRIX3",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data=_field_data_mtrix
)


MODEL = Record(
    index=42,
    name="MODEL",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={"serial": fields.Column((10, 14), fields.Integer)}
)


ATOM = Record(
    index=43,
    name="ATOM",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((6, 11), fields.Integer),
        "atom_name": fields.Column((12, 16), fields.Atom),
        "alt_loc": fields.Column((16, 17), fields.Character),
        "res_name": fields.Column((17, 20), fields.ResidueName),
        "chain_id": fields.Column((21, 22), fields.Character),
        "res_num": fields.Column((22, 26), fields.Integer),
        "res_icode": fields.Column((26, 27), fields.AChar),
        "x": fields.Column((30, 38), fields.Real),
        "y": fields.Column((38, 46), fields.Real),
        "z": fields.Column((46, 54), fields.Real),
        "occupancy": fields.Column((54, 60), fields.Real),
        "temp_factor": fields.Column((60, 66), fields.Real),
        "element": fields.Column((76, 78), fields.Element),
        "charge": fields.Column(((79, 80), (78, 79)), fields.LString)
    },
)


ANISOU = Record(
    index=44,
    name="ANISOU",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((6, 11), fields.Integer),
        "atom_name": fields.Column((12, 16), fields.Atom),
        "alt_loc": fields.Column((16, 17), fields.Character),
        "res_name": fields.Column((17, 20), fields.ResidueName),
        "chain_id": fields.Column((21, 22), fields.Character),
        "res_num": fields.Column((22, 26), fields.Integer),
        "res_icode": fields.Column((26, 27), fields.AChar),
        "u00": fields.Column((28, 35), fields.Integer),
        "u11": fields.Column((35, 42), fields.Integer),
        "u22": fields.Column((42, 49), fields.Integer),
        "u01": fields.Column((49, 56), fields.Integer),
        "u02": fields.Column((56, 63), fields.Integer),
        "u12": fields.Column((63, 70), fields.Integer),
        "element": fields.Column((76, 78), fields.Element),
        "charge": fields.Column(((79, 80), (78, 79)), fields.LString)
    }
)


TER = Record(
    index=45,
    name="TER",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((6, 11), fields.Integer),
        "res_name": fields.Column((17, 20), fields.ResidueName),
        "chain_id": fields.Column((21, 22), fields.Character),
        "res_num": fields.Column((22, 26), fields.Integer),
        "res_icode": fields.Column((26, 27), fields.AChar),
    }
)


HETATM = Record(
    index=46,
    name="HETATM",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((6, 11), fields.Integer),
        "atom_name": fields.Column((12, 16), fields.Atom),
        "alt_loc": fields.Column((16, 17), fields.Character),
        "res_name": fields.Column((17, 20), fields.ResidueName),
        "chain_id": fields.Column((21, 22), fields.Character),
        "res_num": fields.Column((22, 26), fields.Integer),
        "res_icode": fields.Column((26, 27), fields.AChar),
        "u00": fields.Column((28, 35), fields.Integer),
        "u11": fields.Column((35, 42), fields.Integer),
        "u22": fields.Column((42, 49), fields.Integer),
        "u01": fields.Column((49, 56), fields.Integer),
        "u02": fields.Column((56, 63), fields.Integer),
        "u12": fields.Column((63, 70), fields.Integer),
        "element": fields.Column((76, 78), fields.Element),
        "charge": fields.Column(((79, 80), (78, 79)), fields.LString)
    },
)


ENDMDL = Record(
    index=47,
    name="ENDMDL",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data=dict()
)


CONECT = Record(
    index=48,
    name="CONECT",
    is_mandatory=False,
    is_one_time=False,
    is_one_line=True,
    field_data={
        "serial": fields.Column((7, 11), fields.Integer),
        "serial_bonded": fields.Column([(i, i+4) for i in range(12, 28, 5)], fields.Integer, cast=False),
    },
)


MASTER = Record(
    index=49,
    name="MASTER",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
    field_data={
        "remark": fields.Column((10, 15), fields.Integer),
        "0": fields.Column((15, 20), fields.Integer),
        "het": fields.Column((20, 25), fields.Integer),
        "helix": fields.Column((25, 30), fields.Integer),
        "sheet": fields.Column((30, 35), fields.Integer),
        "turn": fields.Column((35, 40), fields.Integer),
        "site": fields.Column((40, 45), fields.Integer),
        "xform": fields.Column((45, 50), fields.Integer),
        "coord": fields.Column((50, 55), fields.Integer),
        "ter": fields.Column((55, 60), fields.Integer),
        "conect": fields.Column((60, 65), fields.Integer),
        "seqres": fields.Column((65, 70), fields.Integer),
    },
)


END = Record(
    index=50,
    name="END",
    is_mandatory=True,
    is_one_time=True,
    is_one_line=True,
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






