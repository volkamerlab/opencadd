"""
Parser for PDB files.
"""

from pathlib import Path
import functools
import datetime
import re
from typing import Literal, NoReturn, Optional, Any, Union, Tuple
import warnings

import numpy as np
import pandas as pd

from . import struct, _records, _fields
from opencadd import _exceptions, _typing


__author__ = "Armin Ariamajd"


class PDBParsingError(Exception):
    pass


class PDBParser:
    """
    Parser for PDB files.
    """

    def __init__(self, content: str, strictness: Literal[0, 1, 2, 3] = 0):
        """
        Parameters
        ----------
        content : str
            String content of the PDB file.
        strictness : {0, 1, 2, 3}, optional, default: 0
            Level of strictness for raising exceptions and warnings when encountering mistakes
            in the PDB file:

            * 0: raise only fatal errors and don't show any warnings.
            * 1: raise only fatal errors. All other errors are reported as warnings.
            * 2: raise fatal errors and mistakes resulting in ambiguous data.
                 Inconsequential mistakes are reported as warnings.
            * 3: completely validate the PDB file and raise all errors.
        """
        if not isinstance(content, str):
            _exceptions.raise_type(
                param_name="content",
                parent_name="PDBParser",
                expected_type=str,
                param_arg=content
            )

        self.strictness: Literal[0, 1, 2, 3] = strictness
        self._lines: np.ndarray = np.array(content.splitlines())
        """Lines of the PDB file, as a 1D array of strings"""
        self._lines_chars: np.ndarray = self._lines.view(dtype=(str, 1)).reshape(self._lines.size, -1)
        """Character view of the PDB file, i.e. 2-d array where first axis (rows) is lines and second axis 
        (columns) is characters. If the lines had different lengths, then all other lines will be padded with
        empty characters from the right, so that the shape of the array is (num_lines, max_num_chars_per_line).
        In a standard PDB file, all lines are 80 characters, so the shape will be (num_lines, 80).
        """
        self._lines_rec_names = np.char.strip(self._lines_chars[:, :6].view(dtype=(str, 6)).reshape(-1))
        """Record name of each line in the PDB file, stripped of whitespaces.
        For example: array(['HEADER', 'TITLE', 'TITLE', ..., 'CONECT', 'MASTER', 'END'])"""
        self._lines_rec_content: np.ndarray = self._lines_chars[:, 6:].view(dtype=(str, 74)).reshape(-1)
        """Record content of each line in the PDB file, i.e. with the 6 first characters stripped"""
        self._mask_records: np.ndarray = np.expand_dims(self._lines_rec_names, axis=-1) == _records.names

        record_line_is_valid = np.isin(self._lines_rec_names, _records.names)
        if not np.any(record_line_is_valid):
            raise PDBParsingError("Input file is not a valid PDB file as it contains no records.")

        if self.strictness:
            if not np.all(record_line_is_valid):
                self._raise_or_warn(
                    "Input file contains invalid lines. "
                    "Following lines do not start with a valid record name: "
                    f"{self._lines[np.logical_not(record_line_is_valid)]}",
                    raise_level=2
                )

            # Check if all lines are 80 characters:
            line_widths = np.char.str_len(self._lines)
            lines_not80 = line_widths != 80
            if np.any(lines_not80):
                faulty_lines = [
                    f"{line_num} ({line_width})" for line_num, line_width in zip(
                        np.argwhere(lines_not80).reshape(-1) + 1,
                        line_widths[lines_not80]
                    )
                ]
                self._raise_or_warn(
                    "Following lines are not 80 characters long (lengths given in brackets): "
                    f"{', '.join(faulty_lines)}",
                    raise_level=3
                )

        idx_lines_records = np.argwhere(self._mask_records)
        self._idx__record_lines = dict()
        self._count_records: np.ndarray = np.empty(
            shape=_records.count,
            dtype=_typing.smallest_integer_dtype_for_range(0, self._lines.size)
        )
        self._has_record: np.ndarray = np.empty(shape=_records.count, dtype=np.bool_)
        for record_idx in range(_records.count):
            indices_record_lines = idx_lines_records[idx_lines_records[:, 1] == record_idx, 0]
            count_record = indices_record_lines.size
            self._idx__record_lines[record_idx] = indices_record_lines
            self._count_records[record_idx] = count_record
            self._has_record[record_idx] = count_record != 0

        self._remarks: dict = None
        self._remark_line_indices: dict = None
        self._remark_args: dict = None
        return

    def parse(self, records) -> struct.PDBStructure:
        return struct.PDBStructure(**{record: getattr(self, record)() for record in records})

    def header(self) -> Optional[struct.RecordHeader]:
        """
        Parse the HEADER record of the PDB file.
        """
        def parse_classification(classification: str) -> Tuple[Tuple[str, ...], ...]:
            """Parse the classification field of the HEADER record."""
            # The classification string is left-justified, and can describe dual functions of molecules
            # (when applicable) separated by a comma “,”. Entries with multiple molecules in a complex
            # will list the classifications of each macromolecule separated by slash “/”.
            # Thus, first, split by '/' and then by ',' to get a tuple of tuples:
            class_per_entity = tuple(classification.split("/"))
            class_per_entity_and_function = tuple(
                tuple(entity_function.strip() for entity_function in entity_class.split(","))
                for entity_class in class_per_entity
            )
            return class_per_entity_and_function
        # Extract record fields as a dict
        if not self.has_record(_records.HEADER):
            return
        data = _records.HEADER.extract(record_char_table=self.record_lines(_records.HEADER))
        # Parse the classification string
        data["classification"] = parse_classification(data["classification"])
        return struct.RecordHeader(**data)

    def obslte(self) -> Optional[struct.RecordObslte]:
        """
        Parse the OBSLTE records of the PDB file.
        """
        data = self._extract_record(_records.OBSLTE)
        if data is None:
            return
        return struct.RecordObslte(**data)

    def title(self) -> Optional[str]:
        """
        Parse the TITLE records of the PDB file.
        """
        return self._extract_record(_records.TITLE)

    def split(self) -> Optional[np.ndarray]:
        """
        Parse the SPLIT records of the PDB file.
        """
        return self._extract_record(_records.SPLIT)

    def caveat(self) -> Optional[struct.RecordCaveat]:
        """
        Parse the CAVEAT records of the PDB file.
        """
        data = self._extract_record(_records.CAVEAT)
        if data is None:
            return
        return data["description"]

    def compnd(self) -> Optional[pd.DataFrame]:
        """
        Parse the COMPND records of the PDB file.
        """
        df = self._parse_records_compnd_source(_records.COMPND)
        if df is None:
            return
        engineered_mutation_values = df[["engineered", "mutation"]].dropna()
        value_is_standard = engineered_mutation_values.isin(["YES", "NO"])
        if not value_is_standard.all(axis=None):
            self._raise_or_warn(
                "Non standard values found for 'ENGINEERED'/'MUTATION' tokens of the COMPND record. "
                "Expected values are 'YES' and 'NO', "
                f"but found: {engineered_mutation_values[~value_is_standard]}",
                raise_level=2,
            )
            return df
        df[["engineered", "mutation"]] = df[["engineered", "mutation"]] == "YES"
        return df

    def source(self) -> Optional[pd.DataFrame]:
        """
        Parse the SOURCE records of the PDB file.
        """
        return self._parse_records_compnd_source(_records.SOURCE)

    def keywds(self) -> Optional[np.ndarray]:
        """
        Parse the KEYWDS records of the PDB file.
        """
        return self._extract_record(_records.KEYWDS)

    def expdta(self) -> Optional[np.ndarray]:
        """
        Parse the EXPDTA records of the PDB file.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str] | None
            If the PDB file contains no KEYWDS record, `None` is returned, otherwise an array of strings,
            containing one (for most cases) or several of following allowed values:
            'X-RAY DIFFRACTION', 'FIBER DIFFRACTION', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY',
            'ELECTRON MICROSCOPY', 'SOLID-STATE NMR', 'SOLUTION NMR', 'SOLUTION SCATTERING'.
        """
        techniques = self._extract_record(_records.EXPDTA)
        if techniques is None:
            return
        if self.strictness:
            allowed_values = [
                'X-RAY DIFFRACTION',
                'FIBER DIFFRACTION',
                'NEUTRON DIFFRACTION',
                'ELECTRON CRYSTALLOGRAPHY',
                'ELECTRON MICROSCOPY',
                'SOLID-STATE NMR',
                'SOLUTION NMR',
                'SOLUTION SCATTERING'
            ]
            value_is_unknown = np.isin(techniques, allowed_values, invert=True)
            if np.any(value_is_unknown):
                self._raise_or_warn(
                    "Encountered non-standard values in the 'technique' field of EXPDTA records. "
                    f"Allowed values are {allowed_values}, but found {techniques[value_is_unknown]}.",
                    raise_level=2,
                )
        return techniques

    def nummdl(self) -> Optional[int]:
        """
        Parse the NUMMDL record of the PDB file.

        Returns
        -------
        int | None
            If the PDB file contains no NUMMDL record, `None` is returned,
            otherwise the total number of models in the entry.
        """
        return self._extract_record(_records.NUMMDL)

    def mdltyp(self) -> Optional[np.ndarray]:
        """
        Parse the MDLTYP records of the PDB file.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str] | None
            If the PDB file contains no MDLTYP record, `None` is returned, otherwise
            an array of strings corresponding to a list of annotations.
        """
        return self._extract_record(_records.MDLTYP)

    def author(self) -> Optional[np.ndarray]:
        """
        Parse the AUTHOR records of the PDB file.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str] | None
            If the PDB file contains no MDLTYP record, `None` is returned, otherwise
            an array of strings corresponding to the list of authors.
        """
        return self._extract_record(_records.AUTHOR)

    def revdat(self) -> Optional[pd.DataFrame]:
        """
        Parse the REVDAT records of the PDB file.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no REVDAT record, `None` is returned, otherwise a dataframe with columns:
        """
        df = self._extract_record(_records.REVDAT)
        if df is None:
            return
        if not df["is_init"].isin([1, 0]).all():
            self._raise_or_warn(
                f"Non standard values found for 'modType' field of the REVDAT record, values are: {df.is_init}",
                raise_level=2,
            )
            return df
        df["is_init"] = df["is_init"] == 0
        return df

    def sprsde(self) -> Optional[struct.RecordSPRSDE]:
        """
        Parse the SPRSDE records of the PDB file.
        """
        data = self._extract_record(_records.SPRSDE)
        if data is None:
            return
        return struct.RecordSPRSDE(**data)

    def jrnl(self):
        """
        JRNL record of the PDB file, containing the primary literature citation that describes the experiment
        corresponding to the entry.

        Returns
        -------

        Notes
        -----
        * There is at most one JRNL reference per entry. If there is no primary reference,
        then there is no JRNL reference.
        """
        if not self.has_record(_records.JRNL):
            return
        record_lines = self.record_lines(_records.JRNL)
        return self.records_publication(record_lines)

    @staticmethod
    def records_publication(record_lines: np.ndarray):
        tags = _records.JRNL.fields_dict["tag"].extract(char_table=record_lines)
        ref = dict()
        for sub_record, sub_record_name in (
                (_records.JRNL_AUTH, "author"),
                (_records.JRNL_TITL, "title"),
                (_records.JRNL_EDIT, "editor"),
                (_records.JRNL_PUBL, "pub"),
                (_records.JRNL_PMID, "pm_id"),
                (_records.JRNL_DOI, "doi"),
        ):
            is_sub_record_line = tags == sub_record.name
            if np.any(is_sub_record_line):
                sub_record_val = sub_record.extract(record_char_table=record_lines[is_sub_record_line])
                ref[sub_record_name] = sub_record_val

        is_ref_line = tags == _records.JRNL_REF.name
        if np.any(is_ref_line):
            ref_dict = _records.JRNL_REF.extract(record_char_table=record_lines[is_ref_line])
            if ref_dict["pub_name"] != "TO BE PUBLISHED":
                ref = ref | ref_dict

        is_refn_line = tags == _records.JRNL_REFN.name
        if np.any(is_refn_line):
            refn_dict = _records.JRNL_REFN.extract(record_char_table=record_lines[is_refn_line])
            issn_essn = refn_dict.pop("issn_essn")
            if issn_essn == "ISSN":
                ref = ref | refn_dict
            elif issn_essn == "ESSN":
                refn_dict["essn"] = refn_dict.pop("issn")
                ref = ref | refn_dict
        return struct.RecordJRNL(**ref)

    def remark(self) -> Optional[struct.RecordREMARK]:
        """

        Returns
        -------

        """
        if not self.has_record(_records.REMARK):
            return
        if self._remark_args is not None:
            return struct.RecordREMARK(**self._remark_args)
        remark_lines = self.record_lines(_records.REMARK)
        is_not_empty = np.any(remark_lines[:, 11:] != " ", axis=1)
        non_empty_lines = remark_lines[is_not_empty]
        idx_nonempty_lines = self.indices_record_lines(_records.REMARK)[is_not_empty]
        data = _records.REMARK.extract(record_char_table=non_empty_lines)
        unique_nums = np.unique(data["remark_num"])
        remarks = dict()
        remark_line_indices = dict()
        args = {"full_text": remarks}
        for unique_num in unique_nums:
            is_num = data["remark_num"] == unique_num
            remarks[unique_num] = data["content"][is_num]
            remark_line_indices[unique_num] = idx_nonempty_lines[is_num]
            if unique_num == 1:
                pass
            elif unique_num == 2:
                args["resolution"] = self.record_remark2(remarks[2])
            elif unique_num == 4:
                args["version"] = self.remark4(idx_lines=remark_line_indices[4])
        self._remarks = remarks
        self._remark_line_indices = remark_line_indices
        self._remark_args = args
        return struct.RecordREMARK(**args)

    # def _parse_remark_records(self):
    #     remark_lines = self.record_lines(_records.REMARK)
    #     is_not_empty = np.any(remark_lines[:, 11:] != " ", axis=1)
    #     idx_nonempty_lines = self.indices_record_lines(_records.REMARK)[is_not_empty]
    #     remark_nums = _records.REMARK.fields_dict["remark_num"].extract(char_table=remark_lines[is_not_empty])
    #     self._remarks_records = {"idx_line": idx_nonempty_lines, "remark_num": remark_nums}
    #     return
    #
    # def _extract_remark_record(self, remark_num: int):
    #     if not self.has_record(_records.REMARK):
    #         return
    #     if self._remarks_records is None:
    #         self._parse_remark_records()
    #     hits = self._remarks_records["remark_num"] == remark_num
    #     if not np.any(hits):
    #         return
    #     idx_lines = self._remarks_records["idx_line"][hits]
    #     return idx_lines

    def record_remark2(self, lines: np.ndarray) -> Optional[float]:
        """
        REMARK 2 record of the PDB file, stating the highest resolution that was used in building the model.

        Returns
        -------
        float | Literal[0] | None
            If the PDB file contains no REMARK 2 record, `None` is returned, otherwise
            the resolution of the model, in Angstroms (Å), unless the record contains the tag
            "RESOLUTION. NOT APPLICABLE." (this is the case for experiments other than diffraction studies,
            e.g. NMR, where resolution is not relevant), in which case 0 is returned.
        """
        # idx_lines = self._extract_remark_record(2)
        # if idx_lines is None:
        #     return
        if lines.size > 1:
            self._raise_or_warn(f"Multiple REMARK 2 records found: {lines}", raise_level=2)
        resolution = lines[0][12:19]
        if resolution == "NOT APP":
            return 0
        try:
            return float(resolution)
        except:
            self._raise_or_warn(f"REMARK 2 record has non-standard format: {lines}", raise_level=2)
        return None

    def remark4(self, idx_lines: np.ndarray) -> Optional[dict[str, Union[str, str, datetime.date]]]:
        """
        REMARK 4 record of the PDB file, indicating the version of the PDB file format
        used to generate the file.

        Returns
        -------
        dict[str, str | str | datetime.date] | None
            If the PDB file contains no REMARK 4 record, `None` is returned, otherwise a dictionary with keys:
            pdb_id : str
                PDB ID of the entry.
            version_num : str
                Version number of the PDB file format used, e.g. "3.20".
            version_date : datetime.date
                Release date of the version.
        """
        lines = self._lines[idx_lines]
        if lines.size > 1:
            self._raise_or_warn(
                "REMARK 4 occurs multiple times; only the first occurrence will be parsed.",
                raise_level=2
            )
        data = {
            field_name: lines[0][start:stop] for field_name, (start, stop) in zip(
                ("pdb_id", "ver_num", "ver_date"),
                ((11, 15), (40, 44), (46, 55))
            )
        }
        data["ver_date"] = _fields.Date.from_pdb(data["ver_date"])
        return data

    def dbref(self) -> Optional[pd.DataFrame]:
        """
        Parse the DBREF and DBREF1/DBREF2 records of the PDB file.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains neither DBREF nor DBREF1/DBREF2 records, `None` is returned,
            otherwise a dataframe with columns:
        """
        if not (self.has_record(_records.DBREF) or self.has_record(_records.DBREF1)):
            return
        # Information may be in both DBREF and DBREF1/DBREF2 records.
        #  These are extracted separately, and merged:
        df_main = pd.DataFrame(columns=_records.DBREF.fields_dict.keys())
        if self.has_record(_records.DBREF):
            # Extract DBREF records; no post-extraction processing needed.
            df_main = pd.concat(
                (df_main, self._extract_record(_records.DBREF))
            )
        if self.has_record(_records.DBREF1):
            # Extract DBREF1/DBREF2 records and merge rows to obtain the corresponding DBREF-like dataframe.
            df = pd.merge(
                self._extract_record(_records.DBREF1),
                self._extract_record(_records.DBREF2),
                how="outer",
                on=["pdb_id", "chain_id"]
            )
            df_main = pd.concat((df_main, df))
        return df_main.set_index("chain_id", drop=False)

    def seqadv(self) -> Optional[pd.DataFrame]:
        """
        Parse the SEQADV records of the PDB file.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no SEQADV record, `None` is returned, otherwise a dataframe with columns:
        """
        return self._extract_record(_records.SEQADV)

    def seqres(self) -> Optional[pd.DataFrame]:
        if self.strictness:
            if not self.has_record(_records.SEQRES):
                if self.has_record(_records.ATOM):
                    self._raise_or_warn(
                        "SEQRES records not found, while ATOM records exist. "
                        "SEQRES is a mandatory record when ATOM records exist.",
                        raise_level=3
                    )
            #self._check_routine(_records.SEQRES)
        return self._extract_record(_records.SEQRES)

    def modres(self) -> Optional[pd.DataFrame]:
        return self._extract_record(_records.MODRES)

    def het(self) -> Optional[pd.DataFrame]:
        data = self._extract_record(_records.HET)
        return data.set_index("het_id", drop=False) if data is not None else None

    def hetnam(self):
        hetnam = self._extract_record(_records.HETNAM)
        hetsyn = self._extract_record(_records.HETSYN)
        formul = self.record_formul()
        het_data = pd.DataFrame(
            columns=[
                "comp_num", "het_id", "name", "synonym", "is_water", "formula", "count_in_chain", "count_rest"
            ]
        )
        het_data.index.name = "het_id"
        non_none = [rec for rec in (hetnam, hetsyn, formul) if rec is not None]
        if len(non_none) == 0:
            return het_data
        if len(non_none) == 1:
            df_tot = pd.concat([het_data, non_none[0]]).set_index("het_id", drop=False)
            df_tot.loc[df_tot.is_water, "name"] = "WATER"
            return df_tot.sort_values("comp_num")
        df_tot = non_none[0].drop("het_id", axis=1)
        for df in non_none[1:]:
            df_tot = pd.merge(df_tot, df.drop("het_id", axis=1), how="outer", left_index=True, right_index=True)
        df_tot["het_id"] = df_tot.index
        df_tot.loc[df_tot.is_water, "name"] = "WATER"
        return pd.concat([het_data, df_tot]).sort_values("comp_num")

    def record_formul(self) -> Optional[pd.DataFrame]:

        def parse_formula_field(text: str):
            """

            Parameters
            ----------
            text

            Returns
            -------
            str, int, int
                Chemical formula, number of occurrences within a chain, remaining number of occurrences.
            """
            # All occurrences of the HET group within a chain are grouped together with a
            # multiplier. The remaining occurrences are also grouped with a multiplier, e.g.:
            # 4(2(C21 H28 N7 O17 P3)), 2(C19 H19 N7 O6) , or C10 H22 O6.
            formula_split = re.split(r"[()]+", text)
            num_split = len(formula_split)
            if num_split == 1:
                return formula_split[0], 1, 0
            elif num_split == 3:
                return formula_split[-2], formula_split[0], 0
            elif num_split == 4:
                return formula_split[-2], formula_split[1], formula_split[0]
            else:
                self._raise_or_warn(
                    f"Could not parse the 'text' field of the FORMUL record; got: {text}.",
                    raise_level=2
                )
                return text, 0, 0

        data: Optional[pd.DataFrame] = self._extract_record(_records.FORMUL)
        if data is None:
            return
        formulas, counts_in_chain, counts_rest = [], [], []
        for raw_formula in data["formula"]:
            formula, count_in_chain, count_rest = parse_formula_field(raw_formula)
            formulas.append(formula)
            counts_in_chain.append(count_in_chain)
            counts_rest.append(count_rest)
        data["formula"] = formulas
        data["count_in_chain"] = counts_in_chain
        data["count_rest"] = counts_rest
        data["is_water"] = data["is_water"] == "*"
        return data

    def helix(self) -> Optional[pd.DataFrame]:
        data = self._extract_record(_records.HELIX)
        return data.set_index("helix_id") if data is not None else None

    def sheet(self):
        data = self._extract_record(_records.SHEET)
        return data.set_index("sheet_id") if data is not None else None

    def ssbond(self):
        data = self._extract_record(_records.SSBOND)
        return data.set_index("serial") if data is not None else None

    def link(self):
        return self._extract_record(_records.LINK)

    def cispep(self):
        data = self._extract_record(_records.CISPEP)
        return data.set_index("serial") if data is not None else None

    def site(self):
        data = self._extract_record(_records.SITE)
        return data.explode(["res_name", "chain_id", "res_num"]) if data is not None else None

    def cryst1(self) -> struct.RecordCRYST1:
        data = self._extract_record(_records.CRYST1)
        lengths = np.array([data["a"], data["b"], data["c"]])
        angles = np.array([data["alpha"], data["beta"], data["gamma"]])
        return struct.RecordCRYST1(lengths=lengths, angles=angles, z=data["z"], space_group=data["space_group"])

    def origx(self) -> Optional[struct.RecordXForm]:
        o1 = self._extract_record(_records.ORIGX1)
        o2 = self._extract_record(_records.ORIGX2)
        o3 = self._extract_record(_records.ORIGX3)
        if (o1 is None) or (o2 is None) or (o3 is None):
            return
        transformation_matrix = np.vstack([o["m"] for o in (o1, o2, o3)])
        translation_vector = np.array([o["v"] for o in (o1, o2, o3)])
        return struct.RecordXForm(matrix=transformation_matrix, vector=translation_vector)

    def scale(self) -> Optional[struct.RecordXForm]:
        s1 = self._extract_record(_records.SCALE1)
        s2 = self._extract_record(_records.SCALE2)
        s3 = self._extract_record(_records.SCALE3)
        if (s1 is None) or (s2 is None) or (s3 is None):
            return
        transformation_matrix = np.vstack([o["m"] for o in (s1, s2, s3)])
        translation_vector = np.array([s1["v"], s2["v"], s3["v"]])
        return struct.RecordXForm(matrix=transformation_matrix, vector=translation_vector)

    def mtrix(self):
        ms = [self._extract_record(rec) for rec in (_records.MTRIX1, _records.MTRIX2, _records.MTRIX3)]
        if any(m is None for m in ms):
            return
        shared_serials = ms[0]["serial"]
        for n, m in enumerate(ms[1:]):
            if m["serial"].size != shared_serials.size:
                self._raise_or_warn("Number of MTRIX records doesn't match.", raise_level=2)
            if not np.all(np.isin(m["serial"], shared_serials)):
                self._raise_or_warn("MTIRX records' serial numbers do not match.", raise_level=2)
            shared_serials = np.intersect1d(shared_serials, m["serial"])

        idx = [
            np.argwhere(shared_serials == np.expand_dims(m["serial"], axis=-1))[:, 1]
            for m in ms
        ]
        matrices = np.stack([ms[n]["m"][idx[n]] for n in range(3)], axis=1)
        vectors = np.stack([ms[n]["v"][idx[n]] for n in range(3)], axis=1)
        is_given = np.any([ms[n]["is_given"][idx[n]] == "1" for n in range(3)], axis=0)
        return pd.DataFrame(
            {"serial": shared_serials, "is_given":is_given, "matrix": list(matrices), "vector": list(vectors)}
        ).set_index("serial", drop=False)

    # def atom(self):
    #     coordinates = self.coordinates()
    #     if coordinates is None:
    #         return self.anisou()
    #     anisou = self.anisou()
    #     if anisou is None:
    #         return coordinates
    #     return pd.merge(
    #         coordinates, anisou,
    #         how="outer",
    #         on=("atom_name", "alt_loc", "res_name", "chain_id", "res_num", "res_icode", "element", "charge", "model_num")
    #     )

    def atom(self) -> Optional[pd.DataFrame]:
        """
        ATOM and HETATM records of the PDB file, containing atomistic data (coordinates, occupancy, and
        temperature factor) of all atoms present in the entry.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains neither ATOM nor HETATM records, `None` is returned,
            otherwise a dataframe with columns:

            serial (index) : int
                Serial number of the record.
            model_num : int
                Number of the model containing the atom. If the entry has only one model,
                this will be 1 for all atoms.
            chain_id : str
                Chain identifier of the chain containing the atom.
            res_name : str
                Name of the residue containing the atom.
            res_num : int
                Number of the last residue.
            res_icode : str
                Insertion code of the last residue.

        Notes
        -----
        * In each protein/nucleic-acid chain, atoms are listed from the amino-/5'-terminus
          to the carboxy-/3'-terminus, respectively. No ordering is specified for polysaccharides.

        """
        if not (self.has_record(_records.ATOM) or self.has_record(_records.HETATM)):
            return
        mask_lines = np.logical_or(
            self.mask_record(_records.ATOM), self.mask_record(_records.HETATM)
        )
        idx_lines = np.argwhere(mask_lines).reshape(-1)
        lines = self._lines_chars[mask_lines]
        data = _records.ATOM.extract(record_char_table=lines)
        data["res_std"] = self._lines_rec_names[mask_lines] == "ATOM"
        data["model_num"] = self._model_num_of_lines(idx_lines)
        data["charge"] = self._parse_field_charge(data["charge"])

        if self.has_record(_records.TER):
            idx_ter_lines = self.indices_record_lines(_records.TER)
            if self.has_record(_records.MODEL):
                is_polymer = np.ones(shape=idx_lines.size, dtype=np.bool_)
                indices_model_start = self.indices_record_lines(_records.MODEL)
                indices_model_end = self.indices_record_lines(_records.ENDMDL)
                for idx_start, idx_end in zip(indices_model_start, indices_model_end):
                    mask_ter = np.logical_and(idx_ter_lines > idx_start, idx_ter_lines < idx_end)
                    if mask_ter.size > 0:
                        mask_not_polymer = np.logical_and(
                            idx_lines > idx_ter_lines[mask_ter][-1],
                            idx_lines < idx_end
                        )
                        is_polymer[mask_not_polymer] = False
                data["res_poly"] = is_polymer
            else:
                data["res_poly"] = idx_lines < idx_ter_lines[-1]
        else:
            if self.has_record(_records.ATOM):
                self._raise_or_warn("ATOM without TER", raise_level=2)
                data["res_poly"] = data["res_std"]
            else:
                data["res_poly"] = False

        df = pd.DataFrame(data).set_index("serial", drop=False)[
            [
                "model_num",
                "chain_id",
                "res_name",
                "res_num",
                "res_icode",
                "res_poly",
                "res_std",
                "serial",
                "atom_name",
                "alt_loc",
                "occupancy",
                "element",
                "charge",
                "x",
                "y",
                "z",
                "temp_factor",
            ]
        ]
        return df

    def anisou(self) -> Optional[pd.DataFrame]:
        """
        ANISOU records of the PDB file, containing the anisotropic temperature factors for atoms in the entry.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no ANISOU record, `None` is returned, otherwise a dataframe with columns:

            serial (index) : int
                Serial number of the record. It is usually, but not necessarily, one greater than
                the serial number of the corresponding ATOM/HETATM record.
            model_num : int
                Number of the model containing the atom. If the entry has only one model,
                this will be 1 for all records.
            chain_id : str
                Chain identifier of the chain containing the atom.
            res_name : str
                Name of the residue containing the atom.
            res_num : int
                Number of the last residue.
            res_icode : str
                Insertion code of the last residue.

        Notes
        -----
        * The anisotropic temperature factors are stored in the same coordinate frame as the atomic
          coordinate records.
        * The standard deviations for the anisotropic temperature factors are related to the corresponding
          ATOM/ HETATM ANISOU temperature factors.
        * ANISOU is an optional record; it is only present when provided by the authors, and even then,
          it may only cover a subset of atoms in the entry.
        """
        data = self._extract_record(_records.ANISOU)
        if data is None:
            return
        data["charge"] = self._parse_field_charge(data["charge"])
        data["model_num"] = self._model_num_of_lines(self.indices_record_lines(_records.ANISOU))
        df = pd.DataFrame(data).set_index("serial", drop=False)[
            [
                "model_num",
                "chain_id",
                "res_name",
                "res_num",
                "res_icode",
                "atom_name",
                "alt_loc",
                "element",
                "charge",
                "serial",
                "u00",
                "u11",
                "u22",
                "u01",
                "u02",
                "u12",
            ]
        ]
        return df

    def ter(self) -> Optional[pd.DataFrame]:
        """
        TER records of the PDB file, indicating the last residue of each chain in the entry.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no TER record, `None` is returned, otherwise a dataframe with columns:

            serial (index) : int
                Serial number of the record. It is one greater than the serial number of
                the ATOM/HETATM record preceding it.
            model_num : int
                Number of the model containing the chain. If the entry has only one model,
                this will be 1 for all records.
            chain_id : str
                Chain identifier of the corresponding chain.
            res_name : str
                Name of the last residue.
            res_num : int
                Number of the last residue.
            res_icode : str
                Insertion code of the last residue.

        Notes
        -----
        * The last residue is the carboxy-terminal residue for proteins,
          and the 3'-terminal residue for nucleic acids.
        * Terminal oxygen atoms are presented as OXT for proteins, and as O5’ or OP3 for nucleic acids.
          These atoms are present only if the last residue in the polymer is truly the last residue in the
          SEQRES.
        * The TER record has the same residue name, chain identifier, sequence number and insertion
          code as the terminal residue (i.e. immediately preceding ATOM or non-water HETATM record).
        """
        df = self._extract_record(_records.TER)
        if df is None:
            return
        df["model_num"] = self._model_num_of_lines(self.indices_record_lines(_records.TER))
        return df.set_index("serial", drop=False)[["model_num", "chain_id", "res_name", "res_num", "res_icode", "serial"]]

    def conect(self):
        """
        CONECT records of the PDB file, specifying connectivity between atoms for which coordinates are
        supplied.

        Connectivity is noted for:

        * Intra-residue bonds in non-standard (i.e. HET) residues, excluding water.
        * Inter-residue bonds between a HET group, and either a standard residue, or another HET group.
        * Disulfide bridges specified in the SSBOND records.
        * Hydrogen bonds between a present hydrogen atom and a hydrogen-bond acceptor atom.

        Returns
        -------

        Notes
        -----
        * No differentiation is made between atoms with delocalized charges (excess negative or positive charge).
        * Atoms specified in the CONECT records have the same serial number as given in the coordinate section.
        * CONECT records are not provided for missing atoms (e.g. symmetry-generated, or due to disorder)
          in the entry.
        """
        # if not self.has_record(_records.CONECT):
        #     return
        # record_lines = self.record_lines(_records.CONECT)
        # serial_bonded = _records.CONECT.fields_dict["serial_bonded"].extract(record_lines)
        # serial_ref = _records.CONECT.fields_dict["serial"].extract(record_lines)
        data = self._extract_record(_records.CONECT)
        if data is None:
            return
        unique_serials = np.unique(data["serial"])
        bonds = dict()
        for unique_serial in unique_serials:
            mask = data["serial"] == unique_serial
            bonded_serials = data["serial_bonded"][mask]
            non_empty_fields = bonded_serials[
                bonded_serials != ""
            ]
            bonds[unique_serial] = non_empty_fields.astype(int)

        serial_bonded = [ser_bond[ser_bond != ""] for ser_bond in data["serial_bonded"]]
        return pd.DataFrame({"serial_1": data["serial"], "serial_2": serial_bonded}).explode("serial_2", ignore_index=True)

    def master(self):
        """
        MASTER record of the PDB file, containing checksums (i.e. counts) of lines
        for certain records in the entry.

        It is a control record for bookkeeping and validation of the PDB file.

        Returns
        -------
        dict[str, int] | None
            If the PDB file contains no MASTER record, `None` is returned, otherwise a dictionary with keys:

            remark : int
                Number of REMARK lines.
            0 : Literal[0]
                A control field that is always 0.
            het : int
                Number of HET lines.
            helix : int
                Number of HELIX lines.
            sheet : int
                Number of SHEET lines.
            turn : Literal[0]
                A deprecated field that is always 0.
            site : int
                Number of SITE lines.
            xform : int
                Number of lines for coordinate transformation records, i.e.
                ORIGXn, SCALEn and MTRIXn (n = 1, 2, 3) lines combined.
            coord : int
                Number of lines for atomic coordinate records, i.e. ATOM and HETATM lines combined.
            ter : int
                Number of TER lines.
            conect : int
                Number of CONECT lines.
            seqres : int
                Number of SEQRES lines.

        Notes
        -----
        * When multiple models are present, only the lines in the first model is counted.
        * The correlation between the record's fields and the keys in the returned dictionary is as follows
          (field: key):

          * numRemark: remark
          * 0: 0
          * numHet: het
          * numHelix: helix
          * numSheet: sheet
          * numTurn: turn
          * numSite: site
          * numXform: xform
          * numCoord: coord
          * numTer: ter
          * numConect: conect
          * numSeq: seqres
        """
        return self._extract_record(_records.MASTER)

    def _parse_records_compnd_source(self, record: _records.Record) -> Optional[pd.DataFrame]:
        """
        Parse either the COMPND or the SOURCE records of the PDB file.
        This method is the backend call of both `self.record_compnd` and `self.record_source`,
        since both these records have an identical structure.

        Parameters
        ----------
        record : opencadd.io.pdb.records.Record
            Either `COMPND` or `SOURCE`.

        Returns
        -------
        pandas.DataFrame | None
        """
        def parse_token_value_pairs(start_: int, stop_: int) -> dict[str, Any]:
            spec_dict = dict()
            for token, value in spec_list[start_:stop_]:
                try:
                    col_name, dtype = record.inner_structure_mapping[token]
                except KeyError:
                    self._raise_or_warn(
                        f"Encountered invalid token in {record.name} records: {token}.",
                        raise_level=2
                    )
                    spec_dict[token] = value
                else:
                    spec_dict[col_name] = dtype.from_pdb(value)
            return spec_dict
        # Both COMPND and SOURCE records have only one field, called 'compound' and 'srcName', respectively.
        #  These fields are specification lists (see `opencadd.io.pdb.fields.SpecificationList`) with a given
        #  set of possible token-value pairs. The problem is that each token may exist several times in the
        #  list, once for each molecule, or fragment of a molecule, in the entry. What's worse is that there is
        #  no rules regarding the ordering of the tokens, except that the occurrence of 'MOL_ID' or 'FRAGMENT'
        #  tokens indicates that the subsequent tokens are related to that specific molecule or fragment of the
        #  molecule. Therefore, we first need to identify the positions of the MOL_ID tokens in the list, and
        #  then for each interval (i.e. between two MOL_ID tokens), indentify the positions of FRAGMENT tokens
        #  within that interval. Then, for each fragment of a molecule, a separate row in the dataframe is
        #  created, containing the general data of the molecule containing the fragment (i.e. all tokens
        #  between a MOL_ID token and the first occurrence of a FRAGMENT token before the next MOL_ID token),
        #  plus the data specific to that fragment (i.e. all tokens between the corresponding FRAGMENT token
        #  and either the next FRAGMENT token before the next MOL_ID, or the next MOL_ID token, in case for
        #  that molecule only one fragment has been entered).
        spec_list = self._extract_record(record)
        if spec_list is None:
            return
        ind_mols = np.argwhere(spec_list[:, 0] == "MOL_ID").reshape(-1)
        df = pd.DataFrame(columns=record.inner_structure_names)
        idx_row = 0
        for start, stop in np.lib.stride_tricks.sliding_window_view((*ind_mols, None), 2):
            ind_fragments = np.argwhere(spec_list[start:stop, 0] == "FRAGMENT").reshape(-1)
            if ind_fragments.size <= 1:
                df.loc[idx_row] = parse_token_value_pairs(start, stop)
                idx_row += 1
            else:
                data_mol = parse_token_value_pairs(start, start+ind_fragments[0])
                for start_frag, end_frag in np.lib.stride_tricks.sliding_window_view((*ind_fragments, None), 2):
                    data_frag = parse_token_value_pairs(start + start_frag, start + end_frag)
                    df.loc[idx_row] = data_mol | data_frag
                    idx_row += 1
        return df.set_index("mol_id", drop=False)

    @staticmethod
    def _parse_field_charge(fields: np.ndarray):
        charges = fields.view(dtype=(str, 2)).reshape(-1).astype((str, 3))
        charges[charges == ""] = "nan"
        return charges.astype(np.double)

    def _model_num_of_lines(self, idx_lines: np.ndarray):
        if not self.has_record(_records.MODEL):
            return np.ones(shape=idx_lines.size, dtype=np.uintc)
        idx_lines = np.expand_dims(idx_lines, axis=-1)
        indices_model_start = self.indices_record_lines(_records.MODEL)
        indices_model_end = self.indices_record_lines(_records.ENDMDL)
        mask_in_model = np.logical_and(idx_lines > indices_model_start, idx_lines < indices_model_end)
        counts = np.count_nonzero(mask_in_model, axis=1)
        line_in_no_model = counts == 0
        line_in_multiple_models = counts > 1
        if np.any(line_in_no_model):
            raise PDBParsingError("Line in no model")
        if np.any(line_in_multiple_models):
            raise PDBParsingError("line in multiple models.")
        ind_lines, ind_model = np.nonzero(mask_in_model)
        return ind_model + 1

    def mask_record(self, record: _records.Record) -> np.ndarray:
        """
        A boolean array indicating whether each line in the PDB file is of a given record type.

        Parameters
        ----------
        record

        Returns
        -------
        ndarray
        """
        return self._mask_records[:, record.index]

    def indices_record_lines(self, record: _records.Record) -> np.ndarray:
        return self._idx__record_lines[record.index]

    def count_record(self, record: _records.Record) -> int:
        return self._count_records[record.index]

    def has_record(self, record: _records.Record) -> bool:
        return self._has_record[record.index]

    def record_lines(self, record: _records.Record = None, as_char: bool = True, drop_name: bool = False):
        idx = self.indices_record_lines(record) if record is not None else slice(None)
        if as_char:
            return self._lines_chars[idx, (6 if drop_name else 0):]
        return self._lines_rec_content[idx] if drop_name else self._lines[idx]

    def _extract_record(
            self,
            record: _records.Record,
    ):
        if not self.has_record(record):
            self._raise_or_warn_if_mandatory(record)
            return
        record_lines = self.record_lines(record)
        return record.extract(record_char_table=record_lines)

    @property
    def count_lines(self) -> int:
        """Total number of lines in the PDB file."""
        return self._lines.size

    @property
    def strictness(self) -> Literal[0, 1, 2, 3]:
        """
        Level of strictness for raising exceptions and warnings when encountering mistakes
        in the PDB file:
            0: raise only fatal errors and don't show any warnings.
            1: raise only fatal errors. All other errors are reported as warnings.
            2: raise fatal errors and mistakes resulting in ambiguous data. Inconsequential mistakes are
            reported as warnings.
            3: completely validate the PDB file and raise all errors.

        Returns
        -------
        literal[0, 1, 2, 3]
        """
        return self._strictness

    @strictness.setter
    def strictness(self, value):
        if value not in (0, 1, 2, 3):
            raise ValueError(
                "Parameter `strictness` expects a value of either 0, 1, 2, or 3. "
                f"Input argument was: {value}."
            )
        self._strictness = value
        return

    def _raise_or_warn(
            self,
            msg: str,
            raise_level: Literal[1, 2, 3]
    ) -> NoReturn:
        """
        Given an error message, raise it as a `PDBParsingError` or post a warning message,
        depending on the level of `strictness` and the error level.

        Parameters
        ----------
        msg : str
            Error message.
        raise_level : {1, 2, 3}
            Minimum strictness level where the error should be raised as an exception.

        Raises
        -------
        PDBParsingError
        """
        if self.strictness == 0:
            return
        if self.strictness >= raise_level:
            raise PDBParsingError(msg)
        else:
            warnings.warn(msg)
        return

    def _raise_or_warn_if_mandatory(self, record: _records.Record):
        if record.is_mandatory:
            self._raise_or_warn(
                f"{record.name} record not found; it is a mandatory record.",
                raise_level=3
            )
        return

    def _validate_existence(self, record: _records.Record):
        if not self.has_record(record):
            self._raise_or_warn(
                f"{record.name} record not found. It is a mandatory record.",
                raise_level=3
            )
            return False
        return True

    def _validate_one_time_single_line_record(self, record: _records.Record) -> bool:
        if self.count_record(record) > 1:
            self._raise_or_warn(
                f"Multiple {record.name} records found at lines "
                f"{', '.join(self.indices_record_lines(record))}. "
                f"{record.name} is a one-time/single-line record and must only appear once in the file.",
                raise_level=1
            )
            return False
        return True

    def _validate_pdb_ids(self, pdb_ids: np.ndarray):
        validated = True
        for pdb_id in pdb_ids:
            if not pdb_id[0].isdigit():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. The first character in the ID must be a digit.",
                    raise_level=3
                )
                validated = False
            if not pdb_id[1:].isalpha():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. The last three characters in the ID must be alphabetical.",
                    raise_level=3
                )
                validated = False
            if not pdb_id.isupper():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. All characters in the ID must be upper-case.",
                    raise_level=3
                )
                validated = False
        return validated

    def validate_empty_columns(self, record: _records.Record):
        if self.has_record(record):
            non_spaces = self.record_lines(record)[:, record.empty_columns] != " "
            if np.any(non_spaces):
                positions = ", ".join([f"({x}, {y})" for x, y in np.argwhere(non_spaces)])
                self._raise_or_warn(
                    f"{record.name} records contain non-space characters in columns that must be empty."
                    f"Positions (line number, character number): {positions}",
                    raise_level=1
                )
                return False
        return True
