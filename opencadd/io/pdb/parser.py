import re
from typing import Literal, NoReturn, Optional, Dict, Any, Union
import warnings

import numpy as np
import pandas as pd

from opencadd.io import parsing
from opencadd.io._helpers_sysio import filelike_to_data_string
from opencadd.io.pdb.data import Record
from opencadd.io.pdb import records
from opencadd.io.pdb import sections
from opencadd.io.pdb import fields as pdbfields

class PDBParsingError(Exception):
    pass


class PDBFile:
    def __init__(
            self,
            title: sections.Title,
            remark: sections.Remark,
            primary_structure: sections.PrimaryStructure,
            heterogen: sections.Heterogen,
            secondary_structure: sections.SecondaryStructure,
            connectivity: sections.Connectivity,
            crystallographic: sections.Crystallographic,
            coordinate: sections.Coordinate,
            bookkeeping: sections.Bookkeeping,
    ):
        pass

    def to_pdb(self):
        def record_end() -> str:
            """
            END record of a PDB file.
            The END record marks the end of the PDB file, and must appear as the final record in every
            file.
            """
            return f"{'END':<80}"


class PDBParser:
    """
    rec_names : ndarray, shape: (n,), dtype: <U6
        1-dimensional array of strings, containing the record name of each line in the PDB file,
        stripped of whitespaces.
        For example: array(['HEADER', 'TITLE', 'TITLE', ..., 'CONECT', 'MASTER', 'END'])
    rec_chars : ndarray, shape: (n, m), dtype: <U1
        2-dimensional array of characters, where the first axis is the rows (lines) and second
        axis the columns (characters in a line) in a PDB file. In a standard PDB file,
        the number of columns 'm' must be 80.

    """
    def __init__(self, file, validate: bool = True):
        """
        Parameters
        ----------
        file : FileLike
            PDB file to parse.
        strictness : Literal[0, 1, 2, 3], optional, default: 3
            Level of strictness for raising exceptions and warnings when encountering mistakes
            in the PDB file:
                0: raise only fatal errors and don't show any warnings.
                1: raise only fatal errors.
                2: raise fatal errors and mistakes resulting in ambiguous data.
                3: completely validate the PDB file and raise all errors.
        """
        # Declare attributes
        self._validate: bool = validate
        self._lines: np.ndarray
        self._lines_chars: np.ndarray
        """Character view of lines, i.e. 2-d array where first axis is lines and second axis is 
        characters"""
        self._lines_rec_names: np.ndarray
        """Extract all record types as 1-d array of stripped record-name strings."""
        self._lines_rec_content: np.ndarray
        self._mask_record: dict
        self._idx_record: dict
        self._has_record: dict
        self._count_record: dict

        # Read file and set attributes
        self.validate = validate
        self._lines = np.array(filelike_to_data_string(file=file).splitlines())
        self._lines_chars = self._lines.view(dtype=(str, 1)).reshape(self._lines.size, -1)
        self._lines_rec_names = np.char.strip(self._lines_chars[:, :6].view(dtype=(str, 6)).reshape(-1))
        self._lines_rec_content = self._lines_chars[:, 6:].view(dtype=(str, 74)).reshape(-1)
        self._mask_record = dict()
        self._idx_record = dict()
        self._has_record = dict()
        self._count_record = dict()
        return

    @property
    def file(self) -> PDBFile:
        pdb_file = PDBFile(
            title=self.section_title,
            remark=self.section_remark,
            primary_structure=self.section_primary_structure,
            heterogen=self.section_heterogen,
            secondary_structure=self.section_secondary_structure,
            connectivity=self.section_connectivity,
            feature=self.section_feature,
            crystallographic=self.section_crystallographic,
            coordinate=self.section_coordinate,
            bookkeeping=self.section_bookkeeping,
        )
        return pdb_file

    @property
    def section_title(self) -> sections.Title:
        title = sections.Title(
            header=self.header,
            title=self.title,
            compound=self.compnd,
        )
        return title

    @property
    def section_remark(self) -> sections.Remark:
        pass

    @property
    def section_primary_structure(self) -> sections.PrimaryStructure:
        primary_structure = sections.PrimaryStructure(
            dbref=self.dbref,
            seqadv=self.seqadv,
            seqres=self.seqres,
            modres=self.modres
        )
        return primary_structure

    @property
    def section_heterogen(self) -> sections.Heterogen:
        pass

    def section_secondary_structure(self) -> sections.SecondaryStructure:
        pass

    def section_connectivity_annotation(self) -> sections.Connectivity:
        return

    def section_feature(self):
        pass

    def section_crystallographic(self) -> sections.Crystallographic:
        pass

    def section_coordinate_transformation(self):
        pass

    def section_coordinate(self):
        if not (self.has_record(Record.ATOM) or self.has_record(Record.HETATM)):
            return None
        df_ter_records = self.ter
        mask_lines = np.logical_or(
            self.mask_record(Record.ATOM), self.mask_record(Record.HETATM)
        )
        if not self.has_record(Record.MODEL):
            lines = self._lines_chars[mask_lines]
            atom_data = parsing.extract_columns(lines, Record.ATOM.columns)
            atom_data["model_nr"] = np.full(shape=lines.shape[0], fill_value=1, dtype=np.ubyte)
            atom_data["is_standard_residue"] = self._lines_rec_names[mask_lines] == "ATOM"
            atom_data["is_polymer"] = atom_data["serial"] < df_ter_records.iloc[-1]["serial"]
            return pd.DataFrame(atom_data)

        num_coordinates = np.count_nonzero(mask_lines)

        if self.has_record(Record.MODEL):
            indices_model_start = self.idx_record(Record.MODEL)
            indices_model_end = self.idx_record(Record.ENDMDL)
            for idx_start, idx_end in zip(indices_model_start, indices_model_end):
                pass


        # df_anisou = parsing.dataframe_from_string_array(
        #     self._rec_chars[self._mask_record[Record.ANISOU]],
        #     column_data=Record.ANISOU.columns
        # )

        df_atom_data = parsing.dataframe_from_string_array(
            self._lines_chars[mask_lines],
            column_data=Record.ATOM.columns
        )
        is_polymer = np.zeros(shape=df_atom_data.index.values.size, dtype=np.bool_)
        if self.has_record(Record.MODEL):
            size_models = np.array(
                [
                    np.count_nonzero(mask_lines[start:end])
                    for start, end in zip(
                        self.idx_record(Record.MODEL), self.idx_record(Record.ENDMDL)
                    )
                ]
            )
            df_atom_data["model_nr"] = np.repeat(
                np.arange(1, size_models.size + 1), repeats=size_models
            )
            if self.has_record(Record.TER):
                idx_ter = self.idx_record(Record.TER)
                for idx_model in range(size_models.size):
                    model_start_line = self.idx_record(Record.MODEL)[idx_model]
                    model_end_line = self.idx_record(Record.ENDMDL)[idx_model]
                    ter_mask = np.logical_and(
                        idx_ter > model_start_line,
                        idx_ter < model_end_line
                    )
                    num_polymer_atoms = np.count_nonzero(
                        mask_lines[model_start_line:idx_ter[ter_mask][-1]]
                    )
                    start_idx = np.sum(size_models[:idx_model])
                    is_polymer[start_idx:start_idx+num_polymer_atoms] = True
        else:
            df_atom_data["model_nr"] = 1
            if self.has_record(Record.TER):
                idx_last_ter = self.idx_record(Record.TER)[-1]
                is_polymer[: -np.count_nonzero(mask_lines[idx_last_ter + 1:])] = True
        df_atom_data["is_polymer"] = is_polymer
        return df_atom_data

    def section_connectivity(self):
        pass

    def section_bookkeeping(self):
        pass

    @property
    def header(self) -> Optional[records.Header]:
        """
        Parse the HEADER record of the PDB file.

        Returns
        -------
        """
        if not self.has_record(Record.HEADER):
            return None
        fields = self.extract_record_fields(Record.HEADER)
        return records.Header(
            classification=fields["classification"],
            deposition_date=fields["deposition_date"][0],
            pdb_id=fields["pdb_id"][0]
        )

    @property
    def obslte(self) -> Optional[records.Obsolete]:
        return

    @property
    def title(self) -> Optional[str]:
        if not self.has_record(Record.TITLE):
            return None
        title = parsing.extract_column(
            self.record_lines(Record.TITLE), Record.TITLE.columns["title"]
        )
        return title

    @property
    def split(self) -> Optional[records.Split]:
        return

    @property
    def caveat(self) -> Optional[records.Caveat]:
        return

    @property
    def compnd(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.COMPND):
            return None
        return self.extract_record_fields(Record.COMPND, to_dataframe=True)

    @property
    def source(self):
        return

    @property
    def keywds(self):
        return

    @property
    def expdata(self):
        return

    @property
    def nummdl(self):
        return

    @property
    def mdltyp(self):
        return

    @property
    def author(self):
        return

    @property
    def revdat(self):
        return

    @property
    def sprsde(self):
        return

    @property
    def jrnl(self):
        return

    @property
    def dbref(self) -> Optional[pd.DataFrame]:
        if not (self.has_record(Record.DBREF) or self.has_record(Record.DBREF1)):
            return None
        df_main = pd.DataFrame(columns=Record.DBREF.columns.keys())
        if self.has_record(Record.DBREF):
            df_main = pd.concat(
                (df_main, self.extract_record_fields(Record.DBREF, to_dataframe=True))
            )
        if self.has_record(Record.DBREF1):
            df = pd.concat(
                [
                    self.extract_record_fields(record, to_dataframe=True)
                    for record in (Record.DBREF1, Record.DBREF2)
                ],
                axis=1
            )
            df_main = pd.concat((df_main, df))
        return df_main

    @property
    def seqadv(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.SEQADV):
            return None
        return self.extract_record_fields(Record.SEQADV, to_dataframe=True)

    @property
    def seqres(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.SEQRES):
            return None
        fields = self.extract_record_fields(Record.SEQRES)
        unique_chain_ids, idx_unique_chain_ids = np.unique(fields["chain_id"], return_index=True)
        res_names = np.empty(unique_chain_ids.size, dtype=object)
        for i, unique_chain_id in enumerate(unique_chain_ids):
            mask_chain = fields["chain_id"] == unique_chain_id
            res_names_pre = fields["residue_name"][mask_chain][
                np.argsort(fields["serial"][mask_chain])
            ]
            res_names[i] = res_names_pre[res_names_pre != ""]
        df = pd.DataFrame(
            {
                "residue_count": fields["residue_count"][idx_unique_chain_ids],
                "residue_name": res_names
            }, index=unique_chain_ids
        )
        df.index.name = "chain_id"
        return df

    @property
    def modres(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.MODRES):
            return
        return self.extract_record_fields(Record.MODRES, to_dataframe=True)

    @property
    def het(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HET):
            return None
        return self.extract_record_fields(Record.HET, to_dataframe=True)

    @property
    def hetnam(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HETNAM):
            return None
        fields = self.extract_record_fields(Record.HETNAM)
        _, idx_unique_ids = np.unique(fields["id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["id"][idx_unique_ids]
        names = list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["id"] == unique_id
            names_pre = fields["name"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            names.append(pdbfields.String.from_pdb(names_pre))
        df = pd.DataFrame(
            {
                "name": names
            }, index=unique_ids
        )
        df.index.name = "hetID"
        return df

    @property
    def hetsyn(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HETSYN):
            return None
        fields = self.extract_record_fields(Record.HETSYN)
        _, idx_unique_ids = np.unique(fields["id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["id"][idx_unique_ids]
        synonyms = list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["id"] == unique_id
            synonyms_pre = fields["synonyms"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            synonyms.append(pdbfields.SList.from_pdb(pdbfields.String.from_pdb(synonyms_pre)))
        df = pd.DataFrame(
            {
                "synonyms": synonyms
            }, index=unique_ids
        )
        df.index.name = "hetID"
        return df

    @property
    def formul(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.FORMUL):
            return None
        fields = self.extract_record_fields(Record.FORMUL)
        _, idx_unique_ids = np.unique(fields["id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["id"][idx_unique_ids]
        formula = list()
        num_in_chain = list()
        num_out_chain = list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["id"] == unique_id
            formula_pre = fields["formula"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            formula_single_line = pdbfields.String.from_pdb(formula_pre)
            formula_split = re.split(r'[()]+', formula_single_line)
            num_split = len(formula_split)
            if num_split == 1:
                formula.append(formula_split[0])
                num_in_chain.append(0)
                num_out_chain.append(0)
            elif num_split == 3:
                formula.append(formula_split[-2])
                num_in_chain.append(formula_split[0])
                num_out_chain.append(0)
            elif num_split == 4:
                formula.append(formula_split[-2])
                num_in_chain.append(formula_split[1])
                num_out_chain.append(formula_split[0])
            else:
                raise PDBParsingError(
                    f"Could not parse FORMUL record's field {formula_single_line}."
                )
        df = pd.DataFrame(
            {
                "component_num": fields["component_num"][idx_unique_ids],
                "formula": formula,
                "count_occurrence_within_chain": num_in_chain,
                "count_oocurrence_outside_chain": num_out_chain,
            }, index=unique_ids
        )
        df.index.name = "hetID"
        return df

    @property
    def helix(self):
        return

    @property
    def sheet(self):
        return

    @property
    def ssbond(self):
        return

    @property
    def link(self):
        return

    @property
    def cispep(self):
        return

    @property
    def site(self):
        return

    @property
    def cryst1(self):
        return

    @property
    def origx(self):
        return

    @property
    def scale(self):
        return

    @property
    def mtrix(self):
        return

    @property
    def atom(self):
        return

    @property
    def anisou(self) -> Optional[pd.DataFrame]:
        return self._parse_record_in_model(Record.ANISOU)

    @property
    def ter(self) -> Optional[pd.DataFrame]:
        return self._parse_record_in_model(Record.TER)

    @property
    def hetatm(self):
        return

    @property
    def conect(self):
        return

    @property
    def master(self):
        return

    def _parse_record_in_model(self, record: Record) -> Optional[pd.DataFrame]:
        if not self.has_record(record):
            return
        fields = self.extract_record_fields(record)
        if not self.has_record(Record.MODEL):
            model_nr = np.full(shape=self.count_record(record), fill_value=1, dtype=np.ubyte)
        else:
            count_record_per_model = [
                np.count_nonzero(
                    np.logical_and(
                        self.idx_record(record) > model_start_idx,
                        self.idx_record(record) < model_end_idx,
                    )
                ) for model_start_idx, model_end_idx in zip(
                    self.idx_record(Record.MODEL), self.idx_record(Record.ENDMDL)
                )
            ]
            model_nr = np.repeat(
                np.arange(1, self.count_record(Record.MODEL) + 1), repeats=count_record_per_model
            )
        return pd.DataFrame({"model_nr": model_nr, **fields})

    @property
    def master(self) -> sections.Bookkeeping:
        """
        Parse the MASTER record in a PDB file.

        Returns
        -------
        BookkeepingSection

        Raises
        ------
        PDBParsingError

        Notes
        -----
        MASTER record is a mandatory one-time/single-line record that must appear only once, as the
        second-last record in every PDB file. It is a control record, counting the occurrence
        (checksum) of selected record types in the PDB file.
        """
        # Count relevant records to construct the MASTER record data
        bookkeeping_section = pdbstruct.BookkeepingSection(
            remark=np.count_nonzero(self._lines_rec_names == "REMARK"),
            het=np.count_nonzero(self._lines_rec_names == "REMARK"),
            helix=np.count_nonzero(self._lines_rec_names == "HELIX"),
            sheet=np.count_nonzero(self._lines_rec_names == "SHEET"),
            site=np.count_nonzero(self._lines_rec_names == "SITE"),
            xform=np.count_nonzero(
                np.isin(
                    self._lines_rec_names,
                    [
                        "ORIGX1",
                        "ORIGX2",
                        "ORIGX3",
                        "SCALE1",
                        "SCALE2",
                        "SCALE3",
                        "MTRIX1",
                        "MTRIX2",
                        "MTRIX3",
                    ]
                )
            ),
            coord=np.count_nonzero(np.isin(self._lines_rec_names, ["ATOM", "HETATM"])),
            ter=np.count_nonzero(self._lines_rec_names == "TER"),
            connect=np.count_nonzero(self._lines_rec_names == "CONECT"),
            seqres=np.count_nonzero(self._lines_rec_names == "SEQRES"),
        )
        # Look for existing MASTER records in file
        idx_master_record = np.argwhere(self._lines_rec_names == "MASTER").rehape(-1)
        # Case 1: no existing MASTER record in file
        if idx_master_record.size == 0:
            self._raise_or_warn(
                "MASTER record not found. It is a mandatory record and must appear once in "
                "the file, as the second-last record (before the END record).",
                raise_level=3,
            )
            return bookkeeping_section

        if idx_master_record.size > 1:
            self._raise_or_warn(
                f"Multiple MASTER records found ({idx_master_record.size} at lines "
                f"{idx_master_record + 1}). It is a one-time/single-line record and must appear "
                "only once in the file.",
                raise_level=2,
            )
        # Parse existing MASTER record
        master = parsing.extract_columns(self._lines_chars[idx_master_record[-1]],
                                         Record.MASTER.columns)
        master_record = pdbstruct.BookkeepingSection(*master)

        if master_record != bookkeeping_section:
            self._raise_or_warn(
                "MASTER record values do not match with number of records in file. "
                f"MASTER record values are {master_record}, "
                f"but actual record counts are {bookkeeping_section}.",
                raise_level=2,
            )

        if idx_master_record != self._lines_rec_names.size - 2:
            self._raise_or_warn(
                f"MASTER record found in wrong position (line {idx_master_record[0] + 1}). "
                f"It must be the second-last record in the file, "
                f"i.e. at line {self._lines_rec_names.size - 2}.",
                raise_level=3,
            )
        empty_cols = np.concatenate((np.arange(7, 11), np.arange(71, 81))) - 1
        if np.any(self._lines_chars[idx_master_record[0], empty_cols] != " "):
            self._raise_or_warn(
                "MASTER record line contains characters in columns that must be empty (i.e. "
                f"columns 7–10 and 71–80). "
                f"Found: {self._lines_chars[idx_master_record[0], empty_cols]}",
                raise_level=3,
            )
        return bookkeeping_section

    def validate_record_header(self):
        ind_lines = self._idx_record[Record.HEADER]
        # Check existence; HEADER is a mandatory record.
        if not self._has_record[Record.HEADER]:
            self._raise_or_warn(
                "HEADER record not found. It is a mandatory record.",
                raise_level=3
            )
            # If there is no HEADER record, return empty header object
            return pdbstruct.HEADER()
        # Check number of occurrences; HEADER is a one-time/single-line record.
        if ind_lines.size != 1:
            self._raise_or_warn(
                "Multiple HEADER records found at lines "
                f"{', '.join(ind_lines)}. "
                "HEADER is a one-time/single-line record and must only appear once in the file.",
                raise_level=3
            )
            # If there are multiple HEADER records, return header object containing all
            # unprocessed HEADER lines marked as invalid liens.
            return pdbstruct.HEADER(
                invalid_lines=self._lines_rec_content[ind_lines]
            )
        # Check position in file; HEADER must always be in first line:
        if ind_lines[0] != 0:
            self._raise_or_warn(
                f"HEADER record found at line {ind_lines[0]}. "
                f"It must always be the first line of a PDB file.",
                raise_level=3
            )
        if self._check_empty_columns(Record.HEADER):
            # Record is in standard format
            return pdbstruct.HEADER(
                classification=parsing.extract_column_by_Column(self._record_lines[Record.HEADER], )
            )
        else:
            # Not in standard format
            return pdbstruct.HEADER(
                invalid_lines=self._lines_rec_content[ind_lines]
            )

    def validate_record_seqres(self):
        # If no SEQRES record exists
        if not self.has_record(Record.SEQRES):
            # but ATOM records exist
            if self.has_record(Record.ATOM):
                self._raise_or_warn(
                    "SEQRES records not found, while ATOM records exist. "
                    "SEQRES is a mandatory record when ATOM records exist.",
                    raise_level=3
                )
        self._check_routine(Record.SEQRES)
        return

    def validate_line_width(self):
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
        return

    def validate_record_names(self):
        if not np.all(np.isin(self._lines_rec_names, Record.names)):
            self._raise_or_warn(
                f"Unrecognized record names found.",
                raise_level=3
            )

    def validate_section_coordinate(self):

        if ind_ter.size == 0:
            if ind_atom.size != 0:
                self._raise_or_warn(
                    "ATOM records found without TER record.",
                    raise_level=1,
                )

        indices_record_model = self._idx_record[Record.MODEL]  # line indices of MODEL records
        count_model = indices_record_model.size  # count of MODEL records
        if count_model > 0:
            indices_record_endmdl = self._idx_record[Record.ENDMDL]
            count_endmdl = indices_record_endmdl.size
            if count_model != count_endmdl:
                self._raise_or_warn(
                    "Number of MODEL and ENDMDL records do not match: "
                    f"{count_model} vs {count_endmdl}",
                    raise_level=1
                )
            model_lengths = indices_record_endmdl - indices_record_model
            if model_lengths[0] < 1 or np.any(model_lengths != model_lengths[0]):
                self._raise_or_warn(
                    "Models have different lengths or are empty",
                    raise_level=1
                )

        else:
            # If there is no MODEL record, then there is only one model
            # Extract the records of the single model
            filter_model = np.isin(self._lines_rec_names, ["ATOM", "ANISOU", "TER", "HETATM"])

            model = read_model(rec_names[filter_model], rec_chars[filter_model])
            return CoordinateSection(models=(model,))

        models = []
        for start_idx, end_idx in zip(models_start_ind, models_end_ind):
            models.append(
                read_model(
                    record_types=rec_names[start_idx + 1:end_idx],
                    pdb_lines=rec_chars[start_idx + 1:end_idx]
                )
            )
        return CoordinateSection(models=tuple(models))


    def validate_record_dbref(self):
        if not (
                self._has_record[Record.DBREF]
                or self._has_record[Record.DBREF1]
                or self._has_record[Record.DBREF2]
        ):
            if self._has_record[Record.ATOM] or self._has_record[Record.TER]:
                self._raise_or_warn(
                    "DBREF and DBREF1/DBREF2 records not found, while ATOM records exist. "
                    "At least one must exist when ATOM records exist.",
                    raise_level=3
                )
        if self._has_record[Record.DBREF1]:
            if self._idx_record[Record.DBREF1].size != self._idx_record[Record.DBREF2].size:
                self._raise_or_warn(
                    "Number of DBREF1 and DBREF2 records do not match. DBREF1/DBREF2 records are "
                    "two-line records that always have to be present together.",
                    raise_level=2
                )
            self._check_empty_columns(Record.DBREF1)
            self._check_empty_columns(Record.DBREF2)
        self._check_routine(Record.DBREF)
        return

    def validate_record_end(self) -> NoReturn:
        """
        Parse the END record in a PDB file.

        Returns
        -------
        None
            Since the END record contains no information, nothing is returned. Only errors are
            raised when applicable.

        Raises
        ------
        PDBParsingError

        Notes
        -----
        END record is a mandatory one-time/single-line record that must appear only once,
        as the final record at the end of every PDB file.
        The first 3 characters of the line must be 'END', followed by 77 spaces to make
        an 80-character line. Since it contains no information, it can be completely ignored when
        PDB file-format validation is not a concern.
        """
        idx_end_record = np.argwhere(self._lines_rec_names == "END").rehape(-1)
        if idx_end_record.size == 0:
            self._raise_or_warn(
                "END record not found. It is a mandatory record and must appear once at "
                "the end of the file.",
                raise_level=3,
            )
        if idx_end_record.size > 1:
            self._raise_or_warn(
                f"Multiple END records found ({idx_end_record.size} at lines "
                f"{idx_end_record + 1}). It is a one-time/single-line record and must appear "
                "only once in the file.",
                raise_level=3,
            )
        if idx_end_record != self._lines_rec_names.size - 1:
            self._raise_or_warn(
                f"END record found at the middle of the file (line {idx_end_record[0] + 1}). "
                "It must be the final record in the file.",
                raise_level=3,
            )
        if np.any(self._lines_chars[-1, 6:] != " "):
            self._raise_or_warn(
                "END record line contains extra characters. It must contain only the record name "
                "'END' at beginning of the line. "
                f"Found: {self._lines_chars[-1, 6:].view(dtype=(str, self._lines_chars.shape[1] - 6)).reshape(-1)}",
                raise_level=3,
            )
        return

    def mask_record(self, record: Record) -> np.ndarray:
        """
        A boolean array indicating whether each line in the PDB file is of a given record type.

        Parameters
        ----------
        record

        Returns
        -------
        ndarray
        """
        mask = self._mask_record.get(record)
        if mask is not None:
            return mask
        mask = self._lines_rec_names == record.name
        self._mask_record[record] = mask
        return mask

    def idx_record(self, record: Record) -> np.ndarray:
        idx = self._idx_record.get(record)
        if idx is not None:
            return idx
        idx = np.argwhere(self.mask_record(record)).reshape(-1)
        self._idx_record[record] = idx
        return idx

    def count_record(self, record: Record) -> int:
        count = self._count_record.get(record)
        if count is not None:
            return count
        count = self.idx_record(record).size
        self._count_record[record] = count
        return count

    def has_record(self, record: Record) -> bool:
        has = self._has_record.get(record)
        if has is not None:
            return has
        has = self.count_record(record) != 0
        self._has_record[record] = has
        return has

    def record_lines(self, record: Record, as_char: bool = True, drop_name: bool = False):
        idx = self.idx_record(record)
        if as_char:
            return self._lines_chars[idx, (6 if drop_name else 0):]
        return self._lines_rec_content[idx] if drop_name else self._lines[idx]

    def extract_record_fields(
            self,
            record: Record,
            to_dataframe: bool = False
    ) -> Union[Dict[str, Any], pd.DataFrame]:
        lines = self.record_lines(record=record)
        return parsing.extract_columns(
            char_array=lines, columns=record.columns, to_dataframe=to_dataframe
        )

    @property
    def validate(self):
        return self._validate

    @validate.setter
    def validate(self, value):
        """Raise a ValueError if the input argument `validate` is not a valid input."""
        if not isinstance(value, bool):
            raise TypeError(
                "Parameter `validate` expects a boolean value. "
                f"Input argument {value} had type {type(value)}."
            )
        self._validate = value
        return


    def _check_routine(self, record: Record):
        if self.has_record(record):
            self._check_records_are_consecutive(record)
            self._check_empty_columns(record)

    def _check_empty_columns(self, record: Record):
        if self._has_record[record] and np.any(
                self.record_lines(record)[:, record.empty_columns] != " "
        ):
            self._raise_or_warn(
                f"{record.name} records contain non-space characters in columns that must be "
                f"empty.",
                raise_level=3
            )
            return False
        return True

    def _check_records_are_consecutive(self, record: Record):
        if self._has_record[record] and np.any(np.diff(self._idx_record[record]) != 1):
            self._raise_or_warn(
                f"{record.name} records are not in consecutive lines.",
                raise_level=3
            )
        return

    def _raise_or_warn(
            self,
            msg: str,
            raise_level: Literal[1, 2, 3]
    ) -> NoReturn:
        if self.strictness >= raise_level:
            raise PDBParsingError(msg)
        else:
            warnings.warn(msg)
        return

    def _parse_record_as_df(self, record: Record):
        return parsing.dataframe_from_string_array(
            lines=self._record_lines[record],
            column_data=record.columns
        )

parser = PDBParser("/Users/home/Downloads/3w32as.pdb")
x= parser.formul