from typing import NamedTuple, Literal, NoReturn
import warnings
import datetime

import numpy as np
import pandas as pd

from opencadd.io._helpers_parsing import extract_column_from_string_array
from opencadd.io.pdb import datastruct as pdbstruct
from opencadd.io._helpers_sysio import filelike_to_data_string

class PDBParsingError(Exception):
    pass


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
    def __init__(self, file, strictness: Literal[0, 1, 2, 3] = 3):
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
        self._strictness: Literal[0, 1, 2, 3]
        self._pdb_lines: np.ndarray
        self._line_widths: np.ndarray
        self._rec_chars: np.ndarray
        self._rec_names: np.ndarray

        # Read file and set attributes
        data_str = filelike_to_data_string(file=file)
        self._pdb_lines = np.array(data_str.splitlines())
        self.strictness = strictness
        # Check if all lines are 80 characters:
        self._line_widths = np.char.str_len(self._pdb_lines)
        lines_not80 = self._line_widths != 80
        if np.any(lines_not80):
            faulty_lines = [
                f"{line_num} ({line_width})" for line_num, line_width in zip(
                    np.argwhere(lines_not80).reshape(-1) + 1,
                    self._line_widths[lines_not80]
                )
            ]
            self._raise_or_warn(
                "Following lines are not 80 characters long (lengths given in brackets): "
                f"{', '.join(faulty_lines)}",
                raise_level=3
            )
        # Get a character view of lines, i.e. 2-d array where first axis is lines and second
        # axis is characters:
        self._rec_chars = self._pdb_lines.view(dtype=(str, 1)).reshape(self._pdb_lines.size, -1)
        # Extract all record types as 1-d array of stripped record-name strings:
        self._rec_names = np.char.strip(self._rec_chars[:, :6].view(dtype=(str, 6)).reshape(-1))
        if not np.all(np.isin(self._rec_names, ))
        return

    @property
    def strictness(self):
        return self._strictness

    @strictness.setter
    def strictness(self, value):
        """Raise a ValueError if the input argument `strictness` is not a valid input."""
        if value not in (1, 2, 3):
            raise ValueError(
                "Parameter `strictness` expects a literal value of either 1, 2, or 3. "
                f"Input argument was: {value}."
            )
        self._strictness = value
        return

    def file(self) -> pdbstruct.PDBFile:
        pdb_file = pdbstruct.PDBFile(
            title=self.section_title,
            remark=self.section_remark,
            primary_structure=self.section_primary_structure,
            heterogen=self.section_heterogen,
            secondary_structure=self.section_secondary_structure,
            connectivity=self.section_connectivity,
            site=self.record_site,
            crystallographic=self.section_crystallographic,
            coordinate=self.section_coordinate,
            bookkeeping=self.section_bookkeeping,
        )
        # Intersectional sanity check
        return pdb_file

    @property
    def section_title(self) -> pdbstruct.TitleSection:
        title = pdbstruct.TitleSection(
            header=self.record_header,
        )
        return title

    @property
    def section_remark(self) -> pdbstruct.RemarkSection:
        pass

    @property
    def section_primary_structure(self) -> pdbstruct.PrimaryStructureSection:
        pass

    @property
    def section_heterogen(self) -> pdbstruct.HeterogenSection:
        pass

    @property
    def section_secondary_structure(self) -> pdbstruct.SecondaryStructureSection:
        pass

    @property
    def section_connectivity(self) -> pdbstruct.ConnectivitySection:

        return

    @property
    def section_crystallographic(self) -> pdbstruct.CrystallographicSection:
        pass

    @property
    def section_coordinate(self):
        # Check how many MODELs there are in the file
        filter_models = self._rec_names == "MODEL"
        models_start_ind = np.argwhere(filter_models).reshape(-1)
        # If there are models present:
        if models_start_ind.size > 0:
            models_end_ind = np.argwhere(self._rec_names == "ENDMDL").reshape(-1)
            if models_start_ind.size != models_end_ind.size:
                self._raise_or_warn(
                    "Number of MODEL and ENDMDL records do not match: "
                    f"{models_start_ind.size} vs {models_end_ind.size}",
                    raise_level=1
                )
            model_lengths = models_end_ind - models_start_ind
            if model_lengths[0] < 1 or np.any(model_lengths != model_lengths[0]):
                self._raise_or_warn(
                    "Models have different lengths or are empty",
                    raise_level=1
                )

        # If there is no MODEL record, then there is only one model
        if models_start_ind.size == 0:
            # Extract the records of the single model


            filter_model = np.isin(rec_names, ["ATOM", "ANISOU", "TER", "HETATM"])


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


    def _parse_coordinate_no_model(self):
        mask_atom = self._rec_names == "ATOM"
        mask_hetatm = self._rec_names == "HETATM"
        mask_anisou = self._rec_names == "ANISOU"
        mask_ter = self._rec_names == "TER"
        ind_atom = np.argwhere(mask_atom).reshape(-1)
        ind_hetatm = np.argwhere(mask_hetatm).reshape(-1)
        ind_anisou = np.argwhere(mask_anisou).reshape(-1)
        ind_ter = np.argwhere(mask_ter).reshape(-1)

        if ind_ter.size == 0:
            if ind_atom.size != 0:
                self._raise_or_warn(
                    "ATOM records found without TER record.",
                    raise_level=1,
                )




    def _build_coordinate(self, idx_line_start: int = 0, idx_line_end: int = None):
        recs = self._rec_names[idx_line_start:idx_line_end]
        chars = self._rec_chars[idx_line_start:idx_line_end]
        return

    def _read_model(self, idx_line_model: int, idx_line_endmdl: int) -> Model:
        has_atom = np.isin("ATOM", record_types)
        has_hetatm = np.isin("HETATM", record_types)
        has_ter = np.isin("TER", record_types)
        atom = read_atom_record(pdb_lines[record_types == "ATOM"]) if has_atom else None
        hetatom = read_atom_record(pdb_lines[record_types == "HETATM"]) if has_hetatm else None
        ter = read_ter_record(pdb_lines[record_types == "TER"]) if has_ter else None
        return Model(atom=atom, hetatm=hetatom)

    @property
    def section_bookkeeping(self) -> pdbstruct.BookkeepingSection:
        self.record_end
        return self.record_master

    @property
    def record_header(self) -> pdbstruct.HeaderRecord:
        """
        Parse the HEADER records of a PDB file.

        Returns
        -------
        HeaderRecord
            A NamedTuple with keys `pdb_id`, `classification`, `deposition_date` and `invalid_lines`.
        """
        # HEADER record specifications
        columns = {
            "classification": ((11, 50), (str, 40)),
            "depDate": ((51, 59), (str, 9)),
            "idCode": ((63, 66), (str, 4)),
        }
        columns_empty = [(7, 10), (60, 62), (67, 80)]
        # First record type must always be HEADER:
        if self._rec_names[0] != "HEADER":
            warnings.warn(f"First line is not a HEADER record, instead: {pdb_lines[0]}")
        # HEADER is a one-time/single-line record:
        is_header = self._rec_names == "HEADER"
        header_lines = self._pdb_lines[is_header]
        header_is_one_time = header_lines.size == 1
        if not header_is_one_time:
            self._raise_or_warn(
                f"There must be only one HEADER record in the file, but found {header_lines.size} "
                f"at lines {', '.join(np.argwhere(is_header).reshape(-1))}: {header_lines}",
                raise_level=3
            )
        # Check whether it's a standard HEADER, i.e. only one line with 80 characters in standard
        # format (check if it's in standard format by checking if the white spaces are conserved):
        if not (
            header_is_one_time
            and len(header_lines[0]) == 80
            and np.all(
                [
                    np.char.isspace(
                        extract_column_from_string_array(
                            array=header_lines,
                            char_range=(start_col - 1, end_col)
                        )
                    ) for (start_col, end_col) in columns_empty
                ]
            )
        ):
            # Not in standard format; return raw data:
            warnings.warn(f"HEADER data is not in standard format: {header_lines}")
            return pdbstruct.HeaderRecord(
                invalid_lines=tuple(header_lines)
            )
        # In standard format; extract data:
        classification, dep_date, pdb_id = (
            np.char.strip(
                extract_column_from_string_array(
                    array=header_lines,
                    char_range=(col_range[0] - 1, col_range[1])
                )
            ).astype(col_dtype)[0]
            for col_name, (col_range, col_dtype) in columns.items()
        )
        # Classification string can have '/' characters separating the classification for each
        # entity in the file, and for each entity its different functions may be separated by ','
        # characters:
        class_per_entity = tuple(classification.split("/"))
        class_per_entity_and_function = tuple(
            tuple(entity_function.strip() for entity_function in entity_class.split(","))
            for entity_class in class_per_entity
        )
        # Create Header record
        header = pdbstruct.HeaderRecord(
            pdb_id=pdb_id,
            classification=class_per_entity_and_function,
            deposition_date=datetime.datetime.strptime(dep_date, "%d-%b-%y").date()
        )
        return header

    @property
    def record_site(self) -> pdbstruct.SiteSection:
        pass

    @property
    def record_ter(self):
        columns = {
            "serial": ((7, 11), int),
            "resName": ((18, 20), (str, 3)),
            "chainID": ((22, 22), (str, 1)),
            "resSeq": ((23, 26), int),
            "iCode": ((27, 27), (str, 1)),
        }
        df = pd.DataFrame
        for col_name, (col_range, col_dtype) in columns.items():
            df[col_name] = np.char.strip(
                extract_column_from_string_array(
                    array=pdb_lines,
                    char_range=(col_range[0] - 1, col_range[1])
                )
            ).astype(col_dtype)
        return df

    def record_atom(self) -> pd.DataFrame:
        """
        Parse ATOM records
        """
        columns = {
            "serial": ((7, 11), int),
            "name": ((13, 16), (str, 4)),
            "altLoc": ((17, 17), (str, 1)),
            "resName": ((18, 20), (str, 3)),
            "chainID": ((22, 22), (str, 1)),
            "resSeq": ((23, 26), int),
            "iCode": ((27, 27), (str, 1)),
            "x": ((31, 38), float),
            "y": ((39, 46), float),
            "z": ((47, 54), float),
            "occupancy": ((55, 60), float),
            "tempFactor": ((61, 66), float),
            "element": ((77, 78), (str, 2)),
            "charge": ((79, 80), (str, 2))
        }
        df = pd.DataFrame()
        for col_name, (col_range, col_dtype) in columns.items():
            df[col_name] = np.char.strip(
                extract_column_from_string_array(
                    array=pdb_lines,
                    char_range=(col_range[0] - 1, col_range[1])
                )
            ).astype(col_dtype)
        # Parse the charge column and turn values into float
        # Values, if present, are composed of a digit followed by a '+' or '-' sign.
        # Pandas recasts all string data types into 'object'. Thus, first change back to '<U2' and
        # then take a per-character view (cannot directly be done from object type)
        charges = df.charge.values.astype("<U2").view(dtype=(str, 1)).reshape(
            pdb_lines.size, -1
        )
        # Swap first and second columns to put the +/- sign in front of digit
        charges[:, [0, 1]] = charges[:, [1, 0]]
        # Combine again
        charges = np.frombuffer(charges.tobytes(), dtype=(str, 2))
        # Replace empty strings with 'nan'
        charges_cleaned = np.where(charges != "", charges, "nan")
        # Turn to float and assign back to dataframe
        df.charge = charges_cleaned.astype(np.single)
        # Replace empty strings with None
        df.replace(r'^$', None, regex=True, inplace=True)
        return df

    @property
    def record_connect(self):
        columns = {
            "serial_ref": ((7, 11), int),
            "serial_b1": ((12, 16), int),
            "serial_b2": ((17, 21), int),
            "serial_b3": ((22, 26), int),
            "serial_b4": ((27, 31), int),
        }


    @property
    def record_master(self) -> pdbstruct.BookkeepingSection:
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
                remark=np.count_nonzero(self._rec_names == "REMARK"),
                het=np.count_nonzero(self._rec_names == "REMARK"),
                helix=np.count_nonzero(self._rec_names == "HELIX"),
                sheet=np.count_nonzero(self._rec_names == "SHEET"),
                site=np.count_nonzero(self._rec_names == "SITE"),
                xform=np.count_nonzero(
                    np.isin(
                        self._rec_names,
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
                coord=np.count_nonzero(np.isin(self._rec_names, ["ATOM", "HETATM"])),
                ter=np.count_nonzero(self._rec_names == "TER"),
                connect=np.count_nonzero(self._rec_names == "CONECT"),
                seqres=np.count_nonzero(self._rec_names == "SEQRES"),
            )
        # Look for existing MASTER records in file
        idx_master_record = np.argwhere(self._rec_names == "MASTER").rehape(-1)
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
        master_record_chars = self._rec_chars[idx_master_record[-1]]
        master_record = pdbstruct.BookkeepingSection(
            remark=master_record_chars[10:15].view(dtype=(str, 5)).astype(int),
            het=master_record_chars[20:25].view(dtype=(str, 5)).astype(int),
            helix=master_record_chars[25:30].view(dtype=(str, 5)).astype(int),
            sheet=master_record_chars[30:35].view(dtype=(str, 5)).astype(int),
            site=master_record_chars[40:45].view(dtype=(str, 5)).astype(int),
            xform=master_record_chars[45:50].view(dtype=(str, 5)).astype(int),
            coord=master_record_chars[50:55].view(dtype=(str, 5)).astype(int),
            ter=master_record_chars[55:60].view(dtype=(str, 5)).astype(int),
            connect=master_record_chars[60:65].view(dtype=(str, 5)).astype(int),
            seqres=master_record_chars[65:70].view(dtype=(str, 5)).astype(int),
        )
        if master_record != bookkeeping_section:
            self._raise_or_warn(
                "MASTER record values do not match with number of records in file. "
                f"MASTER record values are {master_record}, "
                f"but actual record counts are {bookkeeping_section}.",
                raise_level=2,
            )

        if idx_master_record != self._rec_names.size - 2:
            self._raise_or_warn(
                f"MASTER record found in wrong position (line {idx_master_record[0] + 1}). "
                f"It must be the second-last record in the file, "
                f"i.e. at line {self._rec_names.size - 2}.",
                raise_level=3,
            )
        empty_cols = np.concatenate((np.arange(7, 11), np.arange(71, 81))) - 1
        if np.any(self._rec_chars[idx_master_record[0], empty_cols] != " "):
            self._raise_or_warn(
                "MASTER record line contains characters in columns that must be empty (i.e. "
                f"columns 7–10 and 71–80). "
                f"Found: {self._rec_chars[idx_master_record[0], empty_cols]}",
                raise_level=3,
            )
        return bookkeeping_section

    @property
    def record_end(self) -> NoReturn:
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
        idx_end_record = np.argwhere(self._rec_names == "END").rehape(-1)
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
        if idx_end_record != self._rec_names.size - 1:
            self._raise_or_warn(
                f"END record found at the middle of the file (line {idx_end_record[0] + 1}). "
                "It must be the final record in the file.",
                raise_level=3,
            )
        if np.any(self._rec_chars[-1, 6:] != " "):
            self._raise_or_warn(
                "END record line contains extra characters. It must contain only the record name "
                "'END' at beginning of the line. "
                f"Found: {self._rec_chars[-1, 6:].view(dtype=(str, self._rec_chars.shape[1] - 6)).reshape(-1)}",
                raise_level=3,
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
