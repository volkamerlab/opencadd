def validate_master(self) -> sections.Bookkeeping:
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
                                     records.MASTER.columns)
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


def validate_section_coordinate(self):
    if ind_ter.size == 0:
        if ind_atom.size != 0:
            self._raise_or_warn(
                "ATOM records found without TER record.",
                raise_level=1,
            )

    indices_record_model = self._idx__record_lines[records.MODEL]  # line indices of MODEL records
    count_model = indices_record_model.size  # count of MODEL records
    if count_model > 0:
        indices_record_endmdl = self._idx__record_lines[records.ENDMDL]
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


def validate_dbref(self):
    if not (
            self._has_record[records.DBREF]
            or self._has_record[records.DBREF1]
            or self._has_record[records.DBREF2]
    ):
        if self._has_record[records.ATOM] or self._has_record[records.TER]:
            self._raise_or_warn(
                "DBREF and DBREF1/DBREF2 records not found, while ATOM records exist. "
                "At least one must exist when ATOM records exist.",
                raise_level=3
            )
    if self._has_record[records.DBREF1]:
        if self._idx__record_lines[records.DBREF1].size != self._idx__record_lines[records.DBREF2].size:
            self._raise_or_warn(
                "Number of DBREF1 and DBREF2 records do not match. DBREF1/DBREF2 records are "
                "two-line records that always have to be present together.",
                raise_level=2
            )
        self.validate_empty_columns(records.DBREF1)
        self.validate_empty_columns(records.DBREF2)
    self._check_routine(records.DBREF)
    return


def validate_end(self) -> NoReturn:
    """
    Validate the END record in the PDB file.

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
    if not self.has_record(records.END):
        self._raise_or_warn(
            "END record not found. It is a mandatory record and must appear once at "
            "the end of the file.",
            raise_level=3,
        )
    idx_record_end = self.indices_record_lines(records.END)
    if self.count_record(records.END) > 1:
        self._raise_or_warn(
            f"Multiple END records found ({idx_record_end.size} at lines "
            f"{idx_record_end + 1}). It is a one-time/single-line record and must appear "
            "only once in the file.",
            raise_level=3,
        )
    if idx_record_end != self.count_lines - 1:
        self._raise_or_warn(
            f"END record found at the middle of the file (line {idx_record_end[0] + 1}). "
            "It must be the final record in the file.",
            raise_level=3,
        )
    lines_record_end = self.record_lines(records.END, drop_name=True)
    mask_non_empty = lines_record_end != " "
    if np.any(mask_non_empty):
        self._raise_or_warn(
            "END record line contains extra characters. It must contain only the record name "
            "'END' at beginning of the line. "
            f"Found characters {lines_record_end[mask_non_empty]} at positions "
            f"{np.argwhere(mask_non_empty)}.",
            raise_level=3,
        )
    return

def _parse_record_in_model(self, record: records.Record) -> Optional[pd.DataFrame]:
    if not self.has_record(record):
        return
    fields = self.extract_record_fields(record)
    if not self.has_record(records.MODEL):
        model_nr = np.full(shape=self.count_record(record), fill_value=1, dtype=np.ubyte)
    else:
        count_record_per_model = [
            np.count_nonzero(
                np.logical_and(
                    self.indices_record_lines(record) > model_start_idx,
                    self.indices_record_lines(record) < model_end_idx,
                )
            ) for model_start_idx, model_end_idx in zip(
                self.indices_record_lines(records.MODEL), self.indices_record_lines(records.ENDMDL)
            )
        ]
        model_nr = np.repeat(
            np.arange(1, self.count_record(records.MODEL) + 1), repeats=count_record_per_model
        )
    return pd.DataFrame({"model_nr": model_nr, **fields})

def _check_routine(self, record: records.Record):
    if self.has_record(record):
        self._check_records_are_consecutive(record)
        self.validate_empty_columns(record)

def _check_records_are_consecutive(self, record: records.Record):
    if self._has_record[record] and np.any(np.diff(self._idx__record_lines[record]) != 1):
        self._raise_or_warn(
            f"{record.name} records are not in consecutive lines.",
            raise_level=3
        )
    return


if self.strictness > 0:
    self._validate_one_time_single_line_record(records.HEADER)
    self.validate_empty_columns(records.HEADER)
    # Check position in file; HEADER must always be in first line:
    indices_lines = self.indices_record_lines(records.HEADER)
    if indices_lines[0] != 0:
        self._raise_or_warn(
            f"HEADER record found at line {indices_lines[0]}. "
            f"It must always be the first line of a PDB file.",
            raise_level=3
        )




if self.strictness:
    if header is not None:
        if obslte is not None:
            if header["pdb_id"] != obslte["pdb_id"]:
                self._raise_or_warn(
                    f"Mismatched data: PDB ID in the HEADER record ({header['pdb_id']}) "
                    f"does not match that of the OBSLTE record ({obslte['pdb_id']}).",
                    raise_level=2
                )
        if caveat is not None:
            if header["pdb_id"] != caveat["pdb_id"]:
                self._raise_or_warn(
                    f"Mismatched data: PDB ID in the HEADER record ({header['pdb_id']}) "
                    f"does not match that of the CAVEAT record ({caveat['pdb_id']}).",
                    raise_level=2
                )