from typing import NamedTuple, List, Tuple
import warnings
import datetime

import numpy as np
import pandas as pd

from opencadd.io.parser.parser import extract_column_from_string_array

RECORDS = [
    "HEADER",
    "OBSLTE",
    "TITLE",
    "SPLIT",
    "CAVEAT",
    "COMPND",
    "SOURCE",
    "KEYWDS",
    "EXPDTA",
    "NUMMDL",
    "MDLTYP",
    "AUTHOR",
    "REVDAT",
    "SPRSDE",
    "JRNL",
    "REMARK",
    "DBREF",
    "SEQADV",
    "SEQRES",
    "MODRES",
    "HET",
    "HETNAM",
    "HETSYN",
    "FORMUL",
    "HELIX",
    "SHEET",
    "SSBOND",
    "LINK",
    "CISPEP",
    "SITE",
    "CRYST1",
    "ORIGX1",
    "ORIGX2",
    "ORIGX3",
    "SCALE1",
    "SCALE2",
    "SCALE3",
    "MTRIX1",
    "MTRIX2",
    "MTRIX3",
    "MODEL",
    "ATOM",
    "ANISOU",
    "TER",
    "HETATM",
    "ENDMDL",
    "CONECT",
    "MASTER",
    "END",
]


class HeaderRecord(NamedTuple):
    pdb_id: str = None
    classification: Tuple[Tuple[str]] = None
    deposition_date: datetime.date = None
    invalid_lines: Tuple = None


class TitleSection(NamedTuple):
    header: HeaderRecord

class Model(NamedTuple):
    atom: pd.DataFrame = None
    anisou: pd.DataFrame = None
    hetatm: pd.DataFrame = None
    ter: Tuple = None

class CoordinateSection(NamedTuple):
    models: Tuple[Model]

class PDBFile(NamedTuple):
    title: TitleSection
    coordinate: CoordinateSection


def parse_pdb(filepath, strict: bool = True):
    with open(filepath, "r") as f:
        pdb_lines = np.array(f.read().splitlines())

    # Check if all lines are 80 characters and raise a warning if not:
    line_widths = np.char.str_len(pdb_lines)
    lines_not80 = line_widths != 80
    if np.any(lines_not80):
        faulty_lines = [
            f"{line_num} ({line_width})" for line_num, line_width in zip(
                np.argwhere(lines_not80).reshape(-1),
                line_widths[lines_not80]
            )
        ]
        warnings.warn(
            "Following lines are not 80 characters long (lengths given in brackets): "
            f"{', '.join(faulty_lines)}"
        )

    # Extract all record types
    record_types = np.char.strip(extract_column_from_string_array(pdb_lines, (0, 6)))

    # Parse Title section
    title = read_title_section(record_types, pdb_lines)
    coordinate = read_coordinate_section(record_types, pdb_lines)

    pdbfile = PDBFile(
        title=title,
        coordinate=coordinate
    )
    return pdbfile


def read_title_section(record_types: np.ndarray, pdb_lines: np.ndarray) -> TitleSection:
    header = read_header_record(record_types, pdb_lines)

    title = TitleSection(
        header=header
    )
    return title


def read_coordinate_section(record_types: np.ndarray, pdb_lines: np.ndarray):
    filter_models = record_types == "MODEL"
    models_start_ind = np.argwhere(filter_models).reshape(-1)
    if models_start_ind.size == 0:
        filter_model = np.isin(record_types, ["ATOM", "ANISOU", "TER", "HETATM"])
        model = read_model(record_types[filter_model], pdb_lines[filter_model])
        return CoordinateSection(models=(model,))
    # There are multiple models present:
    models_end_ind = np.argwhere(record_types == "ENDMDL").reshape(-1)
    if models_start_ind.size != models_end_ind.size:
        warnings.warn(
            "Number of MODEL and ENDMDL records do not match: "
            f"{models_start_ind.size} vs {models_end_ind.size}"
        )
    models = []
    for start_idx, end_idx in zip(models_start_ind, models_end_ind):
        models.append(
            read_model(
                record_types=record_types[start_idx+1:end_idx],
                pdb_lines=pdb_lines[start_idx+1:end_idx]
            )
        )
    return CoordinateSection(models=tuple(models))


def read_model(record_types: np.ndarray, pdb_lines: np.ndarray) -> Model:
    has_atom = np.isin("ATOM", record_types)
    atom = read_atom_record(pdb_lines[record_types == "ATOM"]) if has_atom else None
    return Model(atom=atom)


def read_atom_record(pdb_lines: np.ndarray) -> pd.DataFrame:
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
    # Parse the charge column if present:
    #print(df.charge.values, df.charge.values.dtype)
    charges = df.charge.values.astype("<U2").view(dtype=(str, 1)).reshape(
        pdb_lines.size, -1
    )
    # Swap first and second columns to put the +/- sign in front of digit
    charges[:, [0, 1]] = charges[:, [1, 0]]
    # Combine again
    charges = np.frombuffer(charges.tobytes(), dtype=(str, 2))
    # Replace empty strings with 'nan'
    charges_cleaned = np.where(charges != "", charges, "nan")
    # Turn to float and asign to dataframe
    df.charge = charges_cleaned.astype(np.single)
    return df


def read_header_record(record_types: np.ndarray, pdb_lines: np.ndarray) -> HeaderRecord:
    """
    Parse the HEADER records of a PDB file.

    Parameters
    ----------
    record_types : ndarray, shape: (n, ), dtype: <U6
        1-dimensional array of 6-character strings specifying the record type of each line in
        the PDB file.
    pdb_lines : ndarray, shape: (n, ), dtype: str
        1-dimensional array of strings representing the full lines of the PDB file.

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
    if record_types[0] != "HEADER":
        warnings.warn(f"First line is not a HEADER record, instead: {pdb_lines[0]}")
    # HEADER is a one-time/single-line record:
    is_header = record_types == "HEADER"
    header_lines = pdb_lines[is_header]
    header_is_one_time = header_lines.size == 1
    if not header_is_one_time:
        warnings.warn(
            f"There must be only one HEADER record in the file, but found {header_lines.size} "
            f"at lines {', '.join(np.argwhere(is_header).reshape(-1))}: {header_lines}"
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
        return HeaderRecord(
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
    header = HeaderRecord(
        pdb_id=pdb_id,
        classification=class_per_entity_and_function,
        deposition_date=datetime.datetime.strptime(dep_date, "%d-%b-%y").date()
    )
    return header


def write_header_record(
        header: HeaderRecord,
        valid_only: bool = True,
        default_pdb_id: str = "XXXX",
        default_classification: str = "N/A",
        default_deposition_date: str = "??-???-??",
) -> str:
    pdb_id = header.pdb_id if header.pdb_id is not None else default_pdb_id
    classification = "/".join(
        ", ".join(entity_functions) for entity_functions in header.classification
    ) if header.classification is not None else default_classification
    deposition_date = header.deposition_date.strftime("%d-%b-%y").upper(
    ) if header.deposition_date is not None else default_deposition_date
    valid_header = f"HEADER{'':4}{classification:<40}{deposition_date}{'':3}{pdb_id}{'':14}"
    return valid_header if valid_only else valid_header + "\n" + "\n".join(header.invalid_lines)


