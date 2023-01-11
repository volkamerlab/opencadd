from typing import Tuple
from abc import ABC, abstractmethod
from typing import Union, Tuple, Sequence
import datetime

import numpy as np
import pandas as pd


class RecordFieldDataType(ABC):
    @staticmethod
    @abstractmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        ...

    @staticmethod
    @abstractmethod
    def to_pdb(fields):
        ...


class AChar(RecordFieldDataType):
    """
    An alphabetic character (A-Z, a-z).
    """


class Atom(RecordFieldDataType):
    """
    Atom name.
    """


class Character(RecordFieldDataType):
    """
    Any non-control character in the ASCII character set or a space.
    """


class Continuation(RecordFieldDataType):
    """
    A two-character field that is either blank (for the first record of a set) or contains a
    two-digit number right-justified and blank-filled which counts continuation records starting
    with 2. The continuation number must be followed by a blank.
    """
    @staticmethod
    def from_pdb(fields: np.ndarray):
        fields[fields == "  "] = "1"
        return fields.astype(int)

    @staticmethod
    def to_pdb(num_lines: int):
        return ["  "] + [f"{line_num:>2}" for line_num in range(2, num_lines+1)]


class Date(RecordFieldDataType):
    """
    A 9 character string in the form DD-MMM-YY where DD is the
    day of the month, zero-filled on the left (e.g., 04); MMM is
    the common English 3-letter abbreviation of the month; and
    YY is the last two digits of the year. This must represent
    a valid date.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]) -> Union[datetime.date, Tuple[datetime.date]]:
        if isinstance(fields, str):
            return datetime.datetime.strptime(fields, "%d-%b-%y").date()
        return tuple(datetime.datetime.strptime(field, "%d-%b-%y").date() for field in fields)

    @staticmethod
    def to_pdb(dates: Union[datetime.date, Sequence[datetime.date]]) -> Union[str, Tuple[str]]:
        if isinstance(dates, datetime.date):
            return dates.strftime("%d-%b-%y").upper()
        return tuple(date.strftime("%d-%b-%y").upper() for date in dates)


class IDcode(RecordFieldDataType):
    """
    A PDB identification code which consists of 4 characters,
    the first of which is a digit in the range 0 - 9; the
    remaining 3 are alphanumeric, and letters are upper case
    only. Entries with a 0 as the first character do not
    contain coordinate data.
    """
    @staticmethod
    def process(fields: np.ndarray):
        return fields


class Integer(RecordFieldDataType):
    """
    Right-justified blank-filled integer value.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if isinstance(fields, str):
            return int(fields)
        return np.asarray(fields).astype(int)


class Token(RecordFieldDataType):
    """
    A sequence of non-space characters followed by a colon and a space.
    """


class List(RecordFieldDataType):
    """
    A String that is composed of text separated with commas.
    """


class LString(RecordFieldDataType):
    """
    A literal string of characters. All spacing is significant and must be preserved.
    """


class Real(RecordFieldDataType):
    """
    Real (floating point) number.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if isinstance(fields, str):
            return float(fields)
        return np.asarray(fields).astype(np.double)


class RecordName(RecordFieldDataType):
    """
    The name of the record: 6 characters, left-justified and blank-filled.
    """


class ResidueName(RecordFieldDataType):
    """
    One of the standard amino acid or nucleic acids, as listed
    below, or the non-standard group designation as defined in
    the HET dictionary. Field is right-justified.
    """


class SList(RecordFieldDataType):
    """
    A String that is composed of text separated with semicolons.
    """


class Specification(RecordFieldDataType):
    """
    A String composed of a token and its associated value separated by a colon.
    """


class SpecificationList(RecordFieldDataType):
    """
    A sequence of Specifications, separated by semi-colons.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if not isinstance(fields, str):
            fields = String.from_pdb(fields)
        specification_list = fields.split(";")
        token_value_pairs = np.char.strip([spec.split(":") for spec in specification_list])
        return token_value_pairs


class String(RecordFieldDataType):
    """
    A sequence of characters. These characters may have
    arbitrary spacing. To interpret a String, concatenate the contents of all continued fields
    together, collapse all sequences of multiple blanks to a single blank, and remove any
    leading and trailing blanks.
    """
    @staticmethod
    def from_pdb(fields: Sequence):
        return " ".join("".join(fields).split())


class SymOP(RecordFieldDataType):
    """
    An integer field of from 4 to 6 digits, right-justified, of
    the form nnnMMM where nnn is the symmetry operator number and
    MMM is the translation vector.
    """


def header_classification(fields: np.ndarray) -> Tuple[Tuple[str]]:
    """
    Parse the classification field of a HEADER record.

    Parameters
    ----------
    fields : ndarray, shape: (num_fields, num_chars), dtype: str
        Array of field values.

    Returns
    -------
    Tuple[Tuple[str]]
        Each sub-tuple corresponds to one molecule in the entry, with each element
        describing one classification/function of that molecule.
    """
    # HEADER is a one-time/single-line record. Take the first occurrence:
    field = fields[0]
    # The classification string is left-justified, and can describe dual functions of
    # molecules (when applicable) separated by a comma “,”. Entries with multiple
    # molecules in a complex will list the classifications of each macromolecule
    # separated by slash “/”.
    # First, split by '/' and then by ',' to get a tuple of tuples
    class_per_entity = tuple(field.split("/"))
    class_per_entity_and_function = tuple(
        tuple(entity_function.strip() for entity_function in entity_class.split(","))
        for entity_class in class_per_entity
    )
    return class_per_entity_and_function


def compound_compound(token_value_pairs: np.ndarray):
    df = pd.DataFrame(
        columns=[
            "MOL_ID"
            "MOLECULE",
            "CHAIN",
            "FRAGMENT",
            "SYNONYM",
            "EC",
            "ENGINEERED",
            "MUTATION",
            "OTHER_DETAILS"
        ]
    )
    df.index.name = "MOL_ID"
    ind_mols = np.argwhere(token_value_pairs[:, 0] == "MOL_ID").reshape(-1)
    for start, stop in np.lib.stride_tricks.sliding_window_view((*ind_mols, None), 2):
        mol_id = token_value_pairs[start, 1]
        data = dict(token_value_pairs[start + 1:stop])
        if "CHAIN" in data:
            data["CHAIN"] = np.char.strip(data["CHAIN"].split(","))
        if "SYNONYM" in data:
            data["SYNONYM"] = np.char.strip(data["SYNONYM"].split(","))
        if "EC" in data:
            data["EC"] = np.char.strip(data["EC".split(",")])
        if "ENGINEERED" in data:
            if not data["ENGINEERED"] in ("YES", "NO"):
                raise
            data["ENGINEERED"] = data["ENGINEERED"] == "YES"
        if "MUTATION" in data:
            if not data["MUTATION"] in ("YES", "NO"):
                raise
            data["MUTATION"] = data["MUTATION"] == "YES"
        df.loc[mol_id] = data
    return df


def atom_charge(fields: np.ndarray):
    """
    Parse the charge column and turn values into float.
    Values, if present, are composed of a digit followed by a '+' or '-' sign.
    """
    # Swap first and second columns to put the +/- sign in front of digit
    fields_reversed = np.array([charge[::-1] for charge in fields], dtype=(str, 3))
    empty_vals = fields_reversed == "  "
    if np.any(empty_vals):
        fields_reversed[empty_vals] = "nan"
        return fields_reversed.astype(np.single)
    return fields_reversed.astype(np.byte)