from typing import Tuple
from abc import ABC, abstractmethod
from typing import Union, Tuple, Sequence, Optional
import datetime

import numpy as np
import numpy.typing as npt
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
    def from_pdb(fields: Union[str, Sequence[str]]) -> Union[datetime.date, np.ndarray]:
        def str_to_date(date_str) -> Optional[datetime.date]:
            try:
                return datetime.datetime.strptime(date_str, "%d-%b-%y").date()
            except ValueError as e:
                if date_str.isspace():
                    return
                raise
        if isinstance(fields, str):
            return str_to_date(fields)
        return np.array([str_to_date(field) for field in fields])

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
    def from_pdb(fields: Union[str, Sequence[str]], dtype: npt.DTypeLike = np.int_):
        dtype = np.dtype(dtype)
        if isinstance(fields, str):
            if fields.isspace():
                return np.iinfo(dtype).min
            return dtype.type(fields)
        fields = np.asarray(fields)
        mask_empty_fields = np.char.isspace(fields)
        if not np.any(mask_empty_fields):
            return fields.astype(dtype)
        ints = np.empty(shape=fields.size, dtype=dtype)
        ints[mask_empty_fields] = np.iinfo(dtype).min
        ints[~mask_empty_fields] = fields[~mask_empty_fields].astype(dtype)
        return ints


class Token(RecordFieldDataType):
    """
    A sequence of non-space characters followed by a colon and a space.
    """


class List(RecordFieldDataType):
    """
    A String that is composed of text separated with commas.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if not isinstance(fields, str):
            fields = String.from_pdb(fields)
        return np.char.strip(fields.split(","))



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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if not isinstance(fields, str):
            fields = String.from_pdb(fields)
        return np.char.strip(fields.split(";"))


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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if isinstance(fields, str):
            return int(fields[:3]), (int(fields[3]), int(fields[4]), int(fields[5]))
        fields = np.asarray(fields)

        char_view = np.frombuffer(
            fields.tobytes(), dtype=(str, 1)
        ).reshape(-1, 6)

        #char_view = fields.view(dtype=(str, 1)).reshape(-1, 6)
        sym_op_num = char_view[:, :3].view(dtype=(str, 3)).astype(int).reshape(-1)
        transl_vec = char_view[:, 3:].astype(int)
        return sym_op_num, transl_vec



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