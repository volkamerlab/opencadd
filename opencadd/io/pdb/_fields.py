from typing import Tuple
from abc import ABC, abstractmethod
from typing import Union, Tuple, Sequence, Optional, Callable, Any, Type
import datetime
import re
import functools

import numpy as np
import numpy.typing as npt
import pandas as pd

from opencadd.io import _parsing


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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


class Atom(RecordFieldDataType):
    """
    Atom name.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


class Character(RecordFieldDataType):
    """
    Any non-control character in the ASCII character set or a space.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


class Continuation(RecordFieldDataType):
    """
    A two-character field that is either blank (for the first record of a set) or contains a
    two-digit number right-justified and blank-filled which counts continuation records starting
    with 2. The continuation number must be followed by a blank.
    """
    @staticmethod
    def from_pdb(fields: np.ndarray):
        fields = np.char.strip(fields)
        # idx_first_row = np.argwhere(fields == "")[0, 0]
        # fields[idx_first_row] = "1"
        fields[fields == ""] = "1"
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

    @staticmethod
    def verify(values):

        def verify_single(value):
            if not isinstance(value, datetime.date):
                raise TypeError(
                    "`deposition_date` must be of type `datetime.date`, "
                    f"but the type of input ({value}) was {type(value)}."
                )
            if value >= datetime.date.today():
                raise ValueError(
                    "`deposition_date` must be a date in the past, "
                    f"but the input date ({value}) is in the future; today is {datetime.date.today()}."
                )
            return

        if isinstance(values, datetime.date):
            verify_single(values)
        else:
            for value in values:
                verify_single(value)
        return


class Element(RecordFieldDataType):
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.capitalize(np.char.strip(fields))


class IDcode(RecordFieldDataType):
    """
    A PDB identification code which consists of 4 characters,
    the first of which is a digit in the range 0 - 9; the
    remaining 3 are alphanumeric, and letters are upper case
    only. Entries with a 0 as the first character do not
    contain coordinate data.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)

    @staticmethod
    def verify(fields: Union[str, Sequence[str]]):

        def verify_pdb_id(pdb_id: str):
            if len(pdb_id) != 4:
                raise ValueError(f"PDB ID must have 4 characters, but input {pdb_id} had {len(pdb_id)}.")
            if not pdb_id[0].isnumeric() or pdb_id[0] == "0":
                raise ValueError("First character of `pdb_id` must be a non-zero digit.")
            return

        if isinstance(fields, str):
            verify_pdb_id(fields)
        else:
            for pdb_id in fields:
                verify_pdb_id(pdb_id)
        return


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
        mask_empty_fields = fields == ""
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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return fields


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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return fields


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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


class ResidueName(RecordFieldDataType):
    """
    One of the standard amino acid or nucleic acids, as listed
    below, or the non-standard group designation as defined in
    the HET dictionary. Field is right-justified.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


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
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        return np.char.strip(fields)


class SpecificationList(RecordFieldDataType):
    """
    A sequence of Specifications, separated by semi-colons.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if not isinstance(fields, str):
            fields = String.from_pdb(fields)
        specification_list = fields.split(";")
        # TODO: some entries (e.g. 1Y34, 1PZY, 5SYD) contain errors in that they don't escape other uses of
        #  semicolon. For example in COMPND record of 1Y34:
        #  'COMPND   4 SYNONYM: SUBTILISIN NOVO; SUBTILISIN DFE; ALKALINE PROTEASE;         '
        #  This causes extra splits that should not be there. We should find a way to not split these extra
        #  non-escaped delimeters.
        token_value_pairs = np.char.strip([spec.split(":", maxsplit=1) for spec in specification_list])
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
        concatenated = " ".join("".join(fields).split())
        # When one line ends with a hyphen, there should be no space between that and the next line.
        return re.sub(r"- ", "-", concatenated)


class SymOP(RecordFieldDataType):
    """
    An integer field of from 4 to 6 digits, right-justified, of
    the form nnnMMM where nnn is the symmetry operator number and
    MMM is the translation vector.
    """
    @staticmethod
    def from_pdb(fields: Union[str, Sequence[str]]):
        if isinstance(fields, str):
            return np.array(list(fields), dtype=np.ubyte)
        return list(np.asarray(fields).astype(np.ubyte))
        # char_view = np.asarray(fields, order="C").view(dtype=(str, 1)).reshape(-1, 6)
        # #char_view = fields.view(dtype=(str, 1)).reshape(-1, 6)
        # sym_op_num = char_view[:, :3].view(dtype=(str, 3)).astype(int).reshape(-1)
        # transl_vec = char_view[:, 3:].astype(int)
        # print(sym_op_num, transl_vec)
        # return sym_op_num, transl_vec


class Column:
    def __init__(
            self,
            intervals: Union[Tuple[int, int], Sequence[Tuple[int, int]]],
            field_dtype: Type[RecordFieldDataType],
            strip: bool = True,
            cast: bool = True,
            only_first: bool = False,
            only_non_empty: bool = False,
    ):
        self._intervals = np.asarray(intervals)
        if not np.issubdtype(self._intervals.dtype, np.integer):
            raise ValueError(
                "Parameter `intervals` expects an array of integer types, "
                f"but input argument elements had type {self._intervals.dtype}. Input was:\n{intervals}"
            )
        if self._intervals.ndim == 1:
            if self._intervals.size != 2:
                raise ValueError(
                    "Parameter `interval` expects a 1D array of size 2, "
                    f"but input argument had size {self._intervals.size}. Input was:\n{intervals}."
                )
            extract = functools.partial(
                _parsing.extract_columns_by_interval,
                interval=self._intervals,
                strip=strip
            )
            self._indices = np.arange(self._intervals[0], self._intervals[1])
        elif self._intervals.ndim == 2:
            if self._intervals.shape[1] != 2:
                raise ValueError(
                    "Parameter `intervals` expects a 2D array of shape (n, 2), "
                    f"but input argument had shape {self._intervals.shape}. Input was:\n{intervals}."
                )
            col_lengths = self._intervals[:, 1] - self._intervals[:, 0]
            if np.any(col_lengths != col_lengths[0]):
                raise ValueError(
                    "Parameter `intervals` expects a sequence of intervals with equal lengths, "
                    f"but input argument had intervals of lengths {col_lengths}. Input was:\n{intervals}"
                )
            self._indices = np.concatenate(
                [np.arange(interval[0], interval[1]) for interval in self._intervals]
            )
            extract = functools.partial(
                _parsing.extract_columns_by_index,
                indices=self._indices,
                column_len=col_lengths[0],
                strip=strip
            )
        else:
            raise ValueError(
                f"Parameter `intervals` expects either a 1D or 2D array, "
                f"but input argument had {self._intervals.ndim} dimensions. Input was:\n{intervals}."
            )
        extract1 = extract if not cast else (lambda char_table: self.cast_to_dtype(extract(char_table=char_table)))
        if only_first:
            self.extract = lambda char_table: extract1(char_table)[0]
        elif only_non_empty:
            self.extract = lambda char_table: self._only_non_empty(extract1(char_table))
        else:
            self.extract = extract1
        self._dtype = field_dtype
        return

    @property
    def indices(self) -> np.ndarray:
        """
        Indices (i.e. all character positions in a line) of the column.
        """
        return self._indices

    def cast_to_dtype(self, fields: np.ndarray) -> np.ndarray:
        return self._dtype.from_pdb(fields)

    @staticmethod
    def _only_non_empty(fields):
        return fields[fields != ""]