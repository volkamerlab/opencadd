"""
Base functionalities to parse an mmCIF datafile or dictionary.
"""

import re
import enum
from typing import Union, Iterator
import itertools
import polars as pl
import pandas as pd

from . import _filestruct

__author__ = "Armin Ariamajd"


TOKENIZER = re.compile(
    r"""(?xmi)  # x (cf. re.X) allows for writing the expression in multiple lines, with comments added; 
                # m (cf. re.M, re.MULTILINE) causes the pattern characters '^' and '$' to also match
                #  the beggining and end of each line, respectively 
                #  (in addition to matching the beggining and end of the whole string, respectively).
                # i (cf. re.I, re.IGNORECASE) performs case-insensitive matching.
    # The following creates different capturing groups (enumerated starting from 1), 
    #  each matching one token type. Notice the order of groups matters, 
    #  since the matching terminates with the first group match.
    ^;([\S\s]*?)(?:\r\n|\s)^;(?:(?=\s)|$)  # 1. Text field, i.e. a non-simple data value 
                                           #    bounded between two '\n;' characters.
    |(?:^|(?<=\s))\#(.*?)\r?$              # 2. Comment
    |(?:^|(?<=\s))(?:
      '(.*?)'                              # 3. Quoted data value
      |"(.*?)"                             # 4. Duble-quoted data value
      |_(\S+)                              # 5. Data name
      |loop_(\S*)                          # 6. Loop header
      |data_(\S*)                          # 7. Block code
      |save_(\S*)                          # 8. Frame code (or terminator)
      |stop_(\S*)                          # 9. STAR-reserved loop terminator
      |global_(\S*)                        # 10. STAR-reserved global block header
      |\$(\S+)                             # 11. STAR-reserved frame reference
      |\[(.+?)]                            # 12. STAR-reserved multi-line value delimeter
      |((?:[^'";_$\[\s]|(?<!^);)\S*)       # 13. Data value
      |(\S+)                               # 14. Bad token (anything else)
    )
    (?:(?=\s)|$)"""
)
"""
Regular expression to capture tokens in an mmCIF file.

This expression can be used on a single multi-line string 
representing the whole content of an mmCIF file.
Used in iterative mode, it will then tokenize the whole file
(tokens are separated by any whitespace character that is 
not encapsulated in a non-simple data value delimiter, as described in the CIF documentation),
and identify the type of each captured token.
"""


class Token(enum.Enum):
    """Types of Tokens in a CIF file.

    The values correspond to the index of capturing groups in `TOKENIZER` above.
    """
    VALUE_FIELD = 1
    COMMENT = 2
    VALUE_QUOTED = 3
    VALUE_DOUBLE_QUOTED = 4
    NAME = 5
    LOOP = 6
    DATA = 7
    SAVE = 8
    STOP = 9
    GLOBAL = 10
    FRAME_REF = 11
    BRACKETS = 12
    VALUE = 13
    BAD_TOKEN = 14
    SAVE_END = 15  # Will be added by the parser after parsing SAVE


class State(enum.Enum):
    """
    Possible states of the parser in a CIF file.

    Attributes
    ----------
    IN_FILE
        The initial state; the parser is in the file
        but has not yet encountered the first data block header.
        In this state, the only expected token (beside comments, which are always expected)
        is a data block header, i.e. `Token.DATA`.
    JUST_IN_DATA
        The parser has just encountered a data block header,
        and is now expecting either a save frame header (i.e. `Token.SAVE`),
        a loop directive (i.e. `Token.LOOP`), or a data name (i.e. `Token.NAME`).
    JUST_IN_SAVE
        The parser has just encountered a save frame header,
        and is now expecting either a loop directive, or a data name.
    JUST_IN_LOOP
        The parser has just encountered a loop directive,
        and is now expecting a data name.
    IN_NAME
        The parser was in a data block and has just encountered a data name.
        It is now expecting a data value, i.e. either
        `Token.VALUE`, `Token.VALUE_QUOTED` or `Token.VALUE_FIELD`.
    JUST_IN_SAVE_LOOP
        The parser was in a save frame and has just encountered a loop directive.
        It is now expecting a data name.
    IN_SAVE_NAME
        The parser was in a save frame and has just encountered a data name.
        It is now expecting a data value.
    IN_LOOP_NAME
        The parser was in a loop and has just encountered the first or n-th data name.
        It is now expecting either another data name, or the first data value.
    IN_DATA
        The parser was either in a save frame and has just encountered a save frame termination directive
        (i.e. `Token.SAVE_END`), or it was in a data block and has just finished parsing a data item
        (i.e. it encountered a data name, followed by a data value). It is now in a data block
        and expecting either another data block header, a save frame header, a loop directive,
        or a data name.
    IN_SAVE_LOOP_NAME
        The parser was in a loop inside a save frame and has just encountered the first or n-th data name.
        It is now expecting either another data name, or the first data value.
    IN_SAVE
        The parser has just encountered a data value (directly after a data name) within a save frame.
        It is now inside the save frame expecting either a save frame termination directive, a loop directive,
        or a data name.
    IN_LOOP_VALUE
        The parser has encountered the first or n-th data value,
        while being previously in state `IN_LOOP_NAME`. It is now expecting either another data value,
        a new data block header, a save frame header, another loop directive, or a data name.
    IN_SAVE_LOOP_VALUE
        The parser has encountered the first or n-th data value,
        while being previously in state `IN_SAVE_LOOP_NAME`. It is now expecting either another data value,
        a save frame termination directive, another loop directive, or a data name.
    """
    IN_FILE = enum.auto()
    JUST_IN_DATA = enum.auto()
    JUST_IN_SAVE = enum.auto()
    JUST_IN_LOOP = enum.auto()
    IN_NAME = enum.auto()
    JUST_IN_SAVE_LOOP = enum.auto()
    IN_SAVE_NAME = enum.auto()
    IN_LOOP_NAME = enum.auto()
    IN_DATA = enum.auto()
    IN_SAVE_LOOP_NAME = enum.auto()
    IN_SAVE = enum.auto()
    IN_LOOP_VALUE = enum.auto()
    IN_SAVE_LOOP_VALUE = enum.auto()


EXPECTED_TOKENS = {
        State.IN_FILE: (Token.DATA, ),
        State.JUST_IN_DATA: (Token.SAVE, Token.LOOP, Token.NAME),
        State.JUST_IN_SAVE: (Token.LOOP, Token.NAME),
        State.JUST_IN_LOOP: (Token.NAME, ),
        State.IN_NAME: (Token.VALUE, Token.VALUE_QUOTED, Token.VALUE_FIELD),
        State.JUST_IN_SAVE_LOOP: (Token.NAME, ),
        State.IN_SAVE_NAME: (Token.VALUE, Token.VALUE_QUOTED, Token.VALUE_FIELD),
        State.IN_LOOP_NAME: (Token.NAME, Token.VALUE, Token.VALUE_QUOTED, Token.VALUE_FIELD),
        State.IN_DATA: (Token.DATA, Token.SAVE, Token.LOOP, Token.NAME),
        State.IN_SAVE_LOOP_NAME: (Token.NAME, Token.VALUE, Token.VALUE_QUOTED, Token.VALUE_FIELD),
        State.IN_SAVE: (Token.SAVE_END, Token.LOOP, Token.NAME),
        State.IN_LOOP_VALUE: (
            Token.DATA,
            Token.SAVE,
            Token.LOOP,
            Token.NAME,
            Token.VALUE,
            Token.VALUE_QUOTED,
            Token.VALUE_FIELD
        ),
        State.IN_SAVE_LOOP_VALUE: (
            Token.SAVE_END,
            Token.LOOP,
            Token.NAME,
            Token.VALUE,
            Token.VALUE_QUOTED,
            Token.VALUE_FIELD
        ),
    }
"""A mapping from each possible state to the expected tokens in that state.

This is only used for generating error messages when unexpected tokens are encountered.
Otherwise, a more verbose mapping is defined in `MMCIFParser.__init__`, used by the parser
to validate and update its state. 
"""


class CIFParser:
    """
    Base class for mmCIF parsers.

    """

    def __init__(self, content: str, strictness):

        do_nothing = lambda: None

        self._state_mapper = {
            (State.IN_FILE, Token.DATA): (do_nothing, State.JUST_IN_DATA),

            (State.JUST_IN_DATA, Token.SAVE): (do_nothing, State.JUST_IN_SAVE),
            (State.JUST_IN_DATA, Token.LOOP): (do_nothing, State.JUST_IN_LOOP),
            (State.JUST_IN_DATA, Token.NAME): (do_nothing, State.IN_NAME),

            (State.JUST_IN_SAVE, Token.LOOP): (do_nothing, State.JUST_IN_SAVE_LOOP),
            (State.JUST_IN_SAVE, Token.NAME): (do_nothing, State.IN_SAVE_NAME),

            (State.JUST_IN_LOOP, Token.NAME): (self._initialize_loop, State.IN_LOOP_NAME),

            (State.IN_NAME, Token.VALUE_FIELD): (self._add_data_item, State.IN_DATA),
            (State.IN_NAME, Token.VALUE_QUOTED): (self._add_data_item, State.IN_DATA),
            (State.IN_NAME, Token.VALUE): (self._add_data_item, State.IN_DATA),

            (State.JUST_IN_SAVE_LOOP, Token.NAME): (self._initialize_loop, State.IN_SAVE_LOOP_NAME),

            (State.IN_SAVE_NAME, Token.VALUE_FIELD): (self._add_data_item, State.IN_SAVE),
            (State.IN_SAVE_NAME, Token.VALUE_QUOTED): (self._add_data_item, State.IN_SAVE),
            (State.IN_SAVE_NAME, Token.VALUE): (self._add_data_item, State.IN_SAVE),

            (State.IN_LOOP_NAME, Token.NAME): (self._add_loop_keyword, State.IN_LOOP_NAME),
            (State.IN_LOOP_NAME, Token.VALUE_FIELD): (self._register_and_fill_loop, State.IN_LOOP_VALUE),
            (State.IN_LOOP_NAME, Token.VALUE_QUOTED): (self._register_and_fill_loop, State.IN_LOOP_VALUE),
            (State.IN_LOOP_NAME, Token.VALUE): (self._register_and_fill_loop, State.IN_LOOP_VALUE),

            (State.IN_DATA, Token.DATA): (do_nothing, State.JUST_IN_DATA),
            (State.IN_DATA, Token.SAVE): (do_nothing, State.JUST_IN_SAVE),
            (State.IN_DATA, Token.LOOP): (do_nothing, State.JUST_IN_LOOP),
            (State.IN_DATA, Token.NAME): (do_nothing, State.IN_NAME),

            (State.IN_SAVE_LOOP_NAME, Token.NAME): (self._add_loop_keyword, State.IN_SAVE_LOOP_NAME),
            (State.IN_SAVE_LOOP_NAME, Token.VALUE_FIELD): (self._register_and_fill_loop, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_NAME, Token.VALUE_QUOTED): (self._register_and_fill_loop, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_NAME, Token.VALUE): (self._register_and_fill_loop, State.IN_SAVE_LOOP_VALUE),

            (State.IN_SAVE, Token.SAVE_END): (do_nothing, State.IN_DATA),
            (State.IN_SAVE, Token.LOOP): (do_nothing, State.JUST_IN_SAVE_LOOP),
            (State.IN_SAVE, Token.NAME): (do_nothing, State.IN_SAVE_NAME),

            (State.IN_LOOP_VALUE, Token.DATA): (self._finalize_loop, State.JUST_IN_DATA),
            (State.IN_LOOP_VALUE, Token.SAVE): (self._finalize_loop, State.JUST_IN_SAVE),
            (State.IN_LOOP_VALUE, Token.LOOP): (self._finalize_loop, State.JUST_IN_LOOP),
            (State.IN_LOOP_VALUE, Token.NAME): (self._finalize_loop, State.IN_NAME),
            (State.IN_LOOP_VALUE, Token.VALUE_FIELD): (self._fill_loop_value, State.IN_LOOP_VALUE),
            (State.IN_LOOP_VALUE, Token.VALUE_QUOTED): (self._fill_loop_value, State.IN_LOOP_VALUE),
            (State.IN_LOOP_VALUE, Token.VALUE): (self._fill_loop_value, State.IN_LOOP_VALUE),

            (State.IN_SAVE_LOOP_VALUE, Token.SAVE_END): (self._finalize_loop, State.IN_DATA),
            (State.IN_SAVE_LOOP_VALUE, Token.LOOP): (self._finalize_loop, State.JUST_IN_SAVE_LOOP),
            (State.IN_SAVE_LOOP_VALUE, Token.NAME): (self._finalize_loop, State.IN_SAVE_NAME),
            (State.IN_SAVE_LOOP_VALUE, Token.VALUE_FIELD): (self._fill_loop_value, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_VALUE, Token.VALUE_QUOTED): (self._fill_loop_value, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_VALUE, Token.VALUE): (self._fill_loop_value, State.IN_SAVE_LOOP_VALUE),
        }

        self._token_processor = {
            Token.DATA: self._process_data,
            Token.SAVE: self._process_save,
            Token.LOOP: self._process_loop,
            Token.NAME: self._process_name,
            Token.VALUE: self._process_value,
            Token.VALUE_QUOTED: self._process_value_quoted,
            Token.VALUE_DOUBLE_QUOTED: self._process_value_double_quoted,
            Token.VALUE_FIELD: self._process_value_field,
        }

        self._tokenizer: Iterator[re.Match] = None

        self._curr_state: State = State.IN_FILE
        self._curr_match: re.Match = None
        self._curr_token_type: Token = None
        self._curr_token_value: str = None
        self._curr_block_code: str = None
        self._curr_frame_code_category: str = None
        self._curr_frame_code_keyword: str = None
        self._curr_data_name_category: str = None
        self._curr_data_name_keyword: str = None
        self._curr_data_value: str = None


        self._cif_dict: dict = dict()

        self._block_codes: list = []
        self._frame_code_categories: list = []
        self._frame_code_keywords: list = []
        self._data_name_categories: list = []
        self._data_name_keywords: list = []
        self._data_values: list = []
        self._value_dimensions: list = []

        self._cif_str = content
        self._tokenizer = TOKENIZER.finditer(content)
        return

    @property
    def _curr_data_block_dict(self) -> dict:
        return self._cif_dict.setdefault(self._curr_block_code, dict())

    @property
    def _curr_save_frame_category_dict(self) -> dict:
        return self._curr_data_block_dict.setdefault(self._curr_frame_code_category, dict())

    @property
    def _curr_save_frame_keyword_dict(self) -> dict:
        return self._curr_save_frame_category_dict.setdefault(self._curr_frame_code_keyword, dict())

    @property
    def _curr_data_category_dict(self) -> dict:
        return self._curr_save_frame_keyword_dict.setdefault(self._curr_data_name_category, dict())

    def parse(self):
        self._curr_state = State.IN_FILE
        for match in self._tokenizer:
            self._curr_match = match
            self._curr_token_type = Token(match.lastindex)
            self._curr_token_value = match.group(match.lastindex)
            self._token_processor.get(self._curr_token_type, lambda: None)()
            update_func, new_state = self._state_mapper.get(
                (self._curr_state, self._curr_token_type),
                (self._unexpected, self._curr_state)
            )
            update_func()
            self._curr_state = new_state
        return self._finalize()

    def _unexpected(self):
        if self._curr_token_type == Token.COMMENT:
            return
        if self._curr_token_type == Token.BAD_TOKEN:
            raise CIFParsingError(
                f"Bad token: got {self._curr_token_value} at position {self._curr_match.span()}"
            )
        if self._curr_token_type in [Token.STOP, Token.GLOBAL, Token.FRAME_REF, Token.BRACKETS]:
            raise
        raise CIFParsingError(
            f"Token out of order: parser is in state {self._curr_state} and expects a token from {EXPECTED_TOKENS[self._curr_state]} but received a {self._curr_token_type}"
        )

    def _process_data(self):
        if self._curr_token_value == "":
            raise
        self._curr_block_code = self._curr_token_value
        return

    def _process_save(self):
        if self._curr_token_value == "":
            self._curr_token_type = Token.SAVE_END
            self._curr_frame_code_category = None
            self._curr_frame_code_keyword = None
            return
        frame_code_components = self._curr_token_value.split(".")
        if frame_code_components[0].startswith("_"):
            frame_code_components[0] = frame_code_components[0][1:]
        num_components = len(frame_code_components)
        if num_components == 1:
            self._curr_frame_code_category, self._curr_frame_code_keyword = (frame_code_components[0], None)
            return
        if num_components == 2:
            self._curr_frame_code_category, self._curr_frame_code_keyword = frame_code_components
            return
        raise

    def _process_loop(self):
        if self._curr_token_value != "":
            raise CIFParsingError
        return

    def _process_name(self):
        name_components = self._curr_token_value.split(".")
        num_components = len(name_components)
        if num_components != 2:
            raise CIFParsingError
        self._curr_data_name_category, self._curr_data_name_keyword = name_components
        return

    def _process_value(self):
        self._curr_data_value = self._curr_token_value
        return

    def _process_value_quoted(self):
        self._curr_data_value = self._curr_token_value.strip()
        return

    def _process_value_double_quoted(self):
        self._curr_data_value = self._curr_token_value.strip()
        self._curr_token_type = Token.VALUE_QUOTED
        return

    def _process_value_field(self):
        lines = self._curr_token_value.strip().splitlines()
        lines_processed = [lines[0].strip()] + [
            f"{' ' if line.startswith(' ') else ''}{line.strip()}" for line in lines[1:]
        ]
        self._curr_data_value = "".join(lines_processed)
        return

    def _add_data_item(self):
        data_category_dict = self._curr_data_category_dict
        if self._curr_data_name_keyword in data_category_dict.keys():
            raise CIFParsingError(self, CIFParsingErrorType.DUPLICATE)
        data_category_dict[self._curr_data_name_keyword] = self._curr_data_value
        self._add_data(
            block_code=self._curr_block_code,
            frame_code_category=self._curr_frame_code_category,
            frame_code_keyword=self._curr_frame_code_keyword,
            data_name_category=self._curr_data_name_category,
            data_name_keyword=self._curr_data_name_keyword,
            data_value=self._curr_data_value,
            value_dimension=0
        )
        return

    def _add_data(
            self,
            block_code: str,
            frame_code_category: str,
            frame_code_keyword: str,
            data_name_category: str,
            data_name_keyword: str,
            data_value: Union[str, list],
            value_dimension: int
    ):
        self._block_codes.append(block_code)
        self._frame_code_categories.append(frame_code_category)
        self._frame_code_keywords.append(frame_code_keyword)
        self._data_name_categories.append(data_name_category)
        self._data_name_keywords.append(data_name_keyword)
        self._data_values.append(data_value)
        self._value_dimensions.append(value_dimension)
        return

    def _initialize_loop(self):
        self._loop_category = self._curr_data_name_category
        self._loop_keywords = [self._curr_data_name_keyword]
        return

    def _add_loop_keyword(self):
        if self._curr_data_name_category != self._loop_category:
            raise CIFParsingError()
        if self._curr_data_name_keyword in self._loop_keywords:
            raise CIFParsingError
        self._loop_keywords.append(self._curr_data_name_keyword)

        return

    def _register_and_fill_loop(self):
        self._register_loop()
        self._fill_loop_value()
        return

    def _register_loop(self):
        num_columns = len(self._loop_keywords)
        dimension = 2 if num_columns >= 2 else 1
        loop_value_lists = list()
        data_category_dict = self._curr_data_category_dict
        for loop_keyword in self._loop_keywords:
            value_list = data_category_dict.setdefault(loop_keyword, list())
            if value_list != list():
                raise CIFParsingError(self, CIFParsingErrorType.DUPLICATE)
            loop_value_lists.append(value_list)
            self._add_data(
                block_code=self._curr_block_code,
                frame_code_category=self._curr_frame_code_category,
                frame_code_keyword=self._curr_frame_code_keyword,
                data_name_category=self._curr_data_name_category,
                data_name_keyword=loop_keyword,
                data_value=value_list,
                value_dimension=dimension
            )
        self._loop_value_lists = itertools.cycle(loop_value_lists)
        self._loop_value_lists_idx = itertools.cycle(range(num_columns))
        return

    def _fill_loop_value(self):
        next(self._loop_value_lists).append(self._curr_data_value)
        next(self._loop_value_lists_idx)
        return

    def _finalize_loop(self):
        if next(self._loop_value_lists_idx) != 0:
            raise CIFParsingError(self, error_type=CIFParsingErrorType.TABLE_INCOMPLETE)
        return

    def _finalize(self):
        if self._curr_state in (State.IN_LOOP_VALUE, State.IN_SAVE_LOOP_VALUE):
            self._finalize_loop()
        elif self._curr_state not in (State.IN_DATA, State.IN_SAVE):
            raise CIFParsingError

        df = pd.DataFrame(
            dict(
                block_code=self._block_codes,
                frame_code_category=self._frame_code_categories,
                frame_code_keyword=self._frame_code_keywords,
                data_name_category=self._data_name_categories,
                data_name_keyword=self._data_name_keywords,
                data_value=self._data_values,
                data_value_dimension=self._value_dimensions,
            ),
            # schema=dict(
            #     block_code=pl.Utf8,
            #     frame_code_category=pl.Utf8,
            #     frame_code_keyword=pl.Utf8,
            #     data_name_category=pl.Utf8,
            #     data_name_keyword=pl.Utf8,
            #     data_value=pl.Object,
            #     data_value_dimension=pl.UInt8,
            # )
        )
        return _filestruct.CIFFile(dic=self._cif_dict, df=df)


class CIFParsingErrorType(enum.Enum):
    DUPLICATE = enum.auto()
    TABLE_INCOMPLETE = enum.auto()


class CIFParsingError(Exception):

    def __init__(self, parser_instance: CIFParser, error_type: CIFParsingErrorType):
        self.parser_instance = parser_instance
        self.error_type = error_type
        self.error_msg = getattr(CIFParsingError, f"_{error_type.name.lower()}")(parser_instance)
        super().__init__(self.error_msg)
        return

    @staticmethod
    def _duplicate(parser: CIFParser):
        return (
            f"Duplicated data item in data block '{parser._curr_block_code}': The data name "
            f"'_{parser._curr_data_name_category}.{parser._curr_data_name_keyword}' "
            f"(i.e. category '{parser._curr_data_name_category}', "
            f"keyword '{parser._curr_data_name_keyword}') "
            f"is already registered with a data value of "
            f"'{parser._curr_data_category_dict[parser._curr_data_name_keyword]}', "
            f"but a second declaration with a value of '{parser._curr_data_value}' "
            f"was encountered at position {parser._curr_match.start()} of the file."
        )

    @staticmethod
    def _table_incomplete(parser: CIFParser):
        return