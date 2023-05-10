"""
Base functionalities to parse an mmCIF datafile or dictionary.
"""


from typing import Union, Iterator, List, Literal, NoReturn
import enum
import re
import itertools

import polars as pl

import opencadd._exceptions as _exception
from . import _filestruct


__author__ = "Armin Ariamajd"


TOKENIZER = re.compile(
    r"""(?xmi)  # `x` (cf. re.X) allows for writing the expression in multiple lines, with comments added; 
                # `m` (cf. re.M, re.MULTILINE) causes the pattern characters '^' and '$' to also match
                #  the beggining and end of each line, respectively 
                #  (in addition to matching the beggining and end of the whole string, respectively).
                # `i` (cf. re.I, re.IGNORECASE) performs case-insensitive matching.
    # The following creates different capturing groups (enumerated starting from 1), 
    #  each matching one token type. Notice the order of groups matters, 
    #  since the matching terminates with the first group match.
    ^;([\S\s]*?)(?:\r\n|\s)^;(?:(?=\s)|$)  # 1. Text field, i.e. a non-simple data value 
                                           #    bounded between two '\n;' characters.
    |(?:^|(?<=\s))\#(.*?)\r?$              # 2. Comment
    |(?:^|(?<=\s))(?:
      '(.*?)'                              # 3. Quoted data value
      |"(.*?)"                             # 4. Duble-quoted data value
      |_(\S*)                              # 5. Data name
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
    VALUE_FIELD = 1  # Will be changed to `VALUE` after processing by parser
    COMMENT = 2
    VALUE_QUOTED = 3  # Will be changed to `VALUE` after processing by parser
    VALUE_DOUBLE_QUOTED = 4  # Will be changed to `VALUE` after processing by parser
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
    SAVE_END = 15  # Will be added by the parser after processing `SAVE`


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


class CIFParsingErrorType(enum.Enum):
    """
    Types of errors that may occur during parsing.
    """
    DUPLICATE = 2
    TABLE_INCOMPLETE = 2
    BLOCK_CODE_EMPTY = 2
    BAD_TOKEN = 2
    RESERVED_TOKEN = 2
    UNEXPECTED_TOKEN = 2
    LOOP_HAS_NAME = 3
    INCOMPLETE_FILE = 2


class CIFParser:
    """
    CIF file parser.
    """

    def __init__(self):

        self._token_processor = {
            Token.DATA: self._process_block_code,
            Token.SAVE: self._process_frame_code,
            Token.LOOP: self._process_loop_directive,
            Token.NAME: self._process_data_name,
            Token.VALUE: self._process_data_value,
            Token.VALUE_QUOTED: self._process_data_value_quoted,
            Token.VALUE_DOUBLE_QUOTED: self._process_data_value_double_quoted,
            Token.VALUE_FIELD: self._process_data_value_text_field,
        }
        """
        A mapping between token type and corresponding processing method.
        """

        self._state_mapper = {
            (State.IN_FILE, Token.DATA): (self._do_nothing, State.JUST_IN_DATA),
            (State.IN_FILE, Token.COMMENT): (self._do_nothing, State.IN_FILE),
            (State.JUST_IN_DATA, Token.SAVE): (self._do_nothing, State.JUST_IN_SAVE),
            (State.JUST_IN_DATA, Token.LOOP): (self._do_nothing, State.JUST_IN_LOOP),
            (State.JUST_IN_DATA, Token.NAME): (self._do_nothing, State.IN_NAME),
            (State.JUST_IN_DATA, Token.COMMENT): (self._do_nothing, State.JUST_IN_DATA),
            (State.JUST_IN_SAVE, Token.LOOP): (self._do_nothing, State.JUST_IN_SAVE_LOOP),
            (State.JUST_IN_SAVE, Token.NAME): (self._do_nothing, State.IN_SAVE_NAME),
            (State.JUST_IN_SAVE, Token.COMMENT): (self._do_nothing, State.JUST_IN_SAVE),
            (State.JUST_IN_LOOP, Token.NAME): (self._initialize_loop, State.IN_LOOP_NAME),
            (State.JUST_IN_LOOP, Token.COMMENT): (self._do_nothing, State.JUST_IN_LOOP),
            (State.IN_NAME, Token.VALUE): (self._add_data_item, State.IN_DATA),
            (State.IN_NAME, Token.COMMENT): (self._do_nothing, State.IN_NAME),
            (State.JUST_IN_SAVE_LOOP, Token.NAME): (self._initialize_loop, State.IN_SAVE_LOOP_NAME),
            (State.JUST_IN_SAVE_LOOP, Token.COMMENT): (self._do_nothing, State.JUST_IN_SAVE_LOOP),
            (State.IN_SAVE_NAME, Token.VALUE): (self._add_data_item, State.IN_SAVE),
            (State.IN_SAVE_NAME, Token.COMMENT): (self._do_nothing, State.IN_SAVE_NAME),
            (State.IN_LOOP_NAME, Token.NAME): (self._add_loop_keyword, State.IN_LOOP_NAME),
            (State.IN_LOOP_NAME, Token.VALUE): (self._register_and_fill_loop, State.IN_LOOP_VALUE),
            (State.IN_LOOP_NAME, Token.COMMENT): (self._do_nothing, State.IN_LOOP_NAME),
            (State.IN_DATA, Token.DATA): (self._do_nothing, State.JUST_IN_DATA),
            (State.IN_DATA, Token.SAVE): (self._do_nothing, State.JUST_IN_SAVE),
            (State.IN_DATA, Token.LOOP): (self._do_nothing, State.JUST_IN_LOOP),
            (State.IN_DATA, Token.NAME): (self._do_nothing, State.IN_NAME),
            (State.IN_DATA, Token.COMMENT): (self._do_nothing, State.IN_DATA),
            (State.IN_SAVE_LOOP_NAME, Token.NAME): (self._add_loop_keyword, State.IN_SAVE_LOOP_NAME),
            (State.IN_SAVE_LOOP_NAME, Token.VALUE): (self._register_and_fill_loop, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_NAME, Token.COMMENT): (self._do_nothing, State.IN_SAVE_LOOP_NAME),
            (State.IN_SAVE, Token.SAVE_END): (self._do_nothing, State.IN_DATA),
            (State.IN_SAVE, Token.LOOP): (self._do_nothing, State.JUST_IN_SAVE_LOOP),
            (State.IN_SAVE, Token.NAME): (self._do_nothing, State.IN_SAVE_NAME),
            (State.IN_SAVE, Token.COMMENT): (self._do_nothing, State.IN_SAVE),
            (State.IN_LOOP_VALUE, Token.DATA): (self._finalize_loop, State.JUST_IN_DATA),
            (State.IN_LOOP_VALUE, Token.SAVE): (self._finalize_loop, State.JUST_IN_SAVE),
            (State.IN_LOOP_VALUE, Token.LOOP): (self._finalize_loop, State.JUST_IN_LOOP),
            (State.IN_LOOP_VALUE, Token.NAME): (self._finalize_loop, State.IN_NAME),
            (State.IN_LOOP_VALUE, Token.VALUE): (self._fill_loop_value, State.IN_LOOP_VALUE),
            (State.IN_LOOP_VALUE, Token.COMMENT): (self._do_nothing, State.IN_LOOP_VALUE),
            (State.IN_SAVE_LOOP_VALUE, Token.SAVE_END): (self._finalize_loop, State.IN_DATA),
            (State.IN_SAVE_LOOP_VALUE, Token.LOOP): (self._finalize_loop, State.JUST_IN_SAVE_LOOP),
            (State.IN_SAVE_LOOP_VALUE, Token.NAME): (self._finalize_loop, State.IN_SAVE_NAME),
            (State.IN_SAVE_LOOP_VALUE, Token.VALUE): (self._fill_loop_value, State.IN_SAVE_LOOP_VALUE),
            (State.IN_SAVE_LOOP_VALUE, Token.COMMENT): (self._do_nothing, State.IN_SAVE_LOOP_VALUE),
        }
        """
        A mapping between (current state, received token) and (action, resulting state).
        """

        self.file_content: str = None
        self._tokenizer: Iterator[re.Match] = None
        self.curr_state: State = None
        self.curr_match: re.Match = None
        self.curr_token_type: Token = None
        self.curr_token_value: str = None
        self.curr_block_code: str = None
        self.curr_frame_code_category: str = None
        self.curr_frame_code_keyword: str = None
        self.curr_data_name_category: str = None
        self.curr_data_name_keyword: str = None
        self.curr_data_value: str = None
        self.errors: List[CIFParsingError] = None
        return

    def _reset_state(self, content: str) -> NoReturn:
        self.file_content: str = content
        self._tokenizer: Iterator[re.Match] = TOKENIZER.finditer(self.file_content)
        self.curr_state = State.IN_FILE
        self.curr_match: re.Match = None
        self.curr_token_type: Token = None
        self.curr_token_value: str = None
        self.curr_block_code: str = None
        self.curr_frame_code_category: str = None
        self.curr_frame_code_keyword: str = None
        self.curr_data_name_category: str = None
        self.curr_data_name_keyword: str = None
        self.curr_data_value: str = None
        self.errors: List[CIFParsingError] = list()
        return

    def parse(self, content: str):
        _exception.raise_for_type(
            f"{self.__class__.__name__}.parse",
            ("content", content, str)
        )
        self._reset_state(content=content)
        for match in self._tokenizer:
            self.curr_match = match
            self.curr_token_type = Token(match.lastindex)
            self.curr_token_value = match.group(match.lastindex)
            self._token_processor.get(self.curr_token_type, self._do_nothing)()
            update_func, new_state = self._state_mapper.get(
                (self.curr_state, self.curr_token_type),
                (self._wrong_token, self.curr_state)
            )
            update_func()
            self.curr_state = new_state
        return self._finalize()

    def _finalize(self) -> NoReturn:
        if self.curr_state in (State.IN_LOOP_VALUE, State.IN_SAVE_LOOP_VALUE):
            self._finalize_loop()
        elif self.curr_state not in (State.IN_DATA, State.IN_SAVE):
            self._register_error(CIFParsingErrorType.INCOMPLETE_FILE)
        return

    def _process_block_code(self) -> NoReturn:
        self.curr_block_code = self.curr_token_value
        return

    def _process_frame_code(self) -> NoReturn:
        if self.curr_token_value == "":
            self.curr_token_type = Token.SAVE_END
            self.curr_frame_code_category = None
            self.curr_frame_code_keyword = None
            return
        if self.curr_token_value.startswith("_"):
            self.curr_token_value = self.curr_token_value[1:]
        frame_code_components = self.curr_token_value.split(".", 1)
        self.curr_frame_code_category = frame_code_components[0]
        try:
            self.curr_frame_code_keyword = frame_code_components[1]
        except IndexError:
            self.curr_frame_code_keyword = None
        return

    def _process_loop_directive(self) -> NoReturn:
        if self.curr_token_value != "":
            self._register_error(CIFParsingErrorType.LOOP_HAS_NAME)
        return

    def _process_data_name(self) -> NoReturn:
        name_components = self.curr_token_value.split(".", 1)
        self.curr_data_name_category = name_components[0]
        try:
            self.curr_data_name_keyword = name_components[1]
        except IndexError:
            self.curr_data_name_keyword = None
        return

    def _process_data_value(self) -> NoReturn:
        self.curr_data_value = self.curr_token_value
        return

    def _process_data_value_quoted(self) -> NoReturn:
        self.curr_data_value = self.curr_token_value.strip()
        self.curr_token_type = Token.VALUE
        return

    def _process_data_value_double_quoted(self) -> NoReturn:
        self.curr_data_value = self.curr_token_value.strip()
        self.curr_token_type = Token.VALUE
        return

    def _process_data_value_text_field(self) -> NoReturn:
        lines = self.curr_token_value.strip().splitlines()
        lines_processed = [lines[0].strip()] + [
            f"{' ' if line.startswith(' ') else ''}{line.strip()}" for line in lines[1:]
        ]
        self.curr_data_value = "".join(lines_processed)
        self.curr_token_type = Token.VALUE
        return

    def _wrong_token(self) -> NoReturn:
        if self.curr_token_type == Token.BAD_TOKEN:
            self._register_error(CIFParsingErrorType.BAD_TOKEN)
        elif self.curr_token_type in [Token.STOP, Token.GLOBAL, Token.FRAME_REF, Token.BRACKETS]:
            self._register_error(CIFParsingErrorType.RESERVED_TOKEN)
        else:
            self._register_error(CIFParsingErrorType.UNEXPECTED_TOKEN)
        return

    def _register_error(self, error_type: CIFParsingErrorType) -> NoReturn:
        """
        Given an error type, raise it as a `CIFParsingError` or post a warning message,
        depending on the level of `strictness` and the error level.

        Parameters
        ----------
        error_type : CIFParsingErrorType
            Error type.
        raise_level : {1, 2, 3}
            Minimum strictness level where the error should be raised as an exception.

        Raises
        -------
        CIFParsingError
        """
        self.errors.append(CIFParsingError(parser_instance=self, error_type=error_type))
        return

    def _do_nothing(self):
        return

    def _add_data_item(self):
        ...

    def _initialize_loop(self):
        ...

    def _add_loop_keyword(self):
        ...

    def _register_and_fill_loop(self):
        ...

    def _fill_loop_value(self):
        ...

    def _finalize_loop(self):
        ...


class CIFParsingError(Exception):

    def __init__(self, parser_instance: CIFParser, error_type: CIFParsingErrorType):

        self.state = parser_instance.curr_state
        self.match = parser_instance.curr_match
        self.token_type = parser_instance.curr_token_type
        self.token_value = parser_instance.curr_token_value
        self.block_code = parser_instance.curr_block_code
        self.frame_code_category = parser_instance.curr_frame_code_category
        self.frame_code_keyword = parser_instance.curr_frame_code_keyword
        self.data_name_category = parser_instance.curr_data_name_category
        self.data_name_keyword = parser_instance.curr_data_name_keyword
        self.data_value = parser_instance.curr_data_value

        self.parser_instance = parser_instance
        self.error_type = error_type
        self.error_msg = self.error_message(parser_instance, error_type)
        super().__init__(self.error_msg)
        return

    @staticmethod
    def error_message(parser_instance: CIFParser, error_type: CIFParsingErrorType):
        return getattr(CIFParsingError, f"_{error_type.name.lower()}")(parser_instance)

    @staticmethod
    def _duplicate(parser: CIFParser):
        return (
            f"Duplicated data item in data block '{parser.curr_block_code}': The data name "
            f"'_{parser.curr_data_name_category}.{parser.curr_data_name_keyword}' "
            f"(i.e. category '{parser.curr_data_name_category}', "
            f"keyword '{parser.curr_data_name_keyword}') "
            f"is already registered with a data value of "
            f"'{parser._curr_data_category_dict[parser.curr_data_name_keyword]}', "
            f"but a second declaration with a value of '{parser.curr_data_value}' "
            f"was encountered at position {parser.curr_match.start()} of the file."
        )

    @staticmethod
    def _table_incomplete(parser: CIFParser):
        return

    @staticmethod
    def _bad_token(parser: CIFParser):
        return f"Bad token: got {parser.curr_token_value} at position {parser.curr_match.span()}"

    @staticmethod
    def _unexpected_token(parser: CIFParser):
        return (
            f"Token out of order: parser is in state {parser.curr_state} "
            f"and expects a token from {EXPECTED_TOKENS[parser.curr_state]}, "
            f"but received a {parser.curr_token_type}."
        )


class CIFToDictParser(CIFParser):
    """
    CIF file parser.
    """

    def __init__(self):

        super().__init__()

        self._cif_dict_horizontal: dict = dict()

        self._block_codes: List[str] = list()
        self._frame_code_categories: List[Union[str, None]] = list()
        self._frame_code_keywords: List[Union[str, None]] = list()
        self._data_name_categories: List[str] = list()
        self._data_name_keywords: List[str] = list()
        self._data_values: List[List[str]] = list()
        self._loop_id: List[int] = list()

        self._loop_value_lists: itertools.cycle = None
        self._loop_value_lists_idx: itertools.cycle = None

        self._curr_loop_id: int = 0
        self._curr_loop_columns: List[List[str]] = list()
        return

    @property
    def _curr_data_block_dict(self) -> dict:
        return self._cif_dict_horizontal.setdefault(self.curr_block_code, dict())

    @property
    def _curr_save_frame_category_dict(self) -> dict:
        return self._curr_data_block_dict.setdefault(self.curr_frame_code_category, dict())

    @property
    def _curr_save_frame_keyword_dict(self) -> dict:
        return self._curr_save_frame_category_dict.setdefault(self.curr_frame_code_keyword, dict())

    @property
    def _curr_data_category_dict(self) -> dict:
        return self._curr_save_frame_keyword_dict.setdefault(self.curr_data_name_category, dict())

    # @property
    # def _curr_data_keyword_list(self) -> list:
    #     return self._curr_data_category_dict.setdefault(self.curr_data_name_keyword, list())

    def _add_data_item(self):
        # data_value_list = self._curr_data_keyword_list
        # if len(data_value_list) != 0:
        #     self._raise_or_warn(CIFParsingErrorType.DUPLICATE)
        #data_value_list.append(self.curr_data_value)
        self._add_data(data_value=[self.curr_data_value], loop_id=0)
        return

    def _initialize_loop(self):
        self._curr_loop_id += 1
        self._curr_loop_columns = list()
        self._add_loop_keyword()
        return

    def _add_loop_keyword(self):
        # column = self._curr_data_keyword_list
        # if len(column) != 0:
        #     self._raise_or_warn(CIFParsingErrorType.DUPLICATE)
        new_column = list()
        self._curr_loop_columns.append(new_column)
        self._add_data(data_value=new_column, loop_id=self._curr_loop_id)
        return

    def _register_and_fill_loop(self):
        self._register_loop()
        self._fill_loop_value()
        return

    def _register_loop(self):
        self._loop_value_lists = itertools.cycle(self._curr_loop_columns)
        self._loop_value_lists_idx = itertools.cycle(range(len(self._curr_loop_columns)))
        return

    def _fill_loop_value(self):
        next(self._loop_value_lists).append(self.curr_data_value)
        next(self._loop_value_lists_idx)
        return

    def _finalize_loop(self):
        if next(self._loop_value_lists_idx) != 0:
            self._register_error(CIFParsingErrorType.TABLE_INCOMPLETE)
        return

    def _finalize(self):
        super()._finalize()

        df = pl.DataFrame(
            dict(
                block_code=self._block_codes,
                frame_code_category=self._frame_code_categories,
                frame_code_keyword=self._frame_code_keywords,
                data_name_category=self._data_name_categories,
                data_name_keyword=self._data_name_keywords,
                data_value=self._data_values,
                loop_id=self._loop_id,
            ),
            dict(
                block_code=pl.Utf8,
                frame_code_category=pl.Utf8,
                frame_code_keyword=pl.Utf8,
                data_name_category=pl.Utf8,
                data_name_keyword=pl.Utf8,
                data_value=pl.List(pl.Utf8),
                loop_id=pl.UInt32,
            )
        )
        return CIFFileValidator(df=df, errors=self.errors)

    def _add_data(
            self,
            data_value: Union[str, list],
            loop_id: int
    ):
        self._curr_data_category_dict[self.curr_data_name_keyword] = data_value

        self._block_codes.append(self.curr_block_code)
        self._frame_code_categories.append(self.curr_frame_code_category)
        self._frame_code_keywords.append(self.curr_frame_code_keyword)
        self._data_name_categories.append(self.curr_data_name_category)
        self._data_name_keywords.append(self.curr_data_name_keyword)
        self._data_values.append(data_value)
        self._loop_id.append(loop_id)
        return


class CIFFileValidator:

    def __init__(self, df: pl.DataFrame, errors: List[CIFParsingError]):
        self._df = df
        self._df_ids = df.select(pl.exclude(["data_value", "loop_id"]))
        self._errors = errors
        self.cif_file = _filestruct.DDL2CIFFile(df=df.select(pl.exclude("loop_id")))
        return

    def validate(self, raise_: bool = True):
        pre_validated = len(self._errors) == 0
        post_validated = self.validate_post()
        validated = pre_validated and post_validated
        if validated:
            return True
        if raise_:
            raise CIFParsingError()
        return False

    def validate_post(self, ddl_version: Literal[1, 2] = 2):
        empty_data_id = not self.has_empty_data_identifier()
        duplicated = not self.has_duplicated_address()
        table_columns = self.table_columns_have_same_length()
        cat_columns = self.data_category_columns_have_same_length()
        null = not self.has_null()
        loops = self.loops_have_same_category()
        ddl = self.is_ddl2_conformant() if ddl_version == 2 else self.is_ddl1_conformant()
        return all([empty_data_id, duplicated, table_columns, cat_columns, null, loops, ddl])

    def has_empty_data_identifier(
            self, dim: Literal["elem", "col", "row", "df"] = "df"
    ) -> Union[pl.DataFrame, pl.Series, bool]:
        if dim == "df":
            return self._df_ids.select((pl.any(pl.all() == "")).any())[0, 0]
        if dim == "row":
            return self._df_ids.select(pl.any(pl.all() == ""))
        if dim == "col":
            return self._df_ids.select((pl.all() == "").any())
        if dim == "elem":
            return self._df_ids.select(pl.all() == "")
        raise ValueError()

    def has_duplicated_address(self, dim: Literal["row", "df"] = "df"):
        mask = self._df_ids.is_duplicated()
        if dim == "row":
            return mask
        if dim == "df":
            return mask.any()
        raise ValueError()

    def has_null(self, per_row: bool = False):
        series: pl.Series = self._df.select(
            pl.any(
                pl.exclude(["frame_code_category", "frame_code_keyword", "data_name_keyword"]).is_null()
            ).alias("mask")
        )["mask"]
        if per_row:
            return series
        return series.any()

    def table_columns_have_same_length(self, per_table: bool = False):
        df_per_loop = self._df.with_columns(
            pl.col("data_value").arr.lengths().alias("list_lengths")
        ).groupby(
            "loop_id"
        ).agg(
            (pl.col("list_lengths").n_unique() == 1).alias("has_same_length")
        )
        if per_table:
            return df_per_loop
        return df_per_loop["has_same_length"].all()

    def data_category_columns_have_same_length(self, per_category: bool = False):
        df_per_category = self._df.with_columns(
            pl.col("data_value").arr.lengths().alias("list_lengths")
        ).groupby(
            ["block_code", "frame_code_category", "frame_code_keyword", "data_name_category"]
        ).agg(
            (pl.col("list_lengths").n_unique() == 1).alias("has_same_length")
        )
        if per_category:
            return df_per_category
        return df_per_category["has_same_length"].all()

    def loops_have_same_category(self, per_loop: bool = False):
        df_per_loop = self._df.filter(
            pl.col("loop_id") > 0
        ).groupby(
            "loop_id"
        ).agg(
            (pl.col("data_name_category").n_unique() == 1).alias("has_same_category")
        )
        if per_loop:
            return df_per_loop
        return df_per_loop["has_same_category"].all()

    def is_ddl2_conformant(self, per_data_name: bool = False) -> [pl.Series, bool]:
        """
        Check whether data names conform to Dictionary Definition Language Version 2 (DDL2) standard.

        According to the DDL2 standard, data names must consist of a category name and a keyword name,
        separated by a '.' character. Notice that there is no constraint on frame codes, i.e. they can
        either have a category name and a keyword name like data names, or just a category name.

        Parameters
        ----------
        per_data_name : bool, default: False
            If True, a polars Series is returned, indicating whether each data name conforms (True) to
            the DDL2 standard or not (False). Otherwise, a single boolean value is returned indicating
            whether all data names are conformant (True) or not (False).

        Returns
        -------
        polars.Series or bool
        """
        # If the data name contains no '.', the parser assigns a None to the keyword name.
        #  Therefore, we can just check whether the `data_name_keyword` column has null values.
        series_per_data_name: pl.Series = self._df_ids.select(
            pl.col("data_name_keyword").is_not_null()
        )["data_name_keyword"]
        if per_data_name:
            return series_per_data_name
        return series_per_data_name.all()

    def is_ddl1_conformant(self, per_data_name: bool = False) -> [pl.Series, bool]:
        series_per_data_name: pl.Series = self._df.select(
            pl.all(pl.col(["frame_code_keyword", "data_name_keyword"]).is_null()).alias("mask")
        )["mask"]
        if per_data_name:
            return series_per_data_name
        return series_per_data_name.all()
