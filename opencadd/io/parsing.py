"""
General functions and routines for data parsing.
"""

from typing import Tuple, Dict, Callable, Sequence, NamedTuple, Union, Any, List, Optional
import numpy as np
import numpy.typing as npt
import pandas as pd


class Column:
    def __init__(
            self,
            interval: Tuple[int, int],
            parser: Optional[Callable[[np.ndarray], Any]] = None
    ):
        self.interval = interval
        self.parser = (lambda array: array) if parser is None else parser
        self.slice = slice(*interval)
        self.length = interval[1] - interval[0]
        return


class MultiColumn:
    def __init__(
            self,
            intervals: Sequence[Tuple[int, int]],
            parser: Optional[Callable[[np.ndarray], Any]] = None,
    ):
        self.intervals = np.asarray(intervals, dtype=int).reshape(-1, 2)
        self.parser = (lambda array: array) if parser is None else parser
        self.columns = np.concatenate(
            [np.arange(interval[0], interval[1]) for interval in self.intervals]
        )
        self.length = self.intervals[0, 1] - self.intervals[0, 0]

    @property
    def is_multi(self) -> bool:
        return self.ranges.shape[0] > 1

    @property
    def count_columns(self) -> int:
        return self.ranges.shape[1]


def extract_column_from_string_array(array: np.ndarray, char_range: Tuple[int, int]) -> np.ndarray:
    """
    From an array of strings representing the rows of a table, extract the sub-string between
    the given start and end column (i.e. character) numbers, in a fast vectorized fashion.

    Parameters
    ----------
    array : numpy.ndarray
        A 1-dimensional array of strings, where each string corresponds to a row of a table.
    char_range : tuple[int, int]
        Range (start, end) of characters in rows, which correspond to a column.

    Returns
    -------
    column : numpy.ndarray
        A 1-dimensional array of strings, where each string now only contains the sub-string in
        the given range, i.e. the extracted column.
    """
    # Create a view with Char data type, reshape array to have the same shape as input,
    # and take only the characters in the given range
    view = array.view(dtype=(str, 1)).reshape(array.size, -1)[:, char_range[0]:char_range[1]]
    # Create array from view buffer
    return np.frombuffer(view.tobytes(), dtype=(str, char_range[1] - char_range[0]))


def extract_column(
        char_array: np.ndarray,
        column: Union[Column, MultiColumn],
) -> Any:
    if isinstance(column, Column):
        return column.parser(
            char_array[:, column.slice].view(dtype=(str, column.length)).reshape(-1)
        )
    cols = np.frombuffer(
        char_array[:, column.columns].tobytes(), dtype=(str, column.length)
    ).reshape(char_array.shape[0], -1)
    #char_array[:, column.columns].view(dtype=(str, column.length))
    return column.parser(cols)


def extract_columns(
        char_array: np.ndarray,
        columns: Union[Sequence[Column], Dict[Any, Column]],
        to_dataframe: bool = False,
) -> Union[List[np.ndarray], Dict[Any, np.ndarray], pd.DataFrame]:

    data = (
        dict(zip(columns.keys(), [extract_column(char_array, col) for col in columns.values()]))
        if isinstance(columns, dict) else [extract_column(char_array, col) for col in columns]
    )
    return pd.DataFrame(data) if to_dataframe else data
