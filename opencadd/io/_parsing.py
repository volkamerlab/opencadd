"""
General functions and routines for data parsing.
"""

from typing import Tuple, Dict, Callable, Sequence, NamedTuple, Union, Any, List, Optional
import numpy as np
import numpy.typing as npt
import pandas as pd



def string_array_to_char_table(array: Sequence[str], view: bool = True) -> np.ndarray:
    arr = np.asarray(array) if view else np.array(array)
    return arr.view(dtype=(str, 1)).reshape(arr.size, -1)


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


def extract_columns_by_interval(
        char_table: np.ndarray,
        interval: Tuple[int, int],
        strip: bool = True
) -> np.ndarray:
    column = char_table[:, interval[0]:interval[1]].view(dtype=(str, interval[1] - interval[0])).reshape(-1)
    return np.char.strip(column) if strip else column


def extract_columns_by_index(
        char_table: np.ndarray,
        indices: Sequence[int],
        column_len: int,
        strip: bool = True
) -> np.ndarray:
    """"""
    column = np.asarray(char_table[:, indices], order="C").view(dtype=(str, column_len))
    return np.char.strip(column) if strip else column


    # if isinstance(column, Column):
    #     return column.parser(
    #         char_table[:, column.slice].view(dtype=(str, column.length)).reshape(-1)
    #     )
    # cols = np.frombuffer(
    #     char_table[:, column.columns].tobytes(), dtype=(str, column.length)
    # ).reshape(char_table.shape[0], -1)
    # #char_array[:, column.columns].view(dtype=(str, column.length))
    # return column.parser(cols)


# def extract_columns(
#         char_array: np.ndarray,
#         columns: Union[Sequence[Column], Dict[Any, Column]],
#         to_dataframe: bool = False,
# ) -> Union[List[np.ndarray], Dict[Any, np.ndarray], pd.DataFrame]:
#
#     data = (
#         dict(zip(columns.keys(), [extract_column(char_array, col) for col in columns.values()]))
#         if isinstance(columns, dict) else [extract_column(char_array, col) for col in columns]
#     )
#     return pd.DataFrame(data) if to_dataframe else data
