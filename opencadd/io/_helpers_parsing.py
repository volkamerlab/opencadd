"""
General functions and routines for data parsing.
"""

from typing import Tuple, Dict, Callable, Sequence, Any
import numpy as np
import numpy.typing as npt
import pandas as pd


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


def dataframe_from_string_array(
        array: np.ndarray,
        column_data: Dict[str, Tuple[Tuple[int, int], npt.DTypeLike]],
        strip: str = "",
        transforms: Sequence[Tuple[str, Any]] = None,
) -> pd.DataFrame:
    df = pd.DataFrame()
    for col_name, (col_range, col_dtype) in column_data.items():
        column = extract_column_from_string_array(
                array=array,
                char_range=(col_range[0] - 1, col_range[1])
            )
        column = np.char.strip(

        ).astype(col_dtype)

    return