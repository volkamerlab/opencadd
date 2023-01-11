
from typing import IO
from opencadd.typing import FileLike
from pathlib import Path
import io


def filelike_to_filepath(file: FileLike) -> Path:
    """
    Turn a file-like object to IO.

    Parameters
    ----------
    file

    Returns
    -------

    """
    # if isinstance(file, (io.BufferedIOBase, io.BytesIO)):
    #     return file
    # if isinstance(file, (io.TextIOBase, io.StringIO)):
    #     return file.buffer
    #
    # if isinstance(file, str):
    #     filepath = Path(file)
    #     try:
    #         if filepath.is_file():
    #             with open(filepath, "rb") as f:
    #                 return f
    return


def filelike_to_data_string(file):
    if isinstance(file, Path):
        with open(file, "r") as f:
            return f.read()
    if isinstance(file, str):
        possible_path = Path(file)
        if possible_path.is_file():
            with open(possible_path, "r") as f:
                return f.read()
        return file
    if isinstance(file, bytes):
        return file.decode()
    return

