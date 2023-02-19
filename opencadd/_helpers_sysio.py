
from typing import IO, Union
from opencadd._typing import FileLike, PathLike, FileContentLike
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


def save_to_file(content: FileContentLike, filename: str, extension: str, path: PathLike) -> Path:
    dir_path = Path(path)
    dir_path.mkdir(parents=True, exist_ok=True)
    ext = extension if extension.startswith(".") else f".{extension}"
    fullpath = (dir_path/filename).with_suffix(ext)
    mode = "xb" if isinstance(content, bytes) else "xt"
    with open(fullpath, mode) as f:
        f.write(content)
    return fullpath.resolve()