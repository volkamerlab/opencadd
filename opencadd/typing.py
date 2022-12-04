from typing import TypeAlias, Union, Sequence
from pathlib import Path
from numpy import ndarray

PathLike: TypeAlias = Union[str, Path]
ArrayLike: TypeAlias = Union[Sequence, ndarray]