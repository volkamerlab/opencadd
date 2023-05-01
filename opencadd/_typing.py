from typing import IO, Union, Sequence, Type, TypeAlias
from pathlib import Path
import numpy as np
import jax.numpy as jnp


__author__ = "Armin Ariamajd"


PathLike: TypeAlias = Union[str, Path]
FileContentLike: TypeAlias = Union[str, bytes, IO]
FileLike: TypeAlias = Union[PathLike, FileContentLike]
ArrayLike: TypeAlias = Union[Sequence, np.ndarray, jnp.ndarray]


def smallest_integer_dtype_for_range(min_val, max_val) -> Type[np.integer]:
    unsigned_types = [np.ubyte, np.ushort, np.uintc, np.uint, np.ulonglong]
    unsigned_maxes = [np.iinfo(dtype).max for dtype in unsigned_types]

    signed_types = [np.byte, np.short, np.intc, np.int_, np.longlong]
    signed_mins = [np.iinfo(dtype).min for dtype in signed_types]
    signed_maxes = [np.iinfo(dtype).max for dtype in signed_types]

    if (
            not np.issubdtype(type(min_val), np.integer)
            or not np.issubdtype(type(max_val), np.integer)
    ):
        raise ValueError(
            "Parameters `min_val` and `max_val` expect integer values. "
            f"Input was {min_val, max_val}"
        )
    if min_val > max_val:
        raise ValueError(
            f"`min_val` must be smaller than `max_val`. Input was: {min_val}, {max_val}"
        )

    if min_val >= 0:
        for idx_dtype, dtype in enumerate(unsigned_types):
            if max_val <= unsigned_maxes[idx_dtype]:
                return dtype
        raise ValueError(
            f"Bit overflow. Bounds for largest type ({[unsigned_types[-1]]}) is "
            f"[0, {unsigned_maxes[-1]}]. Given interval was [{min_val}, {max_val}]."
        )
    else:
        for idx_dtype, dtype in enumerate(signed_types):
            if min_val >= signed_mins[idx_dtype] and max_val <= signed_maxes[idx_dtype]:
                return dtype
        raise ValueError(
            f"Bit overflow. Bounds for largest type ({[signed_types[-1]]}) is "
            f"[{signed_mins[-1]}, {signed_maxes[-1]}]. Given interval was [{min_val}, {max_val}]."
        )
