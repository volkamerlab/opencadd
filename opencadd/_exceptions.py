"""
General exceptions raised within the openCADD library.
"""

from typing import Any, NoReturn, Type, Optional, Sequence, Tuple, Union, Literal
from functools import partial
import operator
import numpy as np
import jax.numpy as jnp

from opencadd import _typing


__author__ = "Armin Ariamajd"


def error_string(
        parent_name: str,
        param_name: str,
        param_arg: Any,
        expectation: str,
        catch: str,
        title: Optional[str] = None
) -> str:
    err = (
        f"Parameter `{param_name}` of `{parent_name}` expects {expectation}, "
        f"but {catch}. Input was: {param_arg}."
    )
    return err if title is None else f"{title}:\n{err}"


def raise_for_type(
        parent_name: str,
        *args: Tuple[str, Any, Type]
) -> NoReturn:
    """
    Check the type of a series of input arguments for a function/method, and raise a TypeError if necessary.

    Parameters
    ----------
    parent_name : str
        Name of the function/method running the type checking.
    *args : tuple of 3-tuples of (str, any, type)
        A number of 3-tuples, each corresponding to one parameter in the function/method.
        Each 3-tuple contains the name of the parameter, its input argument,
        and the expected type, respectively.

    Raises
    ------
    TypeError
    """
    for param_name, param_arg, expected_type in args:
        if not isinstance(param_arg, expected_type):
            raise TypeError(
                error_string(
                    parent_name=parent_name,
                    param_name=param_name,
                    param_arg=param_arg,
                    expectation=f"an input of type {expected_type}",
                    catch=f"the type of input argument was {type(param_arg)}"
                )
            )


def raise_literal(
        parent_name: str,
        *args: Tuple[str, Any, Sequence[Any]]
) -> NoReturn:
    for param_name, param_arg, expected_vals in args:
        if param_arg not in expected_vals:
            raise ValueError(
                error_string(
                    parent_name=parent_name,
                    param_name=param_name,
                    param_arg=param_arg,
                    expectation=f"a value from the set {expected_vals}",
                    catch=f"the input argument was {param_arg}"
                )
            )



def raise_array(
        parent_name: str,
        param_name: str,
        array: Union[np.ndarray, jnp.ndarray],
        **kwargs
        # ndim_lt: Optional[int] = None,
        # ndim_eq: Optional[int] = None,
        # ndim_gt: Optional[int] = None,
        # size: Optional[int] = None,
        # shape: Optional[Tuple[Union[int, Sequence[int]], Optional[slice]]] = None,
        # dtypes: Optional[Sequence[np.dtype]] = None,
        # ge: Optional[float] = None,
        # gt: Optional[float] = None,
        # le: Optional[float] = None,
        # lt: Optional[float] = None,
):
    """
    Check the specifications of an array and raise an error if necessary.

    Parameters
    ----------
    parent_name : str
        Name of the function/method running the checking.
    param_name : str
        Name of the input parameter corresponding to the array.
    array : numpy.ndarray or jax.numpy.ndarray
        The input argument.
    kwargs

    Returns
    -------

    """
    error_string_partial = partial(
        error_string, parent_name=parent_name, param_name=param_name, param_arg=array
    )

    def ndim(expected_value, compare_func, compare_text):
        if compare_func(array.ndim, expected_value):
            raise ValueError(
                error_string_partial(
                    expecatition=f"an array with {compare_text}{expected_value} dimensions",
                    catch=f"the input array was {array.ndim}-dimensional",
                    title="Array Dimension Mismatch"
                )
            )

    def size(expected_value, compare_func, compare_text):
        if compare_func(array.size, expected_value):
            raise ValueError(
                error_string_partial(
                    expecatition=f"an array with {compare_text}{expected_value} elements",
                    catch=f"the size of input array was {array.size}",
                    title="Array Size Mismatch"
                )
            )

    def dtype(types):
        if isinstance(types, str) or not isinstance(types, Sequence):
            types = [types]
        if "real" in types:
            types.remove("real")
            types.extend([np.integer, np.floating])
        if not np.any([np.issubdtype(array.dtype, datatype) for datatype in types]):
            raise TypeError(
                error_string_partial(
                    expectation=f"an array of type(s) {types}",
                    catch=f"the datatype of input array was {array.dtype}"
                )
            )

    def shape():
        pass

    def val(expected_value, compare_func, compare_text):
        has_condition = compare_func(array, expected_value)
        if np.any(has_condition):
            raise ValueError(
                error_string_partial(
                    expecatition=f"an array with values {compare_text}{expected_value}",
                    catch=f"the input array had values {array[has_condition]} at positions {np.argwhere(has_condition)}",
                    title="Array Value Mismatch"
                )
            )

    maps = dict(
        ndim_lt=partial(ndim, compare_func=operator.ge, compare_text="less than "),
        ndim_eq=partial(ndim, compare_func=operator.ne, compare_text=""),
        ndim_gt=partial(ndim, compare_func=operator.le, compare_text="more than "),
        size_lt=partial(size, compare_func=operator.ge, compare_text="less than "),
        size_eq=partial(size, compare_func=operator.ne, compare_text=""),
        size_gt=partial(size, compare_func=operator.le, compare_text="more than "),
        val_lt=partial(val, compare_func=operator.ge, compare_text="less than "),
        val_le=partial(val, compare_func=operator.gt, compare_text="less than, or equal to, "),
        val_eq=partial(val, compare_func=partial(np.isin, invert=True), compare_text="in "),
        val_ge=partial(val, compare_func=operator.lt, compare_text="greater than, or equal to, "),
        val_gt=partial(val, compare_func=operator.le, compare_text="greater than "),
        dtype=dtype

    )

    for key, val in kwargs.items():
        maps[key](val)

    return


def check_number(
        number,
        dtypes: Optional[Union[np.dtype, Sequence[np.dtype], Literal["real"]]] = None,
        ge: Optional[float] = None,
        gt: Optional[float] = None,
        le: Optional[float] = None,
        lt: Optional[float] = None,
):
    if dtypes is None:
        dtypes = [np.number]
    elif isinstance(dtypes, str) and dtypes == "real":
        dtypes = [np.integer, np.floating]
    elif not isinstance(dtypes, Sequence):
        dtypes = [dtypes]
    num_type = type(number)
    if not np.any([np.issubdtype(num_type, dtype) for dtype in dtypes]):
        raise ValueError()
    if ge is not None and number < ge:
        raise ValueError()
    elif gt is not None and number <= gt:
        raise ValueError()
    if le is not None and number > le:
        raise ValueError()
    elif lt is not None and number >= lt:
        raise ValueError()
    return

