"""
General exceptions raised within the openCADD library.
"""

from typing import Any, NoReturn, Type


__author__ = "Armin Ariamajd"


def raise_type_error(param_name: str, func_name: str, expected_type: Type, arg: Any) -> NoReturn:
    """
    Raise a TypeError.

    Parameters
    ----------
    param_name : str
        Name of the input parameter that received an inappropriate type.
    func_name : str
        Name of the function/method to which the parameter belongs.
    expected_type : str
        Expected type of the parameter.
    arg : Any
        Input argument of the parameter.

    Returns
    -------
    None

    Raises
    ------
    TypeError
    """
    raise TypeError(
        f"Parameter `{param_name}` of `{func_name}` expects an input of type {expected_type}, "
        f"but the type of input argument was {expected_type(arg)}. Input was: {arg}."
    )
