"""
Custom decorators used by different functions within the openCADD library.
"""

from typing import Callable, Optional, Union, Tuple, Type, NamedTuple
from functools import wraps
import time


__author__ = "Armin Ariamajd"


class RetryConfig(NamedTuple):
    """
    Configuration for the `retry_on_exception` decorator.

    Attributes
    ----------
    num_tries : int, default: 3
        Maximum number of times the decorated function will be called
        before an exception is reraised.
    sleep_time_init : float, default: 1
        Amount of time (in seconds) to wait before two function calls.
    sleep_time_scale : float, default: 3
        Scaling factor for `sleep_time_init`.
        This can be used to scale down/up the waiting time between function calls.
        After the n-th function call, the waiting time will be equal to:
        `sleep_time_init` * `sleep_time_scale` ^ (n - 1).
    """
    num_tries: int = 3
    sleep_time_init: float = 1
    sleep_time_scale: float = 3


def retry_on_exception(
    function: Optional[Callable] = None,
    *,
    config: RetryConfig = RetryConfig(),
    catch: Union[Type[Exception], Tuple[Type[Exception]]] = Exception,
) -> Callable:
    """
    Decorator to retry a function call for a given number of times (while waiting for a certain
    amount of time between calls), when one of the given exceptions is raised.

    Parameters
    ----------
    function : callable
        The function to be decorated.
    config : opencadd._decorator.RetryConfig, default: RetryConfig(3, 1, 3)
        Retry configuration.
    catch : Type[Exception] | Tuple[Type[Exception]], default: Exception
        Exception type(s) that will be ignored during the retries. 
        All other exceptions will be raised immediately.

    Returns
    -------
    callable
        Decorated function.
    """
    if not isinstance(config.num_tries, int) or config.num_tries < 1:
        raise ValueError("`num_tries` must be a positive integer.")

    def retry_decorator(func):
        @wraps(func)
        def retry_wrapper(*args, **kwargs):
            curr_sleep_seconds = config.sleep_time_init
            for try_count in range(config.num_tries):
                try:
                    return func(*args, **kwargs)
                except catch as e:
                    if try_count == config.num_tries - 1:
                        raise e
                    time.sleep(curr_sleep_seconds)
                    curr_sleep_seconds *= config.sleep_time_scale
        return retry_wrapper
    return retry_decorator if function is None else retry_decorator(function)
