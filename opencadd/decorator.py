"""
Custom decorators used by other functions in the package.
"""

from typing import Callable, Optional, Union, Tuple, Type, NamedTuple
from functools import wraps
import time


class RetryConfig(NamedTuple):
    """
    Configuration for `retry_on_exception` decorator.

    Attributes
    ----------
    num_tries : int, optional, default: 3
        Maximum number of times the decorated function will be called before the exception is
        reraised.
    sleep_time_init : float, optional, default = 1
        Amount of time (in seconds) to wait before two function calls.
    sleep_time_scale : float, optional, default = 3
        Scaling factor for `sleep_time_init`. After the `n`-th function call, the waiting time
        will be equal to `sleep_time_init * (sleep_time_scale ** n-1)`. This can be used to
        scale down/up the waiting time between function calls.
    """
    num_tries: int = 3
    sleep_time_init: float = 1
    sleep_time_scale: float = 3


def retry_on_exception(
    _function: Optional[Callable] = None,
    *,
    config: RetryConfig = RetryConfig(),
    catch: Union[Type[Exception], Tuple[Type[Exception]]] = Exception,
) -> Callable:
    """
    Decorator to retry a function call for a given number of times (while waiting for a certain
    amount between calls), when one of the given exceptions occurs.

    Parameters
    ----------
    config : opencadd.decorator.RetryConfig, optional, default: RetryConfig(3, 1, 3)
        Retry configurations.
    catch : Type[Exception] | Tuple[Type[Exception]], optional, default: Exception
        Exception type(s) that will be ignored during the retries. All other exceptions will be
        raised immediately.

    Returns
    -------
    Callable
        Decorated function.
    """
    if not isinstance(config.num_tries, int) or config.num_tries < 1:
        raise ValueError("`num_tries` must be a positive integer.")

    def retry_decorator(function):
        @wraps(function)
        def retry_wrapper(*args, **kwargs):
            curr_sleep_seconds = config.sleep_time_init
            for try_count in range(config.num_tries):
                try:
                    return function(*args, **kwargs)
                except catch as e:
                    if try_count == config.num_tries - 1:
                        raise e
                    time.sleep(curr_sleep_seconds)
                    curr_sleep_seconds *= config.sleep_time_scale
        return retry_wrapper
    return retry_decorator if _function is None else retry_decorator(_function)
