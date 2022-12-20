from typing import Callable, Optional, Union, Tuple
from functools import wraps
import time


def retry_on_exception(
    _function: Optional[Callable] = None,
    *,
    num_tries: int = 3,
    sleep_seconds: float = 0.1,
    sleep_time_scale: float = 1,
    exception: Union[Exception, Tuple[Exception]] = Exception,
) -> Callable:
    """
    Decorator to recall a function several times, if an exception occurs.

    Parameters
    ----------
    num_tries : int, optional, default: 3
        Maximum number of times the decorated function will be called before the exception is
        reraised.
    sleep_seconds : float, optional, default = 0.1
        Amount of time (in seconds) to wait before two function calls.
    sleep_time_scale : float, optional, default = 1
        Scaling factor for `sleep_seconds`. After the `n`-th function call, the amount of time to
        wait will be equal to `sleep_seconds * (sleep_time_scale ** n-1)`. This can be used to
        scale down/up the waiting time between function calls.
    exception : Exception | Tuple[Exception], optional, default: Exception
        Exception(s) that will be ignored during the retries. All other exceptions will be
        raised immediately.

    Returns
    -------
    Callable
        Decorated function.
    """
    if not isinstance(num_tries, int) or num_tries < 1:
        raise ValueError("`num_tries` must be a positive integer.")

    def retry_decorator(function):
        @wraps(function)
        def retry_wrapper(*args, **kwargs):
            curr_sleep_seconds = sleep_seconds
            for try_count in range(num_tries):
                try:
                    return function(*args, **kwargs)
                except exception as e:
                    if try_count < num_tries - 1:
                        time.sleep(curr_sleep_seconds)
                        curr_sleep_seconds *= sleep_time_scale
                    else:
                        raise e
        return retry_wrapper
    return retry_decorator if _function is None else retry_decorator(_function)
