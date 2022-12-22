"""
Handling HTTP requests and responses.
"""

# Standard library
from typing import Optional, Literal, Union, List, Tuple, NoReturn, Any, Callable, Sequence, NamedTuple
# 3rd-party
import requests
# Self
from opencadd.decorator import RetryConfig, retry_on_exception


class WebAPIError(IOError):
    """Base Exception class for all web API exceptions."""
    pass


class WebAPIStatusCodeError(WebAPIError):
    """
    Base Exception class for web API status code related exceptions.
    Raised by `opencadd.webapi.http_request.raise_for_status_code`.
    By default, raised when status code is in range [400, 600).
    """
    def __init__(self, response: requests.Response):
        self.response: requests.Response = response
        # Decode error reason from server
        # This part is adapted from `requests` library; See PR #3538 on their GitHub
        if isinstance(response.reason, bytes):
            try:
                reason = response.reason.decode("utf-8")
            except UnicodeDecodeError:
                reason = response.reason.decode("iso-8859-1")
        else:
            reason = response.reason
        error_msg = (
            f"HTTP {'client' if response.status_code < 500 else 'server'} error "
            f"(status code: {response.status_code}) from {response.url}: {reason}\n"
            f"HTTP request: {response.request}\n"
            f"HTTP response: {response.text}"
        )
        super().__init__(error_msg)
        return


class WebAPITemporaryStatusCodeError(WebAPIStatusCodeError):
    """
    Exception class for status code errors related to temporary issues.
    Raised by `opencadd.webapi.http_request.raise_for_status_code`.
    By default, raised when status code is in (408, 429, 500, 502, 503, 504).
    """
    pass


class WebAPIPersistentStatusCodeError(WebAPIStatusCodeError):
    """
    Exception class for status code errors related to persistent issues.
    Raised by `opencadd.webapi.http_request.raise_for_status_code`.
    By default, raised when status code is in range [400, 600),
    but not in (408, 429, 500, 502, 503, 504).
    """
    pass


class WebAPIValueError(WebAPIError):
    """
    Exception class for response value errors.
    Raised by `opencadd.webapi.http_request.response_http_request`.

    """
    def __init__(self, response_value: Any, response_verifier: Callable[[Any], bool]):
        self.response_value = response_value
        self.response_verifier = response_verifier
        error_msg = (
            f"Response verifier function {response_verifier} failed to verify {response_value}."
        )
        super().__init__(error_msg)
        return


def raise_for_status_code(
        response: requests.Response,
        error_status_code_range: Tuple[int, int] = (400, 599),
        temporary_error_status_codes: Optional[Sequence[int]] = (408, 429, 500, 502, 503, 504),
        ignored_status_codes: Optional[Sequence[int]] = None,
) -> NoReturn:
    """
    Raise an `opencadd.webapi.http_request.WebAPIStatusCodeError` for certain HTTP status codes.

    Parameters
    ----------
    response : requests.Response
        Response of an HTTP request, i.e. the returned object from `requests.request`.
    error_status_code_range : Tuple[int, int], optional, default: (400, 599)
        Range of HTTP status codes considered to be error codes.
    temporary_error_status_codes : Sequence[int], optional, default: (408, 429, 500, 502, 503, 504)
        Set of status codes related to temporary errors. If the status code of the response is
        one of these, an `opencadd.webapi.http_request.WebAPITemporaryStatusCodeError` is raised,
        otherwise an `opencadd.webapi.http_request.WebAPIPersistentStatusCodeError` is raised if
        the status code is inside `error_status_code_range`.
    ignored_status_codes : Sequence[int], optional, default: None
        Set of status codes inside `error_status_code_range` to ignore, i.e. not raise.

    Raises
    ------
    opencadd.webapi.http_request.WebAPITemporaryStatusCodeError
        When the status code is in `temporary_error_status_codes`.
    opencadd.webapi.http_request.WebAPIPersistentError
        When the satus code is in range `error_status_code_range`
        and not inside `temporary_error_status_codes`.
    """
    if ignored_status_codes is not None and response.status_code in ignored_status_codes:
        return
    if temporary_error_status_codes is not None and response.status_code in temporary_error_status_codes:
        raise WebAPITemporaryStatusCodeError(response)
    if error_status_code_range[0] <= response.status_code <= error_status_code_range[1]:
        raise WebAPIPersistentStatusCodeError(response)
    return


class HTTPRequestRetryConfig(NamedTuple):
    """
    Configurations for retrying an HTTP request
    sent by `opencadd.webapi.http_request.response_http_request`.

    Attributes
    ----------
    status_codes_to_retry: Sequence[int], optional, default: (408, 429, 500, 502, 503, 504)
        Set of HTTP status codes of response that will trigger a retry.
        If set to `None`, all error status codes will immediately raise an
        `opencadd.webapi.http_request.WebAPIPersistentError`.
    retry_config_status: RetryConfig, optional, default: RetryConfig(5, 1, 2)
        Configurations for retrying when the status code of response is in `status_codes_to_retry`.
        If set to `None`, all error status codes will immediately raise an
        `opencadd.webapi.http_request.WebAPIPersistentError`.
    retry_config_response: RetryConfig, optional, default: RetryConfig(5, 1, 2)
        Configurations for retrying when `response_verifier` returns `False`.
        If set to `None`, all response errors will immediately raise an
        `opencadd.webapi.http_request.WebAPIValueError`.
    """
    status_codes_to_retry: Optional[Sequence[int]] = (408, 429, 500, 502, 503, 504)
    config_status: RetryConfig = RetryConfig(5, 1, 2)
    config_response: RetryConfig = RetryConfig(5, 1, 2)


def response_http_request(
        url: str,
        verb: Union[str, Literal["GET", "POST", "PUT", "PATCH", "OPTIONS", "DELETE"]] = "GET",
        params: Optional[Union[dict, List[tuple], bytes]] = None,
        data: Optional[Union[dict, List[tuple], bytes]] = None,
        headers=None,
        cookies=None,
        files=None,
        auth=None,
        timeout: Optional[Union[float, Tuple[float, float]]] = (10, 20),
        allow_redirects=True,
        proxies=None,
        hooks=None,
        stream=None,
        verify=None,
        cert=None,
        json=None,
        response_type: Optional[Literal["str", "json", "bytes"]] = None,
        encoding: Optional[str] = None,
        response_verifier: Optional[Callable[[Any], bool]] = None,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
        ignored_status_codes: Optional[Sequence[int]] = None,
        **json_kwargs,
) -> Union[requests.Response, str, dict, bytes]:
    """
    Send an HTTP request and get the response in specified type.

    Parameters
    ----------
    url : str
        URL of the API request.
    verb : Literal['GET', 'POST', 'PUT', 'PATCH', 'OPTIONS', 'DELETE'], optional, default: 'GET'
        HTTP verb of the API request. Besides the standard verbs, custom verbs are also accepted.
    params : dict | list[tuple] | bytes, optional, default: None
        Additional parameters to send in the query string of the request.
    data : dict | list[tuple] | bytes, optional, default: None
    response_verifier: Callable[[Any], bool], optional, default: None
        A single-parameter function that when called on the response value of the HTTP request,
        returns a boolean describing whether the response value is accepted (`True`),
        or the request must be retried (`False`). The type of response value is defined by user in
        argument `response_type` of `opencadd.webapi.http_request.response_http_request`.
    headers
    cookies
    files
    auth
    timeout
    allow_redirects
    proxies
    hooks
    stream
    verify
    cert
    json
    response_type
    encoding
    json_kwargs

    Returns
    -------

    Raises
    ------
    requests.exceptions.RequestException
        In case of an unsuccessful request.
    ValueError
        If some arguments are not recognized.

    References
    ----------
    * Documentation of the `requests.request` function:
        https://requests.readthedocs.io/en/latest/api/#requests.request
    * Documentation of JSON decoding keyword arguments:
        https://docs.python.org/3/library/json.html#json.loads
    """
    def get_response_value():

        def get_response():
            response = requests.request(
                method=verb,
                url=url,
                params=params,
                data=data,
                headers=headers,
                cookies=cookies,
                files=files,
                auth=auth,
                timeout=timeout,
                allow_redirects=allow_redirects,
                proxies=proxies,
                hooks=hooks,
                stream=stream,
                verify=verify,
                cert=cert,
                json=json,
            )
            raise_for_status_code(
                response=response,
                temporary_error_status_codes=(
                    None if retry_config is None else retry_config.status_codes_to_retry
                ),
                ignored_status_codes=ignored_status_codes
            )
            return response
        # Depending on specifications in argument `retry_config`, either decorate `get_response`
        # with `retry_on_exception`, or leave it as is.
        response_func = get_response if (
                retry_config is None
                or retry_config.status_codes_to_retry is None
                or retry_config.config_status is None
        ) else retry_on_exception(
                get_response,
                config=retry_config.config_status,
                catch=WebAPITemporaryStatusCodeError
        )
        # Call the (decorated or non-decorated) response function.
        response = response_func()
        # Set encoding of response if specified
        if encoding is not None:
            response.encoding = encoding
        # Get the response value based on specified type
        if response_type is None:
            response_value = response
        elif response_type == "str":
            response_value = response.text
        elif response_type == "json":
            response_value = response.json(**json_kwargs)
        elif response_type == "bytes":
            response_value = response.content
        else:
            raise ValueError(f"`response_type` {response_type} not recognized.")
        # If no verifier is specified, or verifier accepts the response value, the return the value
        if response_verifier is None or response_verifier(response_value):
            return response_value
        # otherwise raise
        raise WebAPIValueError(
            response_value=response_value, response_verifier=response_verifier
        )
    # Depending on specifications in argument `retry_config`, either decorate `get_response_value`
    # with `retry_on_exception`, or leave it as is.
    response_val_func = get_response_value if (
            retry_config is None
            or retry_config.config_response is None
            or response_verifier is None
    ) else retry_on_exception(
            get_response_value,
            config=retry_config.config_response,
            catch=WebAPIValueError,
    )
    # Call the (decorated or non-decorated) response-value function and return.
    return response_val_func()

