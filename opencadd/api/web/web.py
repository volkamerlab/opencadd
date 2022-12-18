"""
Handling HTTP requests and responses.
"""

from typing import Optional, Literal, Union, List, Tuple
import requests


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
    response.raise_for_status()

    if encoding is not None:
        response.encoding = encoding
    if response_type is None:
        return response
    elif response_type == "str":
        return response.text
    elif response_type == "json":
        return response.json(**json_kwargs)
    elif response_type == "bytes":
        return response.content
    raise ValueError(f"`response_type` {response_type} not recognized.")
