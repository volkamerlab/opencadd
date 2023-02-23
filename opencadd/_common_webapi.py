"""
REST-API implementation of ProteinsPlus tools

References
----------
Tools:
    * https://proteins.plus/pages/about
REST-API:
    * https://proteins.plus/help/dogsite_rest
"""


from typing import Optional, Tuple
import io

import pandas as pd

from opencadd._typing import FileLike
from opencadd._http_request import response_http_request, HTTPRequestRetryConfig
#from opencadd.io.io import filelike_to_filepath
from opencadd._decorator import RetryConfig




def upload_protein(
        pdb_file: FileLike,
        retry_config: HTTPRequestRetryConfig = HTTPRequestRetryConfig(),
) -> str:
    """
    Upload a PDB file to Proteins.Plus webserver.

    Proteins.Plus webtools at https://proteins.plus allow for uploading PBD files
    to use in different functionalities.
    After uploading, a dummy PDB ID is returned, which can be used in place of a valid PDB ID
    in all available tools.

    Parameters
    ----------
    pdb_file : FileLike
        A file-like object in PDB format.
    retry_config : opencadd.webapi.http_request.HTTPRequestRetryConfig, optional,
                   default: HTTPRequestRetryConfig()
        Retry configurations for HTTP requests.

    Returns
    -------
    dummy_pdb_id : str
        A dummy PDB ID, which can be used in place of a real PDB ID in Proteins.Plus tools.

    References
    ----------
    https://proteins.plus/help/index#REST-help
    """
    # Upload the file; this returns a dict having a "location" key, which holds a web address
    # to retrieve the dummy PDB ID from.
    with open(filelike_to_filepath(pdb_file), "rb") as file:
        url_of_pdb_id = response_http_request(
            url="https://proteins.plus/api/pdb_files_rest",
            verb="POST",
            files={"pdb_file[pathvar]": file},
            response_type="json",
            response_verifier=lambda response_dict: "location" in response_dict.keys(),
            retry_config=retry_config
        )["location"]
    # Try to retrieve the dummy PDB ID from the address.
    dummy_pdb_id = response_http_request(
        url=url_of_pdb_id,
        verb="GET",
        response_type="json",
        response_verifier=lambda response_dict: "id" in response_dict.keys(),
        retry_config=retry_config,
    )["id"]
    return dummy_pdb_id










