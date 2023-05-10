"""
Common web API calls used by different subpackages/modules within the openCADD library.
"""


from typing import Union
import io

from opencadd._http_request import response_http_request, HTTPRequestRetryConfig


__author__ = "Armin Ariamajd"


def proteinsplus_upload_pdb(
        pdb_content: Union[bytes, str],
        retry_config: HTTPRequestRetryConfig = HTTPRequestRetryConfig(),
) -> str:
    """
    Upload a PDB file to the `ProteinsPlus <https://proteins.plus>`_ webserver.

    ProteinsPlus is a web portal for structural analysis of macromolecules, offering a variety of tools
    for processing structural data, such as protonation/tautomerization state prediction and
    binding pocket detection.

    Proteins.Plus webtools at allow for uploading PBD files
    to use in different functionalities.

    Parameters
    ----------
    pdb_content : str or bytes
        Content of the PDB file.
    retry_config : opencadd._http_request.HTTPRequestRetryConfig, optional
        Retry configurations for HTTP requests.

    Returns
    -------
    str
        A dummy PDB ID, which can be used in place of a real PDB ID in all ProteinsPlus tools.

    References
    ----------
    * `About ProteinsPlus <https://proteins.plus/pages/about>`_
    * `REST API <https://proteins.plus/help/index#REST-help>`_
    """
    # Upload the file; this returns a dict having a "location" key, which holds a web address
    # to retrieve the dummy PDB ID from.
    if isinstance(pdb_content, str):
        pdb_content = pdb_content.encode("utf-8")
    with io.BytesIO(pdb_content) as file:
        url_of_pdb_id = response_http_request(
            url="https://proteins.plus/api/pdb_files_rest",
            verb="POST",
            files={"pdb_file[pathvar]": ("dummy_name.pdb", file)},
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
