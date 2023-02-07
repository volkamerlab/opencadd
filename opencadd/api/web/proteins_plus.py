"""
REST-API implementation of ProteinsPlus tools at https://proteins.plus/

References
----------
Tools:
    * https://proteins.plus/pages/about
REST-API:
    * https://proteins.plus/help/index#REST-help
    * https://proteins.plus/help/dogsite_rest
"""


from typing import Optional, Tuple
import io

import pandas as pd

from opencadd._typing import FileLike
from opencadd._http_request import response_http_request, HTTPRequestRetryConfig
#from opencadd.io.io import filelike_to_filepath
from opencadd._decorator import RetryConfig


_URL_ROOT: str = "https://proteins.plus/api"
_URL_FILE_UPLOAD: str = f"{_URL_ROOT}/pdb_files_rest"
_URL_DOGSITESCORER: str = f"{_URL_ROOT}/dogsite_rest"


def upload_protein(
        pdb_file: FileLike,
        retry_config: HTTPRequestRetryConfig = HTTPRequestRetryConfig(),
) -> str:
    """
    Upload a PDB file to Proteins.Plus server to use in their tools.

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
    """
    # Upload the file; this returns a dict having a "location" key, which holds a web address
    # to retrieve the dummy PDB ID from.
    with open(filelike_to_filepath(pdb_file), "rb") as file:
        url_of_pdb_id = response_http_request(
            url=_URL_FILE_UPLOAD,
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


def dog_site_scorer(
        pdb_id: str,
        chain_id: Optional[str] = None,
        ligand_id: Optional[str] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> Tuple[pd.DataFrame, Tuple[bytes], Tuple[bytes], str]:
    """
    Detect binding pockets in a protein structure using DoGSiteScorer functionality of
    ProteinsPlus server.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the protein; either valid 4-letter PDB ID, or dummy ID issued by Proteins.Plus
        for an uploaded PDB file (see `opencadd.webapi.proteins_plus.upload_protein`).
    chain_id : str, optional, default: None
        ID of a protein chain present in the PDB file. if provided, binding sites are only
        detected on this chain.
    ligand_id : str, optional, default: None
        ID of a ligand present in the PDB file. If provided, ligand coverage will also be
        calculated for each detected (sub)pocket. Ligand ID must be in DoGSiteScorer format:
        <PDB ligand ID>_<chain ID>_<PDB residue number of the ligand>
        Example: For PDB-ID 3W32, the co-crystallized ligand (ID: W32), which is located on
        chain A and has residue number 1101, the DoGSiteScorer ligand ID is: W32_A_1101
    include_subpockets : bool, optional, default: True
        Whether to divide detected pockets into sub-pockets and return both, or to only return
        pockets.
    calculate_druggability : bool, optional, default: True
        Whether to calculate druggability scores for each (sub)pocket.
    retry_config

    Returns
    -------
    pockets_data, pockets_map_files, pockets_pdb_files, data_description : Tuple
        pockets_data : pandas.DataFrame
            All calculated data for detected binding pockets (& sub-pockets).
        pockets_map_files : Tuple[bytes]
            Binary file contents of pocket map files (one for each detected sub(pocket))
            in CCP4.GZ format (i.e. Gzipped CCP4).
        pockets_pdb_files : Tuple[bytes]
            Binary file contents of PDB files (one for each detected sub(pocket)), containing the
            surrounding atoms of each pocket, plus the position of pocket's center and its radius.
        data_description : str
            Description of each data descriptor in `pockets_data`.

    References
    ----------
    Publications:
        1.  Volkamer, A.; Griewel, A.; Grombacher, T.; Rarey, M.,
            Analyzing the topology of active sites: on the prediction of pockets and subpockets.
            J. Chem. Inf. Model. 2010, 50 (11), 2041-52. DOI: https://doi.org/10.1021/ci100241y.
        2.  Volkamer, A.; Kuhn, D.; Grombacher, T.; Rippmann, F.; Rarey, M.,
            Combining global and local measures for structure-based druggability predictions.
            J. Chem. Inf. Model. 2012, 52 (2), 360-72. DOI: https://doi.org/10.1021/ci200454v.

    """
    # Submit the job and get submission URL
    url_of_job = response_http_request(
        url=_URL_DOGSITESCORER,
        verb="POST",
        json={
            "dogsite": {
                "pdbCode": pdb_id,
                "analysisDetail": str(int(include_subpockets)),
                "bindingSitePredictionGranularity": str(int(calculate_druggability)),
                "ligand": ligand_id if ligand_id is not None else "",
                "chain": chain_id if chain_id is not None else "",
            }
        },
        headers={"Content-type": "application/json", "Accept": "application/json"},
        response_type="json",
        response_verifier=lambda response_dict: "location" in response_dict.keys(),
        retry_config=retry_config
    )["location"]
    # Get the URLs of job results when job is done (takes usually 90 seconds)
    job_result = response_http_request(
        url=url_of_job,
        verb="GET",
        response_type="json",
        response_verifier=lambda response_dict: "result_table" in response_dict.keys(),
        retry_config=HTTPRequestRetryConfig(
            config_status=retry_config.config_status,
            config_response=RetryConfig(num_tries=5, sleep_time_init=55, sleep_time_scale=0.5)
        )
    )
    url_data = job_result["result_table"]
    url_ccp4_files = job_result["pockets"]
    url_pdb_files = job_result["residues"]
    url_descriptor_descriptions = job_result["descriptor_explanation"]
    # Get results data from URL
    pocket_data = response_http_request(
        url=url_data,
        response_type="str",
        retry_config=retry_config
    )
    pocket_data_df = pd.read_csv(io.StringIO(pocket_data), sep="\t").set_index("name")
    # Get files from URLs
    ccp4_files, pdb_files = [
        [
            response_http_request(
                url=url_file,
                response_type="bytes",
                retry_config=retry_config
            ) for url_file in url_files
        ] for url_files in (url_ccp4_files, url_pdb_files)
    ]
    # Get description file
    description_file = response_http_request(
        url=url_descriptor_descriptions,
        response_type="str",
        retry_config=retry_config
    )
    return pocket_data_df, tuple(ccp4_files), tuple(pdb_files), description_file







