
from typing import Optional, Tuple, Sequence, Union, List, Dict

from opencadd._http_request import response_http_request, HTTPRequestRetryConfig
from opencadd._decorator import RetryConfig


def _run(
        pdb_id: str,
        chain_id: Optional[str] = None,
        ligand_id_chain_num: Optional[Tuple[str, str, int]] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> Tuple[pd.DataFrame, pd.DataFrame, List, str]:
    """
    Detect binding pockets for a protein structure from the Protein Data Bank.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the protein; either valid 4-letter PDB ID, or dummy ID issued by Proteins.Plus
        for an uploaded PDB file (see `opencadd.webapi.proteins_plus.upload_protein`).
    retry_config

    Returns
    -------
    DoGSiteScorerDetector
    """
    job_result_urls = _submit_job(
        pdb_id=pdb_id,
        retry_config=retry_config
    )
    protein_file, ligand_file, log = _get_results(
        job_result_urls=job_result_urls,
        retry_config=retry_config
    )
    pocket_data_df, atoms_df, pocket_volumes = _parse_results(
        pockets_table=pockets_table,
        ccp4_files=ccp4_files,
        pdb_files=pdb_files
    )
    return pocket_data_df, atoms_df, pocket_volumes, description


def _submit_job(
        pdb_id: str,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> Dict[str, Union[str, List[str]]]:
    # Submit the job and get submission URL
    url_of_job = response_http_request(
        url="https://proteins.plus/api/protoss_rest",
        verb="POST",
        json={
            "protoss": {
                "pdbCode": pdb_id,
            }
        },
        headers={"Content-type": "application/json", "Accept": "application/json"},
        response_type="json",
        response_verifier=lambda response_dict: "location" in response_dict.keys(),
        retry_config=retry_config
    )["location"]
    # Get the URLs of job results when job is done (takes usually 90 seconds)
    job_result_urls: dict = response_http_request(
        url=url_of_job,
        verb="GET",
        response_type="json",
        response_verifier=lambda response_dict: "protein" in response_dict.keys(),
        retry_config=HTTPRequestRetryConfig(
            config_status=retry_config.config_status,
            config_response=RetryConfig(num_tries=30, sleep_time_init=10, sleep_time_scale=1)
        )
    )
    return job_result_urls


def _get_results(
        job_result_urls: dict,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig()
) -> Tuple[bytes, bytes, str]:
    log = response_http_request(
            url=job_result_urls["log"],
            response_type="str",
            retry_config=retry_config
        )
    protein_file, ligand_file = [
        response_http_request(
            url=url_file,
            response_type="bytes",
            retry_config=retry_config
        ) for url_file in (job_result_urls["protein"], job_result_urls["ligands"])
    ]
    return protein_file, ligand_file, log
