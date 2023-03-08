"""
Detect binding pockets using DoGSiteScorer web API.

at DoGSiteScorer functionality of ProteinsPlus server.

References
----------
1.  Volkamer, A.; Griewel, A.; Grombacher, T.; Rarey, M.,
    Analyzing the topology of active sites: on the prediction of pockets and subpockets.
    J. Chem. Inf. Model. 2010, 50 (11), 2041-52. DOI: https://doi.org/10.1021/ci100241y.
2.  Volkamer, A.; Kuhn, D.; Grombacher, T.; Rippmann, F.; Rarey, M.,
    Combining global and local measures for structure-based druggability predictions.
    J. Chem. Inf. Model. 2012, 52 (2), 360-72. DOI: https://doi.org/10.1021/ci200454v.
"""

from pathlib import Path
from typing import Optional, Tuple, Sequence, Union, List, Dict
import io
import gzip

import numpy as np
import pandas as pd
from opencadd._http_request import response_http_request, HTTPRequestRetryConfig
#from opencadd.io.io import filelike_to_filepath
from opencadd._decorator import RetryConfig
import opencadd as oc
import opencadd._regex
import opencadd._common_webapi


class DoGSiteScorerDetector:
    """
    DoGSiteScorer binding site detector.
    """

    def __init__(
            self,
            ensemble,
            pockets,
            pocket_data,
            pocket_atoms,
    ):
        """

        Parameters
        ----------
        ensemble
        pocket_data
        pocket_atoms
        pocket_volumes
        data_description : str
            Description of data descriptors used by DoGSiteScorer.
        """
        self._ensemble = ensemble
        self._pockets = pockets
        self._pocket_data = pocket_data
        self._pocket_atoms = pocket_atoms
        return

    @property
    def ensemble(self):
        return self._ensemble

    @property
    def pockets(self):
        return self._pockets

    @property
    def pocket_data(self):
        """

        Returns
        -------
        pandas.DataFrame
            Columns:

            * name (index) : str
                Name of the pocket, as annotated by DoGSiteScorer.
            * lig_cov : float
                Percentage of ligand volume covered by the pocket.
            * poc_cov : float
                Percentage of pocket volume covered by the ligand.
            * lig_name : str
                PDB name, chain ID, and residue number of the ligand in the PDB file,
                for which `lig_cov` and `poc_cov` are calculated.
            * volume : float
                Pocket volume in cubic angstroms (Å^3), calculated from number of grid points.
            * enclosure : float
                Ratio of number of surface to hull grid points.
            * surface : float
                Pocket surface in square angstroms (Å^2), calculated from number of grid points.
            * depth : float
                Depth of the pocket in angstroms (Å).
            * surf/vol: float
                Surface to volume ratio.
            * lid/hull
            * ellVol
                Ellipsoid volume.
            * ell c/a : float
                Ellipsoid main axis ratio c/a (with a > b > c).
            * ell b/a : float
                Ellipsoid main axis ratio b/a (with a > b > c).
            * siteAtms : int
                Number of surface atoms lining the pocket.
            * accept : int
                Number of hydrogen-bond acceptor atoms.
            * donor : int
                Number of hydrogen-bond donor atoms.
            * hydrophobic_interactions : int
                Number of hydrophobic contacts.
            * hydrophobicity : float
            * metal : int
                Number of metal atoms.
            * Cs : int
                Number of carbon atoms.
            * Ns : int
                Number of nitrogen atoms.
            * Os : int
                Number of oxygen atoms.
            * Ss : int
                Number of sulfur atoms.
            * Xs : int
                Number of other atoms.
            * negAA : float
                Percentage of negatively charged amino acids.
            * posAA : float
                Percentage of positively charged amino acids.
            * polarAA : float
                Percentage of polar amino acids.
            * apolarAA : float
                Percentage of apolar amino acids.
            * ALA : int
            * ARG : int
            * ASN : int
            * ASP : int
            * CYS : int
            * GLN : int
            * GLU : int
            * GLY : int
            * HIS : int
            * ILE : int
            * LEU : int
            * LYS : int
            * MET : int
            * PHE : int
            * PRO : int
            * SER : int
            * THR : int
            * TRP : int
            * TYR : int
            * VAL : int
            * simpleScore : float
                Simple druggability score (in range [0, 1]), based on a linear combination
                of volume, hydrophobicity and enclosure values.
            * drugScore : float
                Druggability score (in range [0, 1]), predicted by a support vector machine (libsvm)
                trained on a subset of meaningful descriptors.
        """
        return self._pocket_data

    @property
    def pocket_atoms(self):
        return self._pocket_atoms

    @property
    def best_coverage(self) -> str:
        """
        Name of the pocket with the highest total coverage for the selected ligand.

        Total coverage is calculated as the product of ligand coverage and pocket coverage.

        Returns
        -------
        str
            Name of the pocket.
        """
        return (self.pocket_data.lig_cov * self.pocket_data.poc_cov).idxmax()


def _from_ensemble(
        ensemble: oc.chem.ensemble.ChemicalEnsemble,
        model: Optional[int] = 0,
        chain_id: Optional[str] = None,
        ligand_id_chain_num: Optional[Tuple[str, str, int]] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> DoGSiteScorerDetector:

    dummy_pdb_id = oc._common_webapi.proteinsplus_upload_pdb(pdb_content=ensemble.to_pdb(model=model))
    pocket_data_df, atoms_df, pocket_volumes, description = _run(
        pdb_id=dummy_pdb_id,
        chain_id=chain_id,
        ligand_id_chain_num=ligand_id_chain_num,
        include_subpockets=include_subpockets,
        calculate_druggability=calculate_druggability,
        retry_config=retry_config
    )
    pockets = {
        pocket_name: oc.pocket.BindingPocket(
            ensemble=ensemble,
            volume=volume,
            atoms=atoms_df.loc[pocket_name, "atom_serial"].to_numpy()
        ) for volume, pocket_name in zip(pocket_volumes, pocket_data_df.index)
    }
    return DoGSiteScorerDetector(ensemble=ensemble, pockets=pockets, pocket_data=pocket_data_df, pocket_atoms=atoms_df)


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
    chain_id : str, optional, default: None
        Chain ID of a polymer instance in the PDB file,
        so that binding pockets are only detected for that chain.
        if not provided (i.e. when set to `None`; default), pockets are detected for all chains.
    ligand_id_chain_num : (str, str, str | int), optional, default: None
        Identifiers for a ligand instance in the PDB file, to calculate its coverage for each pocket.
        The identifiers are: ligand ID (het ID), chain ID, and residue number of the specific instance.
        If not provided (i.e. when set to `None`; default), coverage data will not be calculated.
    include_subpockets : bool, optional, default: True
        Whether to divide detected pockets into sub-pockets and return both (True; default),
        or to only return pockets (False).
    calculate_druggability : bool, optional, default: True
        Whether to calculate druggability scores for each (sub)pocket (True; default).
    retry_config

    Returns
    -------
    DoGSiteScorerDetector
    """
    job_result_urls = _submit_job(
        pdb_id=pdb_id,
        chain_id=chain_id,
        ligand_id_chain_num=ligand_id_chain_num,
        include_subpockets=include_subpockets,
        calculate_druggability=calculate_druggability,
        retry_config=retry_config
    )
    pockets_table, description, ccp4_files, pdb_files = _get_results(
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
        chain_id: Optional[str] = None,
        ligand_id_chain_num: Optional[Tuple[str, str, int]] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> Dict[str, Union[str, List[str]]]:
    # Submit the job and get submission URL
    url_of_job = response_http_request(
        url="https://proteins.plus/api/dogsite_rest",
        verb="POST",
        json={
            "dogsite": {
                "pdbCode": pdb_id,
                "analysisDetail": str(int(include_subpockets)),
                "bindingSitePredictionGranularity": str(int(calculate_druggability)),
                "ligand": "_".join(
                    [str(i) for i in ligand_id_chain_num]) if ligand_id_chain_num is not None else "",
                "chain": chain_id if chain_id is not None else "",
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
        response_verifier=lambda response_dict: "result_table" in response_dict.keys(),
        retry_config=HTTPRequestRetryConfig(
            config_status=retry_config.config_status,
            config_response=RetryConfig(num_tries=30, sleep_time_init=10, sleep_time_scale=1)
        )
    )
    return job_result_urls


def _get_results(
        job_result_urls: dict,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig()
) -> Tuple[str, str, List[bytes], List[bytes]]:
    pockets_table, description = [
        response_http_request(
            url=job_result_urls[data],
            response_type="str",
            retry_config=retry_config
        ) for data in ["result_table", "descriptor_explanation"]
    ]
    ccp4_files, pdb_files = [
        [
            response_http_request(
                url=url_file,
                response_type="bytes",
                retry_config=retry_config
            ) for url_file in url_files
        ] for url_files in (job_result_urls["pockets"], job_result_urls["residues"])
    ]
    return pockets_table, description, ccp4_files, pdb_files


def _parse_results(
        pockets_table: str,
        ccp4_files: Sequence[bytes],
        pdb_files: Sequence[bytes]
):
    """
    Parse the PDB files generated by DoGSiteScorer for each detected pocket.

    Parameters
    ----------
    contents : sequence of bytes
        A sequence of PDB file contents, each as bytes.
    pocket_names : sequence of str
        Name of each pocket, corresponding to `contents`.

    Returns
    -------
    atoms_df, pocket_centers, pocket_radii : pandas.DataFrame, numpy.ndarray, numpy.ndarray
        * atoms_df: Table of atom serials for each pocket.
        * pocket_centers: Coordinates of the geometric center of each pocket.
        * pocket_radii: Maximum radius of each pocket.
    """
    # Create a pandas DataFrame from the pockets' data.
    pocket_data_df = pd.read_csv(io.StringIO(pockets_table), sep="\t").set_index("name")
    # Parse the PDB files; they contain ATOM records from the input PDB file,
    #  corresponding to a specific pocket's residues.
    #  They also contain the geometric center and the max. radius of the pockets.
    pocket_names_repeated, atom_serials, pocket_centers, pocket_radii = [], [], [], []
    for name, content in zip(pocket_data_df.index, pdb_files):
        lines = content.decode().splitlines()
        info_line = lines[5]
        center_and_radius = [float(elem) for elem in info_line.split() if oc._regex.ANY_NUM.match(elem)]
        pocket_centers.append(center_and_radius[:3])
        pocket_radii.append(center_and_radius[3])
        serials = [int(line[6:11]) for line in lines[6:]]
        atom_serials.extend(serials)
        pocket_names_repeated.extend([name] * len(serials))
    # Add centers and radii to the dataframe
    pocket_data_df[["center_x", "center_y", "center_z"]] = np.array(pocket_centers)
    pocket_data_df["max_radius"] = pocket_radii
    # Create another dataframe for atom serial numbers of each pocket
    atoms_df = pd.DataFrame(data={"atom_serial": atom_serials}, index=pocket_names_repeated)
    atoms_df.index.name = "name"
    # Parse CCP4 map files to volumes
    pocket_volumes = [
        oc.io.ccp4.mrc.from_file_content(gzip.decompress(ccp4)) for ccp4 in ccp4_files
    ]
    return pocket_data_df, atoms_df, pocket_volumes
