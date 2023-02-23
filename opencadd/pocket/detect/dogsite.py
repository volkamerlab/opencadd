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
from typing import Optional, Tuple, Sequence, Union
import io
import gzip

import numpy as np
import pandas as pd
from opencadd._http_request import response_http_request, HTTPRequestRetryConfig
#from opencadd.io.io import filelike_to_filepath
from opencadd._decorator import RetryConfig
import opencadd as oc
import opencadd._regex


class DoGSiteScorerDetector:
    """
    DoGSiteScorer binding site detector.
    """

    def __init__(
            self,
            receptor,
            pocket_data,
            pocket_atoms,
            pocket_volumes,
            data_description: str,
    ):
        """

        Parameters
        ----------
        receptor
        pocket_data
        pocket_atoms
        pocket_volumes
        data_description : str
            Description of data descriptors used by DoGSiteScorer.
        """
        self._receptor = receptor
        self._data = pocket_data
        self._atoms = pocket_atoms
        self._volumes = pocket_volumes
        self._description = data_description
        return
    
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
        return self._data

    def create_pocket(self, pocket_name: str):
        pocket_idx = self._df.index.get_loc(pocket_name)

        toxel_vol = oc.io.ccp4.mrc.from_file_content(gzip.decompress(self._ccp4[pocket_idx]))
        return oc.pocket.BindingPocket(volume=toxel_vol)

    def select_best_pocket(binding_site_df, selection_method, selection_criteria, ascending=False):
        """
        Select the best binding site from the table of all detected binding sites,
        either by sorting the binding sites based on a set of properties in the table,
        or by applying a function on the property values.

        Parameters
        ----------
        binding_site_df : pandas.DataFrame
            Binding site data retrieved from the DoGSiteScorer webserver.
        selection_method : str
            Selection method for selecting the best binding site.
            Either 'sorting' or 'function'.
        selection_criteria : str or list
            If 'selection_method' is 'sorting':
                List of one or more property names.
            If 'selection_method' is 'function':
                Any valid python syntax that generates a list-like object
                with the same length as the number of detected binding sites.
        ascending : bool
            Optional; default: False.
            If set to True, the binding site with the lowest value will be selected,
            otherwise, the binding site with the highest value is selected.

        Returns
        -------
        str
            Name of the selected binding site.
        """
        df = binding_site_df

        if selection_method == "sorting":
            sorted_df = df.sort_values(by=selection_criteria, ascending=ascending)
        elif selection_method == "function":
            df["function_score"] = eval(selection_criteria)
            sorted_df = df.sort_values(by="function_score", ascending=ascending)
        else:
            raise ValueError(f"Binding site selection method unknown: {selection_method}")

        selected_pocket_name = sorted_df.iloc[0].name
        return selected_pocket_name


def from_pdb_id(
        pdb_id: str,
        chain_id: Optional[str] = None,
        ligand_id_chain_num: Optional[Tuple[str, str, int]] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
) -> DoGSiteScorerDetector:
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
    # Submit the job and get submission URL
    url_of_job = response_http_request(
        url="https://proteins.plus/api/dogsite_rest",
        verb="POST",
        json={
            "dogsite": {
                "pdbCode": pdb_id,
                "analysisDetail": str(int(include_subpockets)),
                "bindingSitePredictionGranularity": str(int(calculate_druggability)),
                "ligand": "_".join([str(i) for i in ligand_id_chain_num]) if ligand_id_chain_num is not None else "",
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
    # Parse all PDB files
    atoms_df, pocket_centers, pocket_radii = _parse_dogsite_pdb_files(
        contents=pdb_files, pocket_names=pocket_data_df.index
    )
    # Add centers and radii to the dataframe
    pocket_data_df[["center_x", "center_y", "center_z"]] = pocket_centers
    pocket_data_df["max_radius"] = pocket_radii
    # Parse all CCP4 files
    pocket_volumes = {
        pocket_name: oc.io.ccp4.mrc.from_file_content(gzip.decompress(ccp4))
        for pocket_name, ccp4 in zip(pocket_data_df.index, ccp4_files)
    }
    # Get description file
    data_description = response_http_request(
        url=url_descriptor_descriptions,
        response_type="str",
        retry_config=retry_config
    )
    # Get the parsed PDB file
    receptor_pdb_file = oc.io.pdb.from_pdb_id(pdb_id=pdb_id)
    # Return DoGSiteScorerDetector
    dogsite = DoGSiteScorerDetector(
        receptor=receptor_pdb_file,
        pocket_data=pocket_data_df,
        pocket_atoms=atoms_df,
        pocket_volumes=pocket_volumes,
        data_description=data_description
    )
    return dogsite


def _parse_dogsite_pdb_files(contents: Sequence[bytes], pocket_names: Sequence[str]):
    """
    Parse the PDB files generated by DoGSiteScorer for each detected pocket.

    The PDB files contains ATOM records from the input PDB file,
    corresponding to a specific pocket's residues.
    They also contain the geometric center and the max. radius of the pockets.

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
    pocket_names_repeated, atom_serials, pocket_centers, pocket_radii = [], [], [], []
    for name, content in zip(pocket_names, contents):
        lines = content.decode().splitlines()
        info_line = lines[5]
        center_and_radius = [float(elem) for elem in info_line.split() if oc._regex.ANY_NUM.match(elem)]
        pocket_centers.append(center_and_radius[:3])
        pocket_radii.append(center_and_radius[3])
        serials = [int(line[6:11]) for line in lines[6:]]
        atom_serials.extend(serials)
        pocket_names_repeated.extend([name]*len(serials))
    atoms_df = pd.DataFrame(data={"atom_serial": atom_serials}, index=pocket_names_repeated)
    atoms_df.index.name = "pocket_name"
    pocket_centers = np.array(pocket_centers)
    pocket_radii = np.array(pocket_radii)
    return atoms_df, pocket_centers, pocket_radii
