"""
opencadd.databases.klifs.utils

Defines utility functions.
"""

from pathlib import Path

import pandas as pd

RENAME_COLUMNS_LOCAL_KLIFS_EXPORT = {
    "NAME": "kinase.name",
    "FAMILY": "kinase.family",
    "GROUPS": "kinase.group",
    "PDB": "structure.pdb",
    "CHAIN": "structure.chain",
    "ALTERNATE_MODEL": "structure.alternate_model",
    "SPECIES": "species",
    "LIGAND": "ligand.orthosteric.name",
    "PDB_IDENTIFIER": "ligand.orthosteric.pdb",
    "ALLOSTERIC_NAME": "ligand.allosteric.name",
    "ALLOSTERIC_PDB": "ligand.allosteric.pdb",
    "DFG": "structure.dfg",
    "AC_HELIX": "structure.ac_helix",
}
RENAME_COLUMNS_LOCAL_KLIFS_OVERVIEW = {
    "species": "species",
    "kinase": "kinase.name",
    "pdb": "structure.pdb",
    "alt": "structure.alternate_model",
    "chain": "structure.chain",
    "orthosteric_PDB": "ligand.orthosteric.pdb",
    "allosteric_PDB": "ligand.allosteric.pdb",
    "rmsd1": "structure.rmsd1",
    "rmsd2": "structure.rmsd2",
    "qualityscore": "structure.qualityscore",
    "pocket": "kinase.pocket",
    "resolution": "structure.resolution",
    "missing_residues": "structure.missing_residues",
    "missing_atoms": "structure.missing_atoms",
    "full_ifp": "structure.ifp",
    "fp_i": "structure.fp_i",
    "fp_ii": "structure.fp_ii",
    "bp_i_a": "structure.bp_i_a",
    "bp_i_b": "structure.bp_i_b",
    "bp_ii_in": "structure.bp_ii_in",
    "bp_ii_a_in": "structure.bp_ii_a_in",
    "bp_ii_b_in": "structure.bp_ii_b_in",
    "bp_ii_out": "structure.bp_ii_out",
    "bp_ii_b": "structure.bp_ii_b",
    "bp_iii": "structure.bp_iii",
    "bp_iv": "structure.bp_iv",
    "bp_v": "structure.bp_v",
}
RENAME_COLUMNS_REMOTE_KINASE = {
    "kinase_ID": "kinase.id",
    "name": "kinase.name",
    "HGNC": "kinase.hgnc",
    "family": "kinase.family",
    "group": "kinase.group",
    "kinase_class": "kinase.class",
    "species": "species",
    "full_name": "kinase.name_full",
    "uniprot": "kinase.uniprot",
    "iuphar": "kinase.iuphar",
    "pocket": "kinase.pocket",
}


def file_path(
    path_to_klifs_download,
    species,
    kinase_name,
    structure_pdb,
    structure_alternate_model,
    structure_chain,
    entity="complex",
    format="mol2",
    in_dir=False,
):
    """
    Get file path.

    Parameters
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder for file destination (if in_dir=False) or KLIFS_downlaod folder destination (if in_dir=True).
    species : str
        Species.
    kinase_name : str
        Kinase name.
    structure_pdb : str
        PDB ID.
    structure_alternate_model : str
        Alternate model ID.
    structure_chain : str
        Chain ID.
    entitiy : str
        Structural entity: complex (default), ligand, pocket, protein, or water (only in local module).
    format : str
        File format: mol2 (default) or pdb (only for entity=complex).
    in_dir : bool
        Use KLIFS directory structure (default: False).

    Returns
    -------
    pathlib.Path
        File path.
    """

    path_to_klifs_download = Path(path_to_klifs_download)
    species = species.upper()
    structure_alternate_model = structure_alternate_model.replace("-", "")
    structure_chain = structure_chain.replace("-", "")
    structure = f"{structure_pdb}{f'_alt{structure_alternate_model}' if bool(structure_alternate_model) else ''}{f'_chain{structure_chain}' if bool(structure_chain) else ''}"

    if in_dir:
        path = (
            path_to_klifs_download
            / "KLIFS_download"
            / species
            / kinase_name
            / structure
            / f"{entity}.{format}"
        )
    else:
        path = (
            path_to_klifs_download
            / f"{species}_{kinase_name}_{structure}_{entity}.{format}"
        )

    print(path)

    return path


def _abc_idlist_to_dataframe(abc_idlist):
    """
    Transform ABC IDList object into DataFrame.

    Parameters
    ----------
    abc_idlist : list of acb.IDList
        List of labeled list objects from abstract base classes module.

    Returns
    -------
    pandas.DataFrame
        Table with list labels as column names.
    """

    result = abc_idlist

    keys = list(result[0])

    results_dict = {key: [] for key in keys}

    for result in abc_idlist:
        for key in keys:
            results_dict[key].append(result[key])

    return pd.DataFrame(results_dict)
