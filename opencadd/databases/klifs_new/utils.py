"""
opencadd.databases.klifs.utils

Defines utility functions.
"""

import logging
from pathlib import Path

_logger = logging.getLogger(__name__)

RENAME_COLUMNS_LOCAL = {
    "klifs_export": {
        "NAME": "kinase.name",
        "FAMILY": "kinase.family",
        "GROUPS": "kinase.group",
        "PDB": "structure.pdb",
        "CHAIN": "structure.chain",
        "ALTERNATE_MODEL": "structure.alternate_model",
        "SPECIES": "species.klifs",
        "LIGAND": "ligand.name",
        "PDB_IDENTIFIER": "ligand.pdb",
        "ALLOSTERIC_NAME": "ligand.name_allosteric",
        "ALLOSTERIC_PDB": "ligand.pdb_allosteric",
        "DFG": "structure.dfg",
        "AC_HELIX": "structure.ac_helix",
    },
    "klifs_overview": {
        "species": "species.klifs",
        "kinase": "kinase.name",
        "pdb": "structure.pdb",
        "alt": "structure.alternate_model",
        "chain": "structure.chain",
        "orthosteric_PDB": "ligand.pdb",
        "allosteric_PDB": "ligand.pdb_allosteric",
        "rmsd1": "structure.rmsd1",
        "rmsd2": "structure.rmsd2",
        "qualityscore": "structure.qualityscore",
        "pocket": "kinase.pocket",
        "resolution": "structure.resolution",
        "missing_residues": "structure.missing_residues",
        "missing_atoms": "structure.missing_atoms",
        "full_ifp": "interaction.fingerprint",
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
    },
}

RENAME_COLUMNS_REMOTE = {
    "kinases": {
        "kinase_ID": "kinase.id",
        "name": "kinase.name",
        "HGNC": "kinase.hgnc",
        "family": "kinase.family",
        "group": "kinase.group",
        "kinase_class": "kinase.class",
        "species": "species.klifs",
        "full_name": "kinase.name_full",
        "uniprot": "kinase.uniprot",
        "iuphar": "kinase.iuphar",
        "pocket": "kinase.pocket",
    },
    "structures": {
        "structure_ID": "structure.id",
        "kinase": "kinase.name",
        "species": "species.klifs",
        "kinase_ID": "kinase.id",
        "pdb": "structure.pdb",
        "alt": "structure.alternate_model",
        "chain": "structure.chain",
        "rmsd1": "structure.rmsd1",
        "rmsd2": "structure.rmsd2",
        "pocket": "kinase.pocket",
        "resolution": "structure.resolution",
        "quality_score": "structure.qualityscore",
        "missing_residues": "structure.missing_residues",
        "missing_atoms": "structure.missing_atoms",
        "ligand": "ligand.pdb",
        "allosteric_ligand": "ligand.name_allosteric",
        "DFG": "structure.dfg",
        "aC_helix": "structure.ac_helix",
        "Grich_distance": "structure.grich_distance",
        "Grich_angle": "structure.grich_angle",
        "Grich_rotation": "structure.grich_rotation",
        "front": "structure.front",
        "gate": "structure.gate",
        "back": "structure.back",
        "fp_I": "structure.fp_i",
        "fp_II": "structure.fp_ii",
        "bp_I_A": "structure.bp_i_a",
        "bp_I_B": "structure.bp_i_b",
        "bp_II_in": "structure.bp_ii_in",
        "bp_II_A_in": "structure.bp_ii_a_in",
        "bp_II_B_in": "structure.bp_ii_b_in",
        "bp_II_out": "structure.bp_ii_out",
        "bp_II_B": "structure.bp_ii_b",
        "bp_III": "structure.bp_iii",
        "bp_IV": "structure.bp_iv",
        "bp_V": "structure.bp_v",
        "index": "structure.pocket_klifs_numbering",  # Interactions.get_interactions_match_residues()
        "Xray_position": "structure.pocket_pdb_numbering",  # Interactions.get_interactions_match_residues()
        "KLIFS_position": "structure.pocket_regions_klifs",  # Interactions.get_interactions_match_residues()
    },
    "ligands": {
        "ligand_ID": "ligand.id",
        "PDB-code": "ligand.pdb",
        "Name": "ligand.name",
        "SMILES": "ligand.smiles",
        "InChIKey": "ligand.inchikey",
    },
    "interactions": {
        "position": "interaction.id",  # Interactions.get_interactions_get_types()
        "name": "interaction.name",  # Interactions.get_interactions_get_types()
        "structure_ID": "structure.id",  # Interactions.get_interactions_get_IFP()
        "IFP": "interaction.fingerprint",  # Interactions.get_interactions_get_IFP()
    },
    "bioactivities": {
        "pref_name": "kinase.pref_name",
        "accession": "kinase.uniprot",
        "organism": "species.chembl",
        "standard_type": "ligand.bioactivity_standard_type",
        "standard_relation": "ligand.bioactivity_standard_relation",
        "standard_value": "ligand.bioactivity_standard_value",
        "standard_units": "ligand.bioactivity_standard_units",
        "pchembl_value": "ligand.bioactivity_pchembl_value",
    },
}


def get_file_path(
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

    return path


def _log_error_empty_query_results():
    """
    Log error: Empty results.
    """
    _logger.error(
        "No query results: One or more input parameter values or the combination thereof are invalid or not in the dataset."
    )
