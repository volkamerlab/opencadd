"""
opencadd.databases.klifs.local.initialize
Utility functions to work with KLIFS data (file)

KLIFS metadata initialization:
Load metadata of KLIFS download by merging details of two KLIFS metadata files, i.e. KLIFS_export.csv and overview.csv.
"""

from pathlib import Path

import pandas as pd


def from_files(klifs_overview_path, klifs_export_path):
    """
    Get KLIFS metadata as DataFrame and save as file in same folder as input files.

    1. Load KLIFS download metadata files.
    2. Unify column names and column cell formatting.
    3. Merge into one DataFrame.

    Parameters
    ----------
    klifs_overview_path : pathlib.Path or str
        Path to KLIFS download file `overview.csv` containing mainly KLIFS alignment-related metadata.
    klifs_export_path : pathlib.Path or str
        Path to KLIFS download file `KLIFS_download/KLIFS_export.csv` containing mainly structure-related metadata.

    Returns
    -------
    pandas.DataFrame
        Metadata of KLIFS download, merged from two KLIFS metadata files.
    """

    klifs_overview_path = Path(klifs_overview_path)
    klifs_export_path = Path(klifs_export_path)

    klifs_overview = _from_klifs_overview_file(klifs_overview_path)
    klifs_export = _from_klifs_export_file(klifs_export_path)

    klifs_metadata = _merge_files(klifs_overview, klifs_export)
    klifs_metadata = _add_filepaths(klifs_metadata)

    klifs_metadata.to_csv(klifs_overview_path.parent / "klifs_metadata.csv", index=False)

    return klifs_metadata


def _from_klifs_export_file(klifs_export_path):
    """
    Read KLIFS_export.csv file from KLIFS database download as DataFrame and unify format with overview.csv format.

    Parameters
    ----------
    klifs_export_path : pathlib.Path or str
        Path to KLIFS_export.csv file from KLIFS database download.

    Returns
    -------
    pandas.DataFrame
        Data loaded and formatted: KLIFS_export.csv file from KLIFS database download.
    """

    klifs_export = pd.read_csv(Path(klifs_export_path))

    # Unify column names with column names in overview.csv
    klifs_export.rename(
        columns={
            "NAME": "kinase",
            "FAMILY": "family",
            "GROUPS": "group",
            "PDB": "pdb_id",
            "CHAIN": "chain",
            "ALTERNATE_MODEL": "alternate_model",
            "SPECIES": "species",
            "LIGAND": "ligand_orthosteric_name",
            "PDB_IDENTIFIER": "ligand_orthosteric_pdb_id",
            "ALLOSTERIC_NAME": "ligand_allosteric_name",
            "ALLOSTERIC_PDB": "ligand_allosteric_pdb_id",
            "DFG": "dfg",
            "AC_HELIX": "ac_helix",
        },
        inplace=True,
    )

    # Unify column 'kinase': Sometimes several kinase names are available, e.g. "EPHA7 (EphA7)"
    # Column "kinase": Retain only first kinase name, e.g. EPHA7
    # Column "kinase_all": Save all kinase names as list, e.g. [EPHA7, EphA7]
    kinase_names = [_format_kinase_name(i) for i in klifs_export.kinase]
    klifs_export.kinase = [i[0] for i in kinase_names]
    klifs_export.insert(loc=1, column="kinase_all", value=kinase_names)

    return klifs_export


def _from_klifs_overview_file(klifs_overview_path):
    """
    Read overview.csv file from KLIFS database download as DataFrame and unify format with KLIFS_export.csv format.

    Parameters
    ----------
    klifs_overview_path : pathlib.Path or str
        Path to overview.csv file from KLIFS database download.

    Returns
    -------
    pandas.DataFrame
        Data loaded and formatted: overview.csv file from KLIFS database download.
    """

    klifs_overview = pd.read_csv(Path(klifs_overview_path))

    # Unify column names with column names in KLIFS_export.csv
    klifs_overview.rename(
        columns={
            "pdb": "pdb_id",
            "alt": "alternate_model",
            "orthosteric_PDB": "ligand_orthosteric_pdb_id",
            "allosteric_PDB": "ligand_allosteric_pdb_id",
        },
        inplace=True,
    )

    # Unify column 'alternate model' with corresponding column in KLIFS_export.csv
    klifs_overview.alternate_model.replace(" ", "-", inplace=True)

    return klifs_overview


def _format_kinase_name(kinase_name):
    """
    Format kinase name(s): One or multiple kinase names (additional names in brackets) are formatted to list of
    kinase names.

    Examples:
    Input: "EPHA7 (EphA7)", output: ["EPHA7", "EphA7"].
    Input: "ITK", output: ["ITK"].

    Parameters
    ----------
    kinase_name : str
        String, here kinase name(s).

    Returns
    -------
    List of str
        List of strings, here list of kinase name(s).
    """

    kinase_name = kinase_name.replace("(", "")
    kinase_name = kinase_name.replace(")", "")
    kinase_name = kinase_name.replace(",", "")
    kinase_name = kinase_name.split()

    return kinase_name


def _merge_files(klifs_export, klifs_overview):
    """
    Merge data contained in overview.csv and KLIFS_export.csv files from KLIFS database download.

    Parameters
    ----------
    klifs_export : pandas.DataFrame
        Metadata contained in KLIFS_export.csv file from KLIFS database download.
    klifs_overview : pandas.DataFrame
        Metadata contained in overview.csv file from KLIFS database download.

    Returns
    -------
    pandas.DataFrame
        Metadata for KLIFS download.
    """

    # Check if PDB IDs occur in one file but not the other
    not_in_export = klifs_export[~klifs_export.pdb_id.isin(klifs_overview.pdb_id)]
    not_in_overview = klifs_overview[~klifs_overview.pdb_id.isin(klifs_export.pdb_id)]

    if not_in_export.size > 0:
        raise ValueError(
            f"Number of PDBs in overview but not in export table: {not_in_export.size}.\n"
        )
    if not_in_overview.size > 0:
        raise (
            f"Number of PDBs in export but not in overview table: {not_in_overview.size}."
            f"PDB codes are probably updated because structures are deprecated."
        )

    # Merge on mutual columns:
    # Species, kinase, PDB ID, chain, alternate model, orthosteric and allosteric ligand PDB ID

    mutual_columns = ["species", "pdb_id", "chain", "alternate_model"]

    klifs_metadata = klifs_export.merge(right=klifs_overview, how="inner", on=mutual_columns)

    klifs_metadata.drop(
        columns=["ligand_orthosteric_pdb_id_y", "ligand_allosteric_pdb_id_y", "kinase_y"],
        inplace=True,
    )

    klifs_metadata.rename(
        columns={
            "ligand_orthosteric_pdb_id_x": "ligand_orthosteric_pdb_id",
            "ligand_allosteric_pdb_id_x": "ligand_allosteric_pdb_id",
            "kinase_x": "kinase",
        },
        inplace=True,
    )

    if not (klifs_overview.shape[1] + klifs_export.shape[1] - 7) == klifs_metadata.shape[1]:
        raise ValueError(
            f"Output table has incorrect number of columns\n"
            f"KLIFS overview table has shape: {klifs_overview.shape}\n"
            f"KLIFS export table has shape: {klifs_export.shape}\n"
            f"KLIFS merged table has shape: {klifs_metadata.shape}"
        )

    if not klifs_overview.shape[0] == klifs_export.shape[0] == klifs_metadata.shape[0]:
        raise ValueError(
            f"Output table has incorrect number of rows:\n"
            f"KLIFS overview table has shape: {klifs_overview.shape}\n"
            f"KLIFS export table has shape: {klifs_export.shape}\n"
            f"KLIFS merged table has shape: {klifs_metadata.shape}"
        )

    return klifs_metadata


def _add_filepaths(klifs_metadata):
    """

    Parameters
    ----------
    klifs_metadata : pandas.DataFrame
        Metadata for KLIFS download.

    Returns
    -------
    pandas.DataFrame
        Metadata for KLIFS download plus column with file paths.
    """

    filepaths = []

    for index, row in klifs_metadata.iterrows():

        # Depending on whether alternate model and chain ID is given build file path:
        mol2_path = Path(".") / row.species.upper() / row.kinase

        if row.alternate_model != "-" and row.chain != "-":
            mol2_path = mol2_path / f"{row.pdb_id}_alt{row.alternate_model}_chain{row.chain}"
        elif row.alternate_model == "-" and row.chain != "-":
            mol2_path = mol2_path / f"{row.pdb_id}_chain{row.chain}"
        elif row.alternate_model == "-" and row.chain == "-":
            mol2_path = mol2_path / f"{row.pdb_id}"
        else:
            raise ValueError(
                f"Incorrect metadata entry {index}: {row.alternate_model}, {row.chain}"
            )

        filepaths.append(mol2_path)

    klifs_metadata["filepath"] = filepaths

    return klifs_metadata
