"""
opencadd.databases.klifs.local

Defines a local KLIFS session.
"""

import logging
from pathlib import Path

import pandas as pd
from bravado_core.exception import SwaggerMappingError

from . import remote
from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    PocketsProvider,
    CoordinatesProvider,
)
from .schema import (
    FIELDS,
    POCKET_KLIFS_REGIONS,
)
from .utils import PATH_DATA, metadata_to_filepath, filepath_to_metadata
from .remote import KLIFS_CLIENT
from .exceptions import KlifsPocketIncompleteError, KlifsPocketUnequalSequenceStructure
from opencadd.io import DataFrame, Rdkit

# Get the newest file version (* = YYYYMMDD)
PATH_TO_KLIFS_IDS = sorted(PATH_DATA.glob("klifs_ids.*.csv.gz"), key=lambda x: x.suffixes[0])[-1]

_logger = logging.getLogger(__name__)


class LocalInitializer:
    """
    Base class used to define __init__ for all local classes.

    Attributes
    ----------
    _database : pandas.DataFrame
        KLIFS metadata.
    _path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download, *args, **kwargs):

        self._database = database
        self._path_to_klifs_download = path_to_klifs_download


class _LocalDatabaseGenerator:
    """
    Class that generates a local database containing metadata for all downloaded structures.

    Attributes
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files, including
        - `overview.csv`, containing mainly KLIFS alignment-related metadata, and
        - `KLIFS_export.csv` containing mainly structure-related metadata.
    """

    def __init__(self):

        self.path_to_klifs_download = None

    @classmethod
    def from_files(cls, path_to_klifs_download):
        """
        Get KLIFS metadata as DataFrame and save as file in same folder as input files.

        1. Load KLIFS download metadata files.
        2. Unify column names and column cell formatting.
        3. Merge into one DataFrame.

        Parameters
        ----------
        path_to_klifs_download : pathlib.Path or str
            Path to folder with KLIFS download files, including
            - `overview.csv`, containing mainly KLIFS alignment-related metadata, and
            - `KLIFS_export.csv` containing mainly structure-related metadata.

        Returns
        -------
        pandas.DataFrame
            Metadata of KLIFS download, merged from two KLIFS metadata files.
        """

        local_database_generator = cls()

        path_to_klifs_download = Path(path_to_klifs_download)
        klifs_overview_path = path_to_klifs_download / "overview.csv"
        klifs_export_path = path_to_klifs_download / "KLIFS_export.csv"
        if not path_to_klifs_download.exists():
            raise FileNotFoundError(f"No such directory: {path_to_klifs_download}")
        if not klifs_overview_path.exists():
            raise FileNotFoundError(f"No such file: {klifs_overview_path}")
        if not klifs_export_path.exists():
            raise FileNotFoundError(f"No such file: {klifs_export_path}")

        _logger.info(f"Load overview.csv...")
        klifs_overview = local_database_generator._from_klifs_overview_file(klifs_overview_path)
        _logger.info(f"Load KLIFS_export.csv...")
        klifs_export = local_database_generator._from_klifs_export_file(klifs_export_path)

        _logger.info(f"Merge both csv files...")
        database = local_database_generator._merge_files(klifs_overview, klifs_export)
        _logger.info(f"Add paths to coordinate folders to structures...")
        database = local_database_generator._add_filepaths(database)
        _logger.info(f"Add KLIFS IDs to structures (uses remote since not available locally!)...")
        database = local_database_generator._add_klifs_ids(database)

        database.to_csv(path_to_klifs_download / "klifs_metadata.csv", index=False)

        return database

    def _from_klifs_export_file(self, klifs_export_path):
        """
        Read KLIFS_export.csv file from KLIFS database download as DataFrame and unify format with
        overview.csv format.

        Parameters
        ----------
        klifs_export_path : pathlib.Path
            Path to `KLIFS_export.csv` file.

        Returns
        -------
        pandas.DataFrame
            Data loaded and formatted: KLIFS_export.csv file from KLIFS database download.
        """

        klifs_export = pd.read_csv(klifs_export_path)

        # Unify column names with column names in overview.csv
        klifs_export.rename(
            columns=FIELDS.local_export_to_oc_name(),
            inplace=True,
        )

        # Unify column 'kinase.hgnc': Sometimes several kinase names are available, e.g. "EPHA7 (EphA7)"
        # Column "kinase.hgnc": Retain only first kinase name, e.g. EPHA7
        # Column "kinase.all_names": Save all kinase names as list, e.g. [EPHA7, EphA7]
        kinase_names = [
            self._format_kinase_names(i)
            for i in klifs_export["kinase.names"]  # pylint: disable=E1136
        ]
        klifs_export["kinase.names"] = kinase_names
        klifs_export.insert(1, "kinase.gene_name", [i[0] for i in kinase_names])
        klifs_export.insert(2, "kinase.klifs_name", [i[-1] for i in kinase_names])

        return klifs_export

    @staticmethod
    def _from_klifs_overview_file(klifs_overview_path):
        """
        Read overview.csv file from KLIFS database download as DataFrame and unify format with
        KLIFS_export.csv format.

        Parameters
        ----------
        klifs_overview_path : pathlib.Path
            Path to `overview.csv` file.

        Returns
        -------
        pandas.DataFrame
            Data loaded and formatted: overview.csv file from KLIFS database download.
        """

        klifs_overview = pd.read_csv(klifs_overview_path)

        # Unify column names with column names in KLIFS_export.csv
        klifs_overview.rename(
            columns=FIELDS.local_overview_to_oc_name(),
            inplace=True,
        )

        # Unify column 'alternate model' with corresponding column in KLIFS_export.csv
        klifs_overview["structure.alternate_model"].replace(  # pylint: disable=E1136
            " ", "-", inplace=True
        )
        # Drop column for kinase name; will be taken from KLIFS_export.csv file upon table merge
        klifs_overview.drop(columns=["kinase.klifs_name"], inplace=True)

        return klifs_overview

    @staticmethod
    def _format_kinase_names(kinase_names):
        """
        Format kinase name(s): One or multiple kinase names (additional names in brackets) are
        formatted to list of kinase names.

        Examples:
        Input: "EPHA7 (EphA7)", output: ["EPHA7", "EphA7"].
        Input: "ITK", output: ["ITK"].

        Parameters
        ----------
        kinase_names : str
            String, here kinase name(s).

        Returns
        -------
        List of str
            List of strings, here list of kinase name(s).
        """

        kinase_names = kinase_names.replace("(", "")
        kinase_names = kinase_names.replace(")", "")
        kinase_names = kinase_names.replace(",", "")
        kinase_names = kinase_names.split()

        if len(kinase_names) > 2:
            _logger.info(
                f"More than two file names are given: {kinase_names}."
                f"Please file an issue at https://github.com/volkamerlab/opencadd/issues"
            )

        return kinase_names

    @staticmethod
    def _merge_files(klifs_export, klifs_overview):
        """
        Merge data contained in overview.csv and KLIFS_export.csv files from KLIFS database
        download.

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
        not_in_export = klifs_export[
            ~klifs_export["structure.pdb_id"].isin(klifs_overview["structure.pdb_id"])
        ]
        not_in_overview = klifs_overview[
            ~klifs_overview["structure.pdb_id"].isin(klifs_export["structure.pdb_id"])
        ]

        if not_in_export.size > 0:
            raise ValueError(
                f"Number of PDBs in overview but not in export table: {not_in_export.size}.\n"
            )
        if not_in_overview.size > 0:
            raise ValueError(
                f"Number of PDBs in export but not in overview table: {not_in_overview.size}."
                f"PDB codes are probably updated because structures are deprecated."
            )

        # Merge on mutual columns
        mutual_columns = [
            "species.klifs",
            "structure.pdb_id",
            "structure.chain",
            "structure.alternate_model",
            "ligand.expo_id",
            "ligand_allosteric.expo_id",
        ]

        klifs_metadata = klifs_export.merge(right=klifs_overview, how="inner", on=mutual_columns)

        if (
            not (klifs_overview.shape[1] + klifs_export.shape[1] - len(mutual_columns))
            == klifs_metadata.shape[1]
        ):
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

    @staticmethod
    def _add_filepaths(klifs_metadata):
        """
        Add paths to structural data files to DataFrame.

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

        for _, row in klifs_metadata.iterrows():

            # Depending on whether alternate model and chain ID is given build file path:
            filepath = metadata_to_filepath(
                ".",
                row["species.klifs"],
                row["kinase.klifs_name"],
                row["structure.pdb_id"],
                row["structure.alternate_model"],
                row["structure.chain"],
                entity="",
                extension="",
                in_dir=True,
            )
            filepaths.append(filepath)

        klifs_metadata["structure.filepath"] = filepaths

        return klifs_metadata

    @staticmethod
    def _add_klifs_ids(klifs_metadata):
        """
        Add kinase and structure KLIFS IDs to KLIFS metadata (from local copy of KLIFS IDs
        and if not found there from remote). Remove local structures that have no structure ID.
        """

        # Load local copy of KLIFS IDs
        klifs_ids = pd.read_csv(PATH_TO_KLIFS_IDS)
        # Merge KLIFS metadata with local copy of KLIFS IDs
        klifs_metadata_with_ids = klifs_metadata.merge(
            klifs_ids.drop(["kinase.klifs_name", "ligand.expo_id"], axis=1),
            on=["structure.pdb_id", "structure.alternate_model", "structure.chain"],
            how="left",
        )
        if klifs_metadata_with_ids.shape[0] != klifs_metadata.shape[0]:
            raise ValueError(f"Adding KLIFS IDs failed: Number of structures changed.")
        # If KLIFS IDs are missing (not in local copy), fetch individual structures from remote
        remote_structures = remote.Structures(KLIFS_CLIENT)
        for index, row in klifs_metadata_with_ids[
            klifs_metadata_with_ids["structure.klifs_id"].isna()
        ].iterrows():
            # Get IDs from remote
            try:
                structure = remote_structures.by_structure_pdb_id(
                    row["structure.pdb_id"],
                    row["structure.alternate_model"],
                    row["structure.chain"],
                )
                structure_klifs_id = structure["structure.klifs_id"][0]
                kinase_klifs_id = structure["kinase.klifs_id"][0]
                # Set IDs locally
                klifs_metadata_with_ids.loc[index, "structure.klifs_id"] = structure_klifs_id
                klifs_metadata_with_ids.loc[index, "kinase.klifs_id"] = kinase_klifs_id
            except SwaggerMappingError as e:
                # If for some reason some structures are available locally but not remotely
                # we will drop them!
                _logger.error(
                    f"Local structure is not part of {KLIFS_CLIENT} (yet? any more?): "
                    f"{row['structure.pdb_id']}-"
                    f"{row['structure.alternate_model']}-"
                    f"{row['structure.chain']}"
                )
                _logger.error(e)
        # Remove structures that have no KLIFS ID
        klifs_metadata_with_ids.dropna(subset=["structure.klifs_id"], inplace=True)

        return klifs_metadata_with_ids


class Kinases(LocalInitializer, KinasesProvider):
    """
    Extends KinasesProvider to provide local kinases requests.
    Refer to KinasesProvider documentation for more information:
    opencadd.databases.klifs.core.KinasesProvider
    """

    def all_kinase_groups(self):

        # Get local database
        kinase_groups = self._database.copy()
        # Standardize DataFrame
        kinase_groups = self._standardize_dataframe(
            kinase_groups, FIELDS.oc_name_to_type("kinase_groups")
        )
        return kinase_groups

    def all_kinase_families(self, groups=None):

        # Get local database and select rows
        kinase_families = self._database.copy()
        groups = self._ensure_list(groups)
        if groups:
            kinase_families = kinase_families[kinase_families["kinase.group"].isin(groups)]
        # Standardize DataFrame
        kinase_families = self._standardize_dataframe(
            kinase_families, FIELDS.oc_name_to_type("kinase_families")
        )
        return kinase_families

    def all_kinases(self, groups=None, families=None, species=None):

        # Get local database and select rows
        kinases = self._database.copy()
        groups = self._ensure_list(groups)
        families = self._ensure_list(families)
        if groups:
            kinases = kinases[kinases["kinase.group"].isin(groups)]
        if families:
            kinases = kinases[kinases["kinase.family"].isin(families)]
        if species:
            kinases = kinases[kinases["species.klifs"] == species.capitalize()]
        # Standardize DataFrame
        kinases = self._standardize_dataframe(kinases, FIELDS.oc_name_to_type("kinases_all"))
        return kinases

    def by_kinase_klifs_id(self, kinase_klifs_ids):

        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Get local database and select rows
        kinases = self._database.copy()
        kinases = kinases[kinases["kinase.klifs_id"].isin(kinase_klifs_ids)]
        # Standardize DataFrame
        kinases = self._standardize_dataframe(kinases, FIELDS.oc_name_to_type("kinases"))
        return kinases

    def by_kinase_name(self, kinase_names, species=None):

        kinase_names = self._ensure_list(kinase_names)
        # Get local database and select rows
        kinases = self._database.copy()
        # Search in HGNC and KLIFS name columns (case insensitive)
        kinase_names = [kinase_name.upper() for kinase_name in kinase_names]
        kinases = kinases[
            kinases["kinase.klifs_name"].str.upper().isin(kinase_names)
            | kinases["kinase.gene_name"].str.upper().isin(kinase_names)
        ]
        # Search for species (case insensitive)
        if species:
            kinases = kinases[kinases["species.klifs"].str.upper() == species.upper()]
        # Standardize DataFrame
        kinases = self._standardize_dataframe(kinases, FIELDS.oc_name_to_type("kinases"))
        return kinases


class Ligands(LocalInitializer, LigandsProvider):
    """
    Extends LigandsProvider to provide local ligands requests.
    Refer to LigandsProvider documentation for more information:
    opencadd.databases.klifs.core.LigandsProvider
    """

    def all_ligands(self):

        # Get local database
        ligands = self._database.copy()
        # Standardize DataFrame
        ligands = self._standardize_dataframe(ligands, FIELDS.oc_name_to_type("ligands"))
        return ligands

    def by_kinase_klifs_id(self, kinase_klifs_ids):

        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Get local database and select rows
        ligands = self._database.copy()
        ligands = ligands[ligands["kinase.klifs_id"].isin(kinase_klifs_ids)]
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands, FIELDS.oc_name_to_type("ligands", {"kinase.klifs_id": "int32"})
        )
        # Rename columns to indicate columns involved in query TODO remove (query) stuff
        # can columns have metadata?
        # https://github.com/pandas-dev/pandas/issues/2485#issuecomment-608227532
        ligands.rename(
            columns={
                "kinase.klifs_id": "kinase.klifs_id (query)",
            },
            inplace=True,
        )
        return ligands

    def by_kinase_name(self, kinase_names):

        kinase_names = self._ensure_list(kinase_names)
        # Get local database and select rows
        ligands = self._database.copy()
        # Search in HGNC and KLIFS name columns (case insensitive)
        kinase_names = [kinase_name.upper() for kinase_name in kinase_names]
        ligands = ligands[
            ligands["kinase.klifs_name"].str.upper().isin(kinase_names)
            | ligands["kinase.gene_name"].str.upper().isin(kinase_names)
        ]
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands,
            FIELDS.oc_name_to_type(
                "ligands",
                {
                    "kinase.klifs_name": "string",
                    "kinase.gene_name": "string",
                    "species.klifs": "string",
                },
            ),
        )
        # Rename columns to indicate columns involved in query
        ligands.rename(
            columns={
                "kinase.klifs_name": "kinase.klifs_name (query)",
                "kinase.gene_name": "kinase.gene_name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )
        return ligands

    def by_ligand_expo_id(self, ligand_expo_ids):

        ligand_expo_ids = self._ensure_list(ligand_expo_ids)
        # Get local database and select rows
        ligands = self._database.copy()
        ligands = ligands[ligands["ligand.expo_id"].isin(ligand_expo_ids)]
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands,
            FIELDS.oc_name_to_type("ligands"),
        )
        return ligands


class Structures(LocalInitializer, StructuresProvider):
    """
    Extends StructuresProvider to provide local structures requests.
    Refer to StructuresProvider documentation for more information:
    opencadd.databases.klifs.core.StructuresProvider
    """

    def all_structures(self):

        # Get local database
        structures = self._database.copy()
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        return structures

    def by_structure_klifs_id(self, structure_klifs_ids):

        structure_klifs_ids = self._ensure_list(structure_klifs_ids)
        # Get local database and select rows
        structures = self._database.copy()
        structures = structures[structures["structure.klifs_id"].isin(structure_klifs_ids)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        # Check: If only one structure ID was given, only one result is allowed
        if len(structure_klifs_ids) == 1:
            if len(structures) != 1:
                raise ValueError(f"More than one structure found for input structure ID.")

        return structures

    def by_kinase_klifs_id(self, kinase_klifs_ids):

        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Get local database and select rows
        structures = self._database.copy()
        structures = structures[structures["kinase.klifs_id"].isin(kinase_klifs_ids)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        return structures

    def by_structure_pdb_id(
        self, structure_pdb_ids, structure_alternate_model=None, structure_chain=None
    ):

        structure_pdb_ids = self._ensure_list(structure_pdb_ids)
        # Get local database and select rows
        structures = self._database.copy()
        structures = structures[structures["structure.pdb_id"].isin(structure_pdb_ids)]
        # If only one structure PDB ID is given, check alternate model and chain filters
        if len(structure_pdb_ids) == 1:
            structures = self._filter_pdb_by_alt_chain(
                structures, structure_alternate_model, structure_chain
            )
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        return structures

    def by_ligand_expo_id(self, ligand_expo_ids):

        ligand_expo_ids = self._ensure_list(ligand_expo_ids)
        # Get local database and select rows
        structures = self._database.copy()
        structures = structures[structures["ligand.expo_id"].isin(ligand_expo_ids)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        return structures

    def by_kinase_name(self, kinase_names):

        kinase_names = self._ensure_list(kinase_names)
        # Get local database and select rows (search in all available kinase names)
        structures = self._database.copy()
        # Search in HGNC and KLIFS name columns (case insensitive)
        kinase_names = [kinase_name.upper() for kinase_name in kinase_names]
        structures = structures[
            structures["kinase.klifs_name"].str.upper().isin(kinase_names)
            | structures["kinase.gene_name"].str.upper().isin(kinase_names)
        ]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
        )
        return structures


class Bioactivities(LocalInitializer, BioactivitiesProvider):
    """
    Extends BioactivitiesProvider to provide local bioactivities requests.
    Refer to BioactivitiesProvider documentation for more information:
    opencadd.databases.klifs.core.BioactivitiesProvider
    """


class Interactions(LocalInitializer, InteractionsProvider):
    """
    Extends InteractionsProvider to provide local kinases requests.
    Refer to InteractionsProvider documentation for more information:
    opencadd.databases.klifs.core.InteractionsProvider
    """

    def all_interactions(self):

        # Get local database
        interactions = self._database.copy()
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions,
            FIELDS.oc_name_to_type("interactions"),
        )
        return interactions

    def by_structure_klifs_id(self, structure_klifs_ids):

        structure_klifs_ids = self._ensure_list(structure_klifs_ids)
        # Get local database and select rows
        interactions = self._database.copy()
        interactions = interactions[interactions["structure.klifs_id"].isin(structure_klifs_ids)]
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions,
            FIELDS.oc_name_to_type("interactions"),
        )
        return interactions

    def by_kinase_klifs_id(self, kinase_klifs_ids):

        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Get local database and select rows
        interactions = self._database.copy()
        interactions = interactions[interactions["kinase.klifs_id"].isin(kinase_klifs_ids)]
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions, FIELDS.oc_name_to_type("interactions", {"kinase.klifs_id": "int32"})
        )
        # Rename columns to indicate columns involved in query
        interactions.rename(
            columns={
                "kinase.klifs_id": "kinase.klifs_id (query)",
            },
            inplace=True,
        )
        return interactions


class Pockets(LocalInitializer, PocketsProvider):
    """
    Extends PocketsProvider to provide local pocket requests.
    Refer to PocketsProvider documentation for more information:
    opencadd.databases.klifs.core.PocketsProvider
    """

    def by_structure_klifs_id(self, structure_klifs_id, extension="mol2"):  # pylint: disable=W0221

        # Get kinase pocket from structure ID
        structures_local = Structures(self._database, self._path_to_klifs_download)
        structure = structures_local.by_structure_klifs_id(structure_klifs_id).squeeze()

        # CHECK: Number of KLIFS pocket sequence equals 85?
        pocket_sequence = structure["structure.pocket"]
        if len(pocket_sequence) != 85:
            raise KlifsPocketIncompleteError(len(pocket_sequence))

        # Get list of KLIFS positions (starting at 1) excluding gap positions
        klifs_ids = [
            index
            for index, residue in enumerate(structure["structure.pocket"], 1)
            if residue != "_"
        ]

        # Load pocket coordinates from file
        pocket_path = (
            self._path_to_klifs_download / structure["structure.filepath"] / f"pocket.{extension}"
        )
        dataframe = DataFrame.from_file(pocket_path)

        # CHECK: Number of residues in KLIFS pocket sequence and structure file are the same?
        structure_residues = dataframe[["residue.name", "residue.id"]].drop_duplicates()
        if len(klifs_ids) != len(structure_residues):
            raise KlifsPocketUnequalSequenceStructure(len(klifs_ids), len(structure_residues))

        # Get number of atoms per residue
        # Note: sort=False important otherwise negative residue IDs will be sorted to the top
        number_of_atoms_per_residue = dataframe.groupby(
            ["residue.name", "residue.id"], sort=False
        ).size()

        # Get KLIFS position IDs for each atom in molecule
        klifs_ids_per_atom = []
        for klifs_id, n in zip(klifs_ids, number_of_atoms_per_residue):
            klifs_ids_per_atom.extend([klifs_id] * n)
        # Add column for KLIFS position IDs to molecule
        dataframe["residue.klifs_id"] = klifs_ids_per_atom
        dataframe = dataframe[["residue.id", "residue.klifs_id"]].drop_duplicates()

        # Add KLIFS IDs that are missing in pocket and fill with "_"
        full_klifs_ids_df = pd.Series(range(1, 86), name="residue.klifs_id").to_frame()
        dataframe = full_klifs_ids_df.merge(dataframe, on="residue.klifs_id", how="left")
        dataframe.fillna("_", inplace=True)

        # Add column for KLIFS regions
        dataframe = dataframe.merge(POCKET_KLIFS_REGIONS, on="residue.klifs_id", how="left")
        dataframe = dataframe.astype({"residue.klifs_id": "Int64"})

        # Standardize DataFrame
        dataframe = self._standardize_dataframe(
            dataframe,
            FIELDS.oc_name_to_type("pockets"),
        )
        # Add KLIFS region and color  TODO not so nice to have this after standardization
        dataframe = self._add_klifs_region_details(dataframe)

        return dataframe


class Coordinates(LocalInitializer, CoordinatesProvider):
    """
    Extends CoordinatesProvider to provide local coordinates requests,
    i.e. loading structural data (coordinates).
    Refer to CoordinatesProvider documentation for more information:
    opencadd.databases.klifs.core.CoordinatesProvider
    """

    def to_text(
        self, structure_klifs_id_or_filepath, entity="complex", extension="mol2"
    ):  # pylint: disable=W0221

        filepath = self._to_filepath(structure_klifs_id_or_filepath, entity, extension)
        with open(filepath, "r") as f:
            text = f.read()
        return text

    def to_dataframe(
        self, structure_klifs_id_or_filepath, entity="complex", extension="mol2"
    ):  # pylint: disable=W0221

        filepath = self._to_filepath(structure_klifs_id_or_filepath, entity, extension)
        dataframe = DataFrame.from_file(filepath)
        dataframe = self._add_residue_klifs_ids(dataframe, filepath)
        return dataframe

    def to_rdkit(
        self, structure_klifs_id_or_filepath, entity="complex", extension="mol2", compute2d=True
    ):  # pylint: disable=W0221

        filepath = self._to_filepath(structure_klifs_id_or_filepath, entity, extension)
        rdkit_mol = Rdkit.from_file(filepath, compute2d)
        return rdkit_mol

    def _to_filepath(self, structure_klifs_id_or_filepath, entity, extension):
        """
        If the input is a structure ID, return the associated filepath.
        If the input is a filepath already, return that filepath.

        Parameters
        ----------
        structure_klifs_id_or_filepath : int / str or pathlib.Path
            Structure KLIFS ID or path to file.
        entity : str
            Structural entity.
        extension : str
            Input file format.

        Returns
        -------
        pathlib.Path
            Path to file.
        """

        if isinstance(structure_klifs_id_or_filepath, int):
            structure_klifs_id = structure_klifs_id_or_filepath
            # Get structure by structure ID
            structures_local = Structures(self._database, self._path_to_klifs_download)
            structure = structures_local.by_structure_klifs_id(structure_klifs_id).squeeze()
            # Get filepath from metadata
            filepath = metadata_to_filepath(
                self._path_to_klifs_download,
                structure["species.klifs"],
                structure["kinase.klifs_name"],
                structure["structure.pdb_id"],
                structure["structure.alternate_model"],
                structure["structure.chain"],
                entity,
                extension,
                in_dir=True,
            )
        else:
            filepath = structure_klifs_id_or_filepath
        return Path(filepath)

    def _add_residue_klifs_ids(self, dataframe, filepath):
        """
        Add KLIFS position IDs from the KLIFS metadata as additional column.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data.

        Returns
        -------
        pandas.DataFrame
            Structural data including KLIFS residue IDs.
        """

        # Get structure ID from file path
        metadata = filepath_to_metadata(filepath)
        structures_local = Structures(self._database, self._path_to_klifs_download)
        structure = structures_local.by_structure_pdb_id(
            metadata["structure_pdb"],
            metadata["structure_alternate_model"],
            metadata["structure_chain"],
        ).squeeze()
        structure_klifs_id = structure["structure.klifs_id"]

        # Get pocket
        pockets_local = Pockets(self._database, self._path_to_klifs_download)
        pocket_dataframe = pockets_local.by_structure_klifs_id(structure_klifs_id)

        # Merge pocket DataFrame with input DataFrame
        dataframe = dataframe.merge(pocket_dataframe, on="residue.id", how="left")

        return dataframe
