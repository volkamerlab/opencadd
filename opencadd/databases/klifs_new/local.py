"""
local.py

Defines local KLIFS session.
"""

import logging
from pathlib import Path

import pandas as pd

from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    CoordinatesProvider,
)
from .utils import (
    RENAME_COLUMNS_LOCAL_KLIFS_EXPORT,
    RENAME_COLUMNS_LOCAL_KLIFS_OVERVIEW,
    file_path,
)

_logger = logging.getLogger(__name__)


class SessionInitializer:
    """Class for local session initialization."""

    def __init__(self, path_to_klifs_download):
        self.path_to_klifs_download = Path(path_to_klifs_download)
        self.klifs_overview_path = self.path_to_klifs_download / "overview.csv"
        self.klifs_export_path = self.path_to_klifs_download / "KLIFS_export.csv"
        self.klifs_metadata = self.from_files()

    def from_files(self):
        """
        Get KLIFS metadata as DataFrame and save as file in same folder as input files.

        1. Load KLIFS download metadata files.
        2. Unify column names and column cell formatting.
        3. Merge into one DataFrame.

        Parameters
        ----------
        path_to_klifs_download : pathlib.Path or str
            Path to folder with KLIFS download files `overview.csv`, containing mainly KLIFS alignment-related metadata, and
            `KLIFS_export.csv` containing mainly structure-related metadata.

        Returns
        -------
        pandas.DataFrame
            Metadata of KLIFS download, merged from two KLIFS metadata files.
        """

        klifs_overview = self._from_klifs_overview_file()
        klifs_export = self._from_klifs_export_file()

        klifs_metadata = self._merge_files(klifs_overview, klifs_export)
        klifs_metadata = self._add_filepaths(klifs_metadata)

        klifs_metadata.to_csv(
            self.path_to_klifs_download / "klifs_metadata.csv", index=False
        )

        return klifs_metadata

    def _from_klifs_export_file(self):
        """
        Read KLIFS_export.csv file from KLIFS database download as DataFrame and unify format with overview.csv format.

        Returns
        -------
        pandas.DataFrame
            Data loaded and formatted: KLIFS_export.csv file from KLIFS database download.
        """

        klifs_export = pd.read_csv(self.klifs_export_path)
        print(klifs_export.columns)

        # Unify column names with column names in overview.csv
        klifs_export.rename(
            columns=RENAME_COLUMNS_LOCAL_KLIFS_EXPORT, inplace=True,
        )

        # Unify column 'kinase': Sometimes several kinase names are available, e.g. "EPHA7 (EphA7)"
        # Column "kinase": Retain only first kinase name, e.g. EPHA7
        # Column "kinase_all": Save all kinase names as list, e.g. [EPHA7, EphA7]
        kinase_names = [
            self._format_kinase_name(i) for i in klifs_export["kinase.name"]
        ]
        klifs_export["kinase.name"] = [i[0] for i in kinase_names]
        klifs_export.insert(loc=1, column="kinase_all", value=kinase_names)

        return klifs_export

    def _from_klifs_overview_file(self):
        """
        Read overview.csv file from KLIFS database download as DataFrame and unify format with KLIFS_export.csv format.

        Returns
        -------
        pandas.DataFrame
            Data loaded and formatted: overview.csv file from KLIFS database download.
        """

        klifs_overview = pd.read_csv(self.klifs_overview_path)
        print(klifs_overview.columns)

        # Unify column names with column names in KLIFS_export.csv
        klifs_overview.rename(
            columns=RENAME_COLUMNS_LOCAL_KLIFS_OVERVIEW, inplace=True,
        )

        # Unify column 'alternate model' with corresponding column in KLIFS_export.csv
        klifs_overview["structure.alternate_model"].replace(" ", "-", inplace=True)

        return klifs_overview

    @staticmethod
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

    @staticmethod
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
        not_in_export = klifs_export[
            ~klifs_export["structure.pdb"].isin(klifs_overview["structure.pdb"])
        ]
        not_in_overview = klifs_overview[
            ~klifs_overview["structure.pdb"].isin(klifs_export["structure.pdb"])
        ]

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

        mutual_columns = [
            "species",
            "structure.pdb",
            "structure.chain",
            "structure.alternate_model",
        ]

        klifs_metadata = klifs_export.merge(
            right=klifs_overview, how="inner", on=mutual_columns
        )

        klifs_metadata.drop(
            columns=[
                "ligand.orthosteric.pdb_y",
                "ligand.allosteric.pdb_y",
                "kinase.name_y",
            ],
            inplace=True,
        )

        klifs_metadata.rename(
            columns={
                "ligand.orthosteric.pdb_x": "ligand.orthosteric.pdb",
                "ligand.allosteric.pdb_x": "ligand.allosteric.pdb",
                "kinase.name_x": "kinase.name",
            },
            inplace=True,
        )

        if (
            not (klifs_overview.shape[1] + klifs_export.shape[1] - 7)
            == klifs_metadata.shape[1]
        ):
            raise ValueError(
                f"Output table has incorrect number of columns\n"
                f"KLIFS overview table has shape: {klifs_overview.shape}\n"
                f"KLIFS export table has shape: {klifs_export.shape}\n"
                f"KLIFS merged table has shape: {klifs_metadata.shape}"
            )

        if (
            not klifs_overview.shape[0]
            == klifs_export.shape[0]
            == klifs_metadata.shape[0]
        ):
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
            mol2_path = file_path(
                ".",
                row["species"],
                row["kinase.name"],
                row["structure.pdb"],
                row["structure.alternate_model"],
                row["structure.chain"],
                entity="",
                format="",
                in_dir=True,
            )
            filepaths.append(mol2_path)

        klifs_metadata["filepath"] = filepaths

        return klifs_metadata


class Kinases(KinasesProvider):
    def __init__(self, database):

        super().__init__()
        self.__database = database

    def all_kinase_groups(self):

        kinase_groups = pd.DataFrame(self.__database["kinase.group"].drop_duplicates())
        return kinase_groups

    def all_kinase_families(self, group=None):

        if group:
            try:
                database = self.__database.groupby(["kinase.group"]).get_group(group)
            except KeyError:
                _logger.error(f"Kinase group {group} not known in local dataset.")
                return None
        else:
            database = self.__database

        kinase_families = pd.DataFrame(
            self.__database["kinase.family"].drop_duplicates()
        )
        return kinase_families

    def all_kinases(self, group=None, family=None, species=None):

        kinases = self.__database.drop_duplicates(subset=["kinase.name", "species"])[
            ["kinase.name", "species"]
        ]
        return kinases

    def from_kinases_names(self, kinase_names, species=None):

        self.__database.kinase


class Ligands(LigandsProvider):
    def __init__(self, database):
        super().__init__()
        self.__database = database


class Bioactivities(BioactivitiesProvider):
    def __init__(self, database):
        super().__init__()
        self.__database = database


class Structures(StructuresProvider):
    def __init__(self, database):
        super().__init__()
        self.__database = database


class Interactions(InteractionsProvider):
    def __init__(self, database):
        super().__init__()
        self.__database = database


class Coordinates(CoordinatesProvider):
    def __init__(self, database):
        super().__init__()
        self.__database = database
