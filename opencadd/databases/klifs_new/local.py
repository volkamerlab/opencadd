"""
local.py

Defines local KLIFS session.
"""

import logging
from pathlib import Path

from biopandas.mol2 import PandasMol2
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    PocketsProvider,
    CoordinatesProvider,
)
from .schema import LOCAL_COLUMNS_MAPPING, MOL2_COLUMNS, LOCAL_REMOTE_COLUMNS
from .utils import get_file_path

_logger = logging.getLogger(__name__)


class SessionInitializer:
    """
    Class for local session initialization.

    Attributes
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files, including 
        - `overview.csv`, containing mainly KLIFS alignment-related metadata, and
        - `KLIFS_export.csv` containing mainly structure-related metadata.
    klifs_overview_path : pathlib.Path
        Path to `overview.csv` file.
    klifs_export_path : pathlib.Path
        Path to `KLIFS_export.csv` file.
    klifs_metadata : pandas.DataFrame
        Metadata of KLIFS download, merged from two KLIFS metadata files.
    """

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

        klifs_metadata.to_csv(self.path_to_klifs_download / "klifs_metadata.csv", index=False)

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

        # Unify column names with column names in overview.csv
        klifs_export.rename(
            columns=LOCAL_COLUMNS_MAPPING["klifs_export"], inplace=True,
        )

        # Unify column 'kinase.name': Sometimes several kinase names are available, e.g. "EPHA7 (EphA7)"
        # Column "kinase.name": Retain only first kinase name, e.g. EPHA7
        # Column "kinase.name_all": Save all kinase names as list, e.g. [EPHA7, EphA7]
        kinase_names = [self._format_kinase_name(i) for i in klifs_export["kinase.name"]]
        klifs_export["kinase.name"] = [i[0] for i in kinase_names]
        klifs_export.insert(loc=1, column="kinase.name_all", value=kinase_names)

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

        # Unify column names with column names in KLIFS_export.csv
        klifs_overview.rename(
            columns=LOCAL_COLUMNS_MAPPING["klifs_overview"], inplace=True,
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
            "species.klifs",
            "structure.pdb",
            "structure.chain",
            "structure.alternate_model",
        ]

        klifs_metadata = klifs_export.merge(right=klifs_overview, how="inner", on=mutual_columns)

        klifs_metadata.drop(
            columns=["ligand.pdb_y", "ligand.pdb_allosteric_y", "kinase.name_y",], inplace=True,
        )

        klifs_metadata.rename(
            columns={
                "ligand.pdb_x": "ligand.pdb",
                "ligand.pdb_allosteric_x": "ligand.pdb_allosteric",
                "kinase.name_x": "kinase.name",
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

        for index, row in klifs_metadata.iterrows():

            # Depending on whether alternate model and chain ID is given build file path:
            mol2_path = get_file_path(
                ".",
                row["species.klifs"],
                row["kinase.name"],
                row["structure.pdb"],
                row["structure.alternate_model"],
                row["structure.chain"],
                entity="",
                format="",
                in_dir=True,
            )
            filepaths.append(mol2_path)

        klifs_metadata["structure.filepath"] = filepaths

        return klifs_metadata


class Kinases(KinasesProvider):
    """
    Extends KinasesProvider to provide local kinases requests.
    Refer to KinasesProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database

    def all_kinase_groups(self):

        kinase_groups = self.__database["kinase.group"].drop_duplicates()
        kinase_groups = pd.DataFrame(kinase_groups)
        kinase_groups.reset_index(drop=True, inplace=True)

        if kinase_groups.shape[0] > 0:
            return kinase_groups

    def all_kinase_families(self, group=None):

        database = self.__database

        if group:
            database = database[database["kinase.group"] == group]

        kinase_families = database["kinase.family"].drop_duplicates()
        kinase_families = pd.DataFrame(kinase_families)
        kinase_families.reset_index(drop=True, inplace=True)

        if kinase_families.shape[0] > 0:
            return kinase_families
        else:
            raise KeyError(f"None of the input values exist, thus no results are returned.")

    def all_kinases(self, group=None, family=None, species=None):

        # Filter database if filtering parameters are set
        database = self.__database
        if group:
            database = database[database["kinase.group"] == group]
        if family:
            database = database[database["kinase.family"] == family]
        if species:
            database = database[database["species.klifs"] == species.capitalize()]

        # From (filtered) database get unique kinase names
        kinases = database.drop_duplicates("kinase.name")[
            LOCAL_REMOTE_COLUMNS["kinases_all"]["local"]
        ]
        kinases.reset_index(drop=True, inplace=True)

        if kinases.shape[0] > 0:
            return kinases
        else:
            raise KeyError(f"None of the input values exist, thus no results are returned.")

    def from_kinase_names(self, kinase_names, species=None):

        if not isinstance(kinase_names, list):
            kinase_names = [kinase_names]

        # Filter database if filtering parameters are set
        database = self.__database
        database = database[database["kinase.name"].isin(kinase_names)]
        if species:
            database = database[database["species.klifs"] == species]

        # From (filtered) database get unique kinase names
        kinases = database.drop_duplicates("kinase.name")[LOCAL_REMOTE_COLUMNS["kinases"]["local"]]
        kinases.reset_index(drop=True, inplace=True)

        if kinases.shape[0] > 0:
            return kinases
        else:
            raise KeyError(f"None of the input values exist, thus no results are returned.")


class Ligands(LigandsProvider):
    """
    Extends LigandsProvider to provide local ligands requests.
    Refer to LigandsProvider documentation for more information.
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database

    def all_ligands(self):

        ligands = self.__database.drop_duplicates("ligand.pdb")[
            LOCAL_REMOTE_COLUMNS["ligands"]["local"]
        ]
        ligands.reset_index(drop=True, inplace=True)

        if ligands.shape[0] > 0:
            return ligands

    def from_kinase_names(self, kinase_names):

        if not isinstance(kinase_names, list):
            kinase_names = [kinase_names]

        database = self.__database

        # Select ligands by kinase names
        ligands = database[database["kinase.name"].isin(kinase_names)]
        # Keep only columns as per class docstring
        ligands = ligands[
            LOCAL_REMOTE_COLUMNS["ligands"]["local"] + ["kinase.name", "species.klifs"]
        ]

        # Rename columns to indicate columns involved in query
        ligands.rename(
            columns={
                "kinase.name": "kinase.name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )
        ligands.drop_duplicates(inplace=True)
        ligands.reset_index(drop=True, inplace=True)

        if ligands.shape[0] > 0:
            return ligands

    def from_ligand_pdbs(self, ligand_pdbs):

        if not isinstance(ligand_pdbs, list):
            ligand_pdbs = [ligand_pdbs]

        database = self.__database

        # Select ligands by ligand PDB IDs
        ligands = database[database["ligand.pdb"].isin(ligand_pdbs)]
        # Keep only columns as per class docstring
        ligands = ligands[LOCAL_REMOTE_COLUMNS["ligands"]["local"]]

        ligands.drop_duplicates(inplace=True)
        ligands.reset_index(drop=True, inplace=True)

        if ligands.shape[0] > 0:
            return ligands


class Structures(StructuresProvider):
    """
    Extends StructuresProvider to provide local structures requests.
    Refer to StructuresProvider documentation for more information.
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database

    def all_structures(self):

        database = self.__database
        structures = database

        if structures.shape[0] > 0:
            return structures

    def from_structure_pdbs(self, structure_pdbs):

        if not isinstance(structure_pdbs, list):
            structure_pdbs = [structure_pdbs]

        database = self.__database

        # Select structures by structure PDB IDs
        structures = database[database["structure.pdb"].isin(structure_pdbs)]

        if structures.shape[0] > 0:
            return structures

    def from_ligand_pdbs(self, ligand_pdbs):

        if not isinstance(ligand_pdbs, list):
            ligand_pdbs = [ligand_pdbs]

        database = self.__database

        # Select structures by ligand PDB IDs
        structures = database[database["ligand.pdb"].isin(ligand_pdbs)]

        if structures.shape[0] > 0:
            return structures

    def from_kinase_names(self, kinase_names):

        if not isinstance(kinase_names, list):
            kinase_names = [kinase_names]

        database = self.__database

        # Select structures by kinase name (search in all available kinase names)
        structures = database[
            database.apply(
                lambda x: any(
                    [kinase_name in kinase_names for kinase_name in x["kinase.name_all"]]
                ),
                axis=1,
            )
        ]

        if structures.shape[0] > 0:
            return structures


class Bioactivities(BioactivitiesProvider):
    """
    Extends BioactivitiesProvider to provide local bioactivities requests.
    Refer to BioactivitiesProvider documentation for more information.
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database


class Interactions(InteractionsProvider):
    """
    Extends InteractionsProvider to provide local kinases requests.
    Refer to InteractionsProvider documentation for more information.
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database

    def all_interactions(self):

        database = self.__database

        # Keep columns as per class docstring
        interactions = database[
            [
                "structure.pdb",
                "structure.alternate_model",
                "structure.chain",
                "interaction.fingerprint",
            ]
        ]
        if interactions.shape[0] > 0:
            return interactions


class Pockets(PocketsProvider):
    """
    Extends PocketsProvider to provide local pocket requests.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client


class Coordinates(CoordinatesProvider):
    """
    Extends CoordinatesProvider to provide local coordinates requests, 
    i.e. loading structural data (coordinates).
    
    Methods
    -------
    load(file_path, output_format, compute2d)
        Load structural data from KLIFS file in different output formats.
    """

    def __init__(self, database):

        super().__init__()
        self.__database = database

    def load(self, file_path, output_format="biopandas", compute2d=True):
        """
        Load structural data from KLIFS file in different output formats.

        Parameters
        ----------
        file_path : pathlib.Path or str
            Path to KLIFS file.
        entity : str
            Structural entity: complex (default), ligand, pocket, protein, or water.
        compute2d : bool
            For entity=ligand only. Compute 2D coordinates (default) or keep 3D coordinates.
        """

        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File does not exist: {file_path}.")

        # Check if parameters are valid
        entity = file_path.stem
        input_format = file_path.suffix[1:]
        self.check_parameter_validity(entity, input_format, output_format)

        # Return different output formats
        if output_format == "rdkit":
            return self._mol2_file_to_rdkit_mol(str(file_path), compute2d)

        elif output_format == "biopandas":
            if input_format == "mol2":
                return self._mol2_file_to_dataframe(file_path).df
            elif input_format == "pdb":
                pass

    @staticmethod
    def _mol2_file_to_dataframe(mol2_file):
        """
        Get structural data from mol2 file.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
        Path to mol2 file.

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        mol2_file = Path(mol2_file)

        pmol = PandasMol2()

        try:
            mol2_df = pmol.read_mol2(str(mol2_file), columns=MOL2_COLUMNS["n_cols_10"])

        except ValueError:
            mol2_df = pmol.read_mol2(str(mol2_file), columns=MOL2_COLUMNS["n_cols_9"])

        return mol2_df

    @staticmethod
    def _mol2_file_to_rdkit_mol(mol2_file, compute2d=True):
        """
        Get structural data from mol2 file.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
        Path to mol2 file.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromMol2File(mol2_file)

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol
