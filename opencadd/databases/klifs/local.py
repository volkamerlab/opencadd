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
    LOCAL_COLUMNS_MAPPING,
    MOL2_COLUMNS,
    LOCAL_REMOTE_COLUMNS,
    POCKET_KLIFS_REGIONS,
)
from .utils import KLIFS_CLIENT, metadata_to_filepath, filepath_to_metadata

PATH_TO_KLIFS_IDS = (
    Path(__file__).parent
    / ".."
    / ".."
    / "data"
    / "klifs_ids.csv"
    # Path(__name__).parent / "opencadd" / "data" / "klifs_ids.csv"
)

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
        if not self.path_to_klifs_download.exists():
            raise FileNotFoundError(f"No such directory: {self.path_to_klifs_download}")

        self.klifs_overview_path = self.path_to_klifs_download / "overview.csv"
        if not self.klifs_overview_path.exists():
            raise FileNotFoundError(f"No such file: {self.klifs_overview_path}")

        self.klifs_export_path = self.path_to_klifs_download / "KLIFS_export.csv"
        if not self.klifs_export_path.exists():
            raise FileNotFoundError(f"No such file: {self.klifs_export_path}")

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

        _logger.info(f"Load overview.csv...")
        klifs_overview = self._from_klifs_overview_file()
        _logger.info(f"Load KLIFS_export.csv...")
        klifs_export = self._from_klifs_export_file()

        _logger.info(f"Merge both csv files...")
        klifs_metadata = self._merge_files(klifs_overview, klifs_export)
        _logger.info(f"Add paths to coordinate folders to structures...")
        klifs_metadata = self._add_filepaths(klifs_metadata)
        _logger.info(f"Add KLIFS IDs to structures (uses remote since not available locally!)...")
        klifs_metadata = self._add_klifs_ids(klifs_metadata)

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
        # Drop column for kinase name; will be taken from KLIFS_export.csv file upon table merge
        klifs_overview.drop(columns=["kinase.name"], inplace=True)

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
            raise ValueError(
                f"Number of PDBs in export but not in overview table: {not_in_overview.size}."
                f"PDB codes are probably updated because structures are deprecated."
            )

        # Merge on mutual columns
        mutual_columns = [
            "species.klifs",
            "structure.pdb",
            "structure.chain",
            "structure.alternate_model",
            "ligand.pdb",
            "ligand.pdb_allosteric",
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

        for index, row in klifs_metadata.iterrows():

            # Depending on whether alternate model and chain ID is given build file path:
            mol2_path = metadata_to_filepath(
                ".",
                row["species.klifs"],
                row["kinase.name"],
                row["structure.pdb"],
                row["structure.alternate_model"],
                row["structure.chain"],
                entity="",
                input_format="",
                in_dir=True,
            )
            filepaths.append(mol2_path)

        klifs_metadata["structure.filepath"] = filepaths

        return klifs_metadata

    @staticmethod
    def _add_klifs_ids(klifs_metadata):
        """
        Add KLIFS kinase and structure IDs to KLIFS metadata (from local copy of KLIFS IDs
        and if not found there from remote). Remove local structures that have no structure ID.
        """

        # Load local copy of KLIFS IDs
        klifs_ids = pd.read_csv(PATH_TO_KLIFS_IDS)
        # Merge KLIFS metadata with local copy of KLIFS IDs
        klifs_metadata_with_ids = klifs_metadata.merge(
            klifs_ids.drop(["kinase.name", "ligand.pdb"], axis=1),
            on=["structure.pdb", "structure.alternate_model", "structure.chain"],
            how="left",
        )
        if klifs_metadata_with_ids.shape[0] != klifs_metadata.shape[0]:
            raise ValueError(f"Adding KLIFS IDs failed: Number of structures changed.")
        # If KLIFS IDs are missing (not in local copy), fetch individual structures from remote
        remote_structures = remote.Structures(KLIFS_CLIENT)
        for index, row in klifs_metadata_with_ids[
            klifs_metadata_with_ids["structure.id"].isna()
        ].iterrows():
            # Get IDs from remote
            structure = remote_structures.from_structure_pdbs(
                row["structure.pdb"], row["structure.alternate_model"], row["structure.chain"],
            )
            structure_id = structure["structure.id"][0]
            kinase_id = structure["kinase.id"][0]
            # Set IDs locally
            klifs_metadata_with_ids.loc[index, "structure.id"] = structure_id
            klifs_metadata_with_ids.loc[index, "kinase.id"] = kinase_id
        # Remove structures that have no KLIFS ID
        klifs_metadata_with_ids.dropna(subset=["structure.id"], inplace=True)

        return klifs_metadata_with_ids


class Kinases(KinasesProvider):
    """
    Extends KinasesProvider to provide local kinases requests.
    Refer to KinasesProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def all_kinase_groups(self):

        # Get local database and select rows
        kinase_groups = self.__database
        # Format DataFrame
        kinase_groups = self._format_dataframe(
            kinase_groups, LOCAL_REMOTE_COLUMNS["kinase_groups"]["local"]
        )
        return kinase_groups

    def all_kinase_families(self, group=None):

        # Get local database and select rows
        kinase_families = self.__database
        if group:
            kinase_families = kinase_families[kinase_families["kinase.group"] == group]
        # Format DataFrame
        kinase_families = self._format_dataframe(
            kinase_families, LOCAL_REMOTE_COLUMNS["kinase_families"]["local"]
        )
        return kinase_families

    def all_kinases(self, group=None, family=None, species=None):

        # Get local database and select rows
        kinases = self.__database
        if group:
            kinases = kinases[kinases["kinase.group"] == group]
        if family:
            kinases = kinases[kinases["kinase.family"] == family]
        if species:
            kinases = kinases[kinases["species.klifs"] == species.capitalize()]
        # Format DataFrame
        kinases = self._format_dataframe(kinases, LOCAL_REMOTE_COLUMNS["kinases_all"]["local"])
        return kinases

    def from_kinase_names(self, kinase_names, species=None):

        kinase_names = self._cast_to_list(kinase_names)
        # Get local database and select rows
        kinases = self.__database
        kinases = kinases[kinases["kinase.name"].isin(kinase_names)]
        if species:
            kinases = kinases[kinases["species.klifs"] == species]
        # Format DataFrame
        kinases = self._format_dataframe(kinases, LOCAL_REMOTE_COLUMNS["kinases"]["local"])
        return kinases


class Ligands(LigandsProvider):
    """
    Extends LigandsProvider to provide local ligands requests.
    Refer to LigandsProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def all_ligands(self):

        # Get local database and select rows
        ligands = self.__database
        # Format DataFrame
        ligands = self._format_dataframe(ligands, LOCAL_REMOTE_COLUMNS["ligands"]["local"])
        return ligands

    def from_kinase_names(self, kinase_names):

        kinase_names = self._cast_to_list(kinase_names)
        # Get local database and select rows
        ligands = self.__database
        ligands = ligands[ligands["kinase.name"].isin(kinase_names)]
        # Format DataFrame
        ligands = self._format_dataframe(
            ligands, LOCAL_REMOTE_COLUMNS["ligands"]["local"] + ["kinase.name", "species.klifs"],
        )
        # Rename columns to indicate columns involved in query
        ligands.rename(
            columns={
                "kinase.name": "kinase.name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )
        return ligands

    def from_ligand_pdbs(self, ligand_pdbs):

        ligand_pdbs = self._cast_to_list(ligand_pdbs)
        # Get local database and select rows
        ligands = self.__database
        ligands = ligands[ligands["ligand.pdb"].isin(ligand_pdbs)]
        # Format DataFrame
        ligands = self._format_dataframe(ligands, LOCAL_REMOTE_COLUMNS["ligands"]["local"],)
        return ligands


class Structures(StructuresProvider):
    """
    Extends StructuresProvider to provide local structures requests.
    Refer to StructuresProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def all_structures(self):

        # Get local database and select rows
        structures = self.__database
        # Format DataFrame
        structures = self._format_dataframe(
            structures, LOCAL_REMOTE_COLUMNS["structures"]["local"],
        )
        return structures

    def from_structure_pdbs(
        self, structure_pdbs, structure_alternate_model=None, structure_chain=None
    ):

        structure_pdbs = self._cast_to_list(structure_pdbs)
        # Get local database and select rows
        structures = self.__database
        structures = structures[structures["structure.pdb"].isin(structure_pdbs)]
        # If only one structure PDB ID is given, check alternate model and chain filters
        if len(structure_pdbs) == 1:
            structures = self._filter_pdb_by_alt_chain(
                structures, structure_alternate_model, structure_chain
            )
        # Format DataFrame
        structures = self._format_dataframe(
            structures, LOCAL_REMOTE_COLUMNS["structures"]["local"],
        )
        return structures

    def from_ligand_pdbs(self, ligand_pdbs):

        ligand_pdbs = self._cast_to_list(ligand_pdbs)
        # Get local database and select rows
        structures = self.__database
        structures = structures[structures["ligand.pdb"].isin(ligand_pdbs)]
        # Format DataFrame
        structures = self._format_dataframe(
            structures, LOCAL_REMOTE_COLUMNS["structures"]["local"],
        )
        return structures

    def from_kinase_names(self, kinase_names):

        kinase_names = self._cast_to_list(kinase_names)
        # Get local database and select rows (search in all available kinase names)
        structures = self.__database
        structures = structures[
            structures.apply(
                lambda x: any(
                    [kinase_name in kinase_names for kinase_name in x["kinase.name_all"]]
                ),
                axis=1,
            )
        ]
        # Format DataFrame
        structures = self._format_dataframe(
            structures, LOCAL_REMOTE_COLUMNS["structures"]["local"],
        )
        return structures

    def from_structure_ids(self, structure_ids):

        structure_ids = self._cast_to_list(structure_ids)
        # Get local database and select rows
        structures = self.__database
        structures = structures[structures["structure.id"].isin(structure_ids)]
        # Format DataFrame
        structures = self._format_dataframe(
            structures, LOCAL_REMOTE_COLUMNS["structures"]["local"],
        )
        # Check: If only one structure ID was given, only one result is allowed
        if len(structure_ids) == 1:
            if len(structures) != 1:
                raise ValueError(f"More than one structure found for input structure ID.")

        return structures


class Bioactivities(BioactivitiesProvider):
    """
    Extends BioactivitiesProvider to provide local bioactivities requests.
    Refer to BioactivitiesProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download


class Interactions(InteractionsProvider):
    """
    Extends InteractionsProvider to provide local kinases requests.
    Refer to InteractionsProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def all_interactions(self):

        # Get local database and select rows
        interactions = self.__database
        # Format DataFrame
        interactions = self._format_dataframe(
            interactions, LOCAL_REMOTE_COLUMNS["interactions"]["local"],
        )
        return interactions

    def from_structure_ids(self, structure_ids):

        structure_ids = self._cast_to_list(structure_ids)
        # Get local database and select rows
        interactions = self.__database
        interactions = interactions[interactions["structure.id"].isin(structure_ids)]
        # Format DataFrame
        interactions = self._format_dataframe(
            interactions, LOCAL_REMOTE_COLUMNS["interactions"]["local"],
        )
        return interactions

    def from_kinase_ids(self, kinase_ids):

        kinase_ids = self._cast_to_list(kinase_ids)
        # Get local database and select rows
        interactions = self.__database
        interactions = interactions[interactions["kinase.id"].isin(kinase_ids)]
        # Format DataFrame
        interactions = self._format_dataframe(
            interactions, LOCAL_REMOTE_COLUMNS["interactions"]["local"],
        )
        return interactions


class Pockets(PocketsProvider):
    """
    Extends PocketsProvider to provide local pocket requests.
    Refer to PocketsProvider documentation for more information.

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def from_structure_id(self, structure_id):

        # Get pocket coordinates
        coordinates_local = Coordinates(self.__database, self.__path_to_klifs_download)
        mol2_df = coordinates_local.from_structure_id(
            structure_id, entity="pocket", input_format="mol2", output_format="biopandas",
        )
        # Format DataFrame
        mol2_df = self._format_dataframe(mol2_df, LOCAL_REMOTE_COLUMNS["pockets"]["local"],)

        return mol2_df


class Coordinates(CoordinatesProvider):
    """
    Extends CoordinatesProvider to provide local coordinates requests, 
    i.e. loading structural data (coordinates).

    Attributes
    ---------- 
    __database : pandas.DataFrame
        KLIFS metadata (set if session type is local).
    __path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    """

    def __init__(self, database, path_to_klifs_download):

        super().__init__()
        self.__database = database
        self.__path_to_klifs_download = path_to_klifs_download

    def from_file(self, file_path, output_format="biopandas", compute2d=True):
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
            rdkit_mol = self._mol2_file_to_rdkit_mol(str(file_path), compute2d)
            return rdkit_mol

        elif output_format == "biopandas":
            if input_format == "mol2":
                mol2_df = self._mol2_file_to_dataframe(file_path).df
                mol2_df = self._split_mol2_subst_names(mol2_df)
                if "pocket" in str(file_path):
                    mol2_df = self._add_residue_klifs_ids(mol2_df, file_path)
                return mol2_df
            elif input_format == "pdb":
                pass

    def from_structure_id(
        self,
        structure_id,
        entity="complex",
        input_format="mol2",
        output_format="biopandas",
        compute2d=True,
    ):
        self.check_parameter_validity(entity, input_format, output_format)

        # Get structure by structure ID
        structures_local = Structures(self.__database, self.__path_to_klifs_download)
        structure = structures_local.from_structure_ids(structure_id).squeeze()
        # Get filepath from metadata
        filepath = metadata_to_filepath(
            self.__path_to_klifs_download,
            structure["species.klifs"],
            structure["kinase.name"],
            structure["structure.pdb"],
            structure["structure.alternate_model"],
            structure["structure.chain"],
            entity,
            input_format,
            in_dir=True,
        )
        # Get coordinates from file
        result = self.from_file(filepath, output_format, compute2d)
        return result

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

    def _add_residue_klifs_ids(self, mol2_df, filepath):
        """
        Add KLIFS position IDs from the KLIFS metadata as additional column.

        Parameters
        ---------- TODO

        Returns
        ------- TODO
        """

        # Get structure ID from file path
        metadata = filepath_to_metadata(filepath)
        structures_local = Structures(self.__database, self.__path_to_klifs_download)
        structure = structures_local.from_structure_pdbs(
            metadata["structure_pdb"],
            metadata["structure_alternate_model"],
            metadata["structure_chain"],
        ).squeeze()
        structure_id = structure["structure.id"]
        # List of KLIFS positions (starting at 1) excluding gap positions
        klifs_ids = [
            index for index, residue in enumerate(structure["kinase.pocket"], 1) if residue != "_"
        ]

        # Number of atoms per residue in molecule (mol2file)
        # Note: sort=False important otherwise negative residue IDs will be sorted to the top
        number_of_atoms_per_residue = mol2_df.groupby(by="residue.subst_name", sort=False).size()

        # Get KLIFS position IDs for each atom in molecule
        klifs_ids_per_atom = []
        for klifs_id, n in zip(klifs_ids, number_of_atoms_per_residue):
            klifs_ids_per_atom.extend([klifs_id] * n)
        # Add column for KLIFS position IDs to molecule
        mol2_df["residue.klifs_id"] = klifs_ids_per_atom
        # Add column for KLIFS regions
        pocket_klifs_regions = (
            pd.Series(POCKET_KLIFS_REGIONS, name="residue.klifs_region")
            .reset_index()
            .rename(columns={"index": "residue.klifs_id"})
        )
        mol2_df = mol2_df.merge(pocket_klifs_regions, on="residue.klifs_id", how="left")

        return mol2_df
