"""
remote.py

Defines remote KLIFS session.
"""

import logging
from pathlib import Path

from bravado_core.exception import SwaggerMappingError
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    CoordinatesProvider,
)
from .schema import REMOTE_COLUMNS_MAPPING, MOL2_COLUMNS
from .utils import get_file_path


_logger = logging.getLogger(__name__)


class Kinases(KinasesProvider):
    """
    Extends KinasesProvider to provide remote kinases requests.
    Refer to KinasesProvider documentation for more information.

    Attributes
    ----------
    __client : bravado.client.SwaggerClient
        KLIFS client (set if session type is remote).
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    def all_kinase_groups(self):

        try:
            # Use KLIFS API
            result = self.__client.Information.get_kinase_groups().response().result
            # Convert list to DataFrame (1 column)
            kinase_groups = pd.DataFrame(result, columns=["kinase.group"])
            return kinase_groups
        except SwaggerMappingError as e:
            _logger.error(e)

    def all_kinase_families(self, group=None):

        try:
            # Use KLIFS API
            result = (
                self.__client.Information.get_kinase_families(kinase_group=group)
                .response()
                .result
            )
            # Convert list to DataFrame (1 column)
            kinase_families = pd.DataFrame(result, columns=["kinase.family"])
            return kinase_families
        except SwaggerMappingError as e:
            _logger.error(e)

    def all_kinases(self, group=None, family=None, species=None):

        try:
            # Use KLIFS API
            result = (
                self.__client.Information.get_kinase_names(
                    kinase_group=group, kinase_family=family, species=species
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            kinases = self._abc_to_dataframe(result)
            kinases = self._format_dataframe(kinases, REMOTE_COLUMNS_MAPPING["kinases"])
            return kinases
        except SwaggerMappingError as e:
            _logger.error(e)

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        try:
            # Use KLIFS API
            result = (
                self.__client.Information.get_kinase_information(kinase_ID=kinase_ids)
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            kinases = self._abc_to_dataframe(result)
            kinases = self._format_dataframe(kinases, REMOTE_COLUMNS_MAPPING["kinases"])
            return kinases
        except SwaggerMappingError as e:
            _logger.error(e)

    def from_kinase_names(self, kinase_names, species=None):

        kinases = self._iterate_over_function(
            self._from_kinase_name,
            kinase_names,
            additional_parameters=[species],
            iterator_isinstance_check=str,
        )

        if kinases is not None:
            # Kinase IDs can occur multiple times if the input kinase names
            # describe the same kinase, thus drop duplicates
            kinases = kinases.drop_duplicates("kinase.id").reset_index(drop=True)

        return kinases

    def _from_kinase_name(self, kinase_name, species=None):
        """
        Get kinases by kinase name.

        Parameters
        ----------
        kinase_name : str
            Kinase name.

        Returns
        -------
        pandas.DataFrame or None
            Kinases (rows) with columns as described in the class docstring.
        """
        try:
            # Use KLIFS API
            result = (
                self.__client.Information.get_kinase_ID(
                    kinase_name=kinase_name, species=species
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            kinases = self._abc_to_dataframe(result)
            kinases = self._format_dataframe(kinases, REMOTE_COLUMNS_MAPPING["kinases"])
            return kinases
        except SwaggerMappingError as e:
            _logger.error(f"Kinase name {kinase_name}: {e}")


class Ligands(LigandsProvider):
    """
    Extends LigandsProvider to provide remote ligands requests.
    Refer to LigandsProvider documentation for more information.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    def all_ligands(self):

        # Use KLIFS API: Get all kinase IDs
        kinases_remote = Kinases(self.__client)
        kinases = kinases_remote.all_kinases()

        try:
            # Use KLIFS API: Get ligands
            result = (
                self.__client.Ligands.get_ligands_list(
                    kinase_ID=kinases["kinase.id"].to_list()
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            ligands = self._abc_to_dataframe(result)
            ligands = self._format_dataframe(ligands, REMOTE_COLUMNS_MAPPING["ligands"])
            return ligands
        except SwaggerMappingError as e:
            _logger.error(e)

    def from_kinase_ids(self, kinase_ids):

        ligands = self._iterate_over_function(
            self._from_kinase_id, kinase_ids, iterator_isinstance_check=int
        )

        return ligands

    def _from_kinase_id(self, kinase_id):
        """
        Get ligands by kinase ID.

        Parameters
        ----------
        kinase_id : int
            Kinase ID.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        try:
            # Use KLIFS API
            result = (
                self.__client.Ligands.get_ligands_list(kinase_ID=[kinase_id])
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            ligands = self._abc_to_dataframe(result)
            ligands = self._format_dataframe(ligands, REMOTE_COLUMNS_MAPPING["ligands"])
            ligands["kinase.id (query)"] = kinase_id
            return ligands
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Kinase ID {kinase_id}: {e}")

    def from_kinase_names(self, kinase_names):

        if isinstance(kinase_names, str):
            kinase_names = [kinase_names]

        # Use KLIFS API: Get kinase IDs for input kinase names (remotely)
        # Note: One kinase name can be linked to multiple kinase IDs (due to multiple species)
        kinases_remote = Kinases(self.__client)
        kinases = kinases_remote.from_kinase_names(kinase_names)

        # Select and rename columns to indicate columns involved in query
        kinases = kinases[["kinase.id", "kinase.name", "species.klifs"]]
        kinases.rename(
            columns={
                "kinase.id": "kinase.id (query)",
                "kinase.name": "kinase.name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )

        if kinases is not None:

            # Use KLIFS API: Get ligands by kinase IDs
            kinase_ids = kinases["kinase.id (query)"].to_list()
            ligands = self.from_kinase_ids(kinase_ids)

            if ligands is not None:

                # Add kinase name and species details to rationalize kinase IDs
                ligands = ligands.merge(kinases, on="kinase.id (query)", how="left")
                return ligands

    def from_ligand_ids(self, ligand_ids):

        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Use KLIFS API: Get all ligands
        ligands = self.all_ligands()

        # Select ligands by ligand IDs
        if ligands is not None:
            ligands = ligands[ligands["ligand.id"].isin(ligand_ids)]

            if ligands.shape[0] > 0:
                return ligands

    def from_ligand_pdbs(self, ligand_pdbs):

        if isinstance(ligand_pdbs, str):
            ligand_pdbs = [ligand_pdbs]

        # Use KLIFS API: Get all ligands
        ligands = self.all_ligands()

        # Select ligands by ligand PDB IDs
        if ligands is not None:
            ligands = ligands[ligands["ligand.pdb"].isin(ligand_pdbs)]

            if ligands.shape[0] > 0:
                return ligands


class Structures(StructuresProvider):
    """
    Extends StructuresProvider to provide remote structures requests.
    Refer to StructuresProvider documentation for more information.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    def all_structures(self):

        # Use KLIFS API: Get all kinase IDs
        kinases_remote = Kinases(self.__client)
        kinases = kinases_remote.all_kinases()

        # Use KLIFS API: Get all structures from these kinase IDs
        if kinases is not None:
            kinase_ids = kinases["kinase.id"].to_list()
            structures = self.from_kinase_ids(kinase_ids)
            return structures

    def from_structure_ids(self, structure_ids):

        if isinstance(structure_ids, int):
            structure_ids = [structure_ids]

        try:
            # Use KLIFS API
            result = (
                self.__client.Structures.get_structure_list(structure_ID=structure_ids)
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            structures = self._abc_to_dataframe(result)
            structures = self._format_dataframe(
                structures, REMOTE_COLUMNS_MAPPING["structures"]
            )
            return structures
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(e)

    def from_ligand_ids(self, ligand_ids):

        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Use KLIFS API: Get ligand PDB IDs for ligand IDs
        remote_ligands = Ligands(self.__client)
        ligands = remote_ligands.from_ligand_ids(ligand_ids)

        # Use KLIFS API: Get structures from ligand PDBs
        if ligands is not None:
            ligand_pdbs = ligands["ligand.pdb"].to_list()
            structures = self.from_ligand_pdbs(ligand_pdbs)
            return structures

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        try:
            # Use KLIFS API
            result = (
                self.__client.Structures.get_structures_list(kinase_ID=kinase_ids)
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            structures = self._abc_to_dataframe(result)
            structures = self._format_dataframe(
                structures, REMOTE_COLUMNS_MAPPING["structures"]
            )
            return structures
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(e)

    def from_structure_pdbs(self, structure_pdbs):

        if isinstance(structure_pdbs, str):
            structure_pdbs = [structure_pdbs]

        try:
            # Use KLIFS API
            result = (
                self.__client.Structures.get_structures_pdb_list(
                    pdb_codes=structure_pdbs
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            structures = self._abc_to_dataframe(result)
            structures = self._format_dataframe(
                structures, REMOTE_COLUMNS_MAPPING["structures"]
            )
            return structures
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(e)

    def from_ligand_pdbs(self, ligand_pdbs):

        if isinstance(ligand_pdbs, str):
            ligand_pdbs = [ligand_pdbs]

        # Use KLIFS API: Get all structures
        structures = self.all_structures()

        # Select structures by ligand PDB IDs
        if structures is not None:
            structures = structures[structures["ligand.pdb"].isin(ligand_pdbs)]

            if structures.shape[0] > 0:
                return structures

    def from_kinase_names(self, kinase_names):

        if isinstance(kinase_names, str):
            kinase_names = [kinase_names]

        # Use KLIFS API: Get all structures
        structures = self.all_structures()

        # Select structures by kinase names
        if structures is not None:
            structures = structures[structures["kinase.name"].isin(kinase_names)]

            if structures.shape[0] > 0:
                return structures


class Bioactivities(BioactivitiesProvider):
    """
    Extends BioactivitiesProvider to provide remote bioactivities requests.
    Refer to BioactivitiesProvider documentation for more information.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    def all_bioactivities(self):

        # Use KLIFS API: Get all kinase IDs
        ligands_remote = Ligands(self.__client)
        ligands = ligands_remote.all_ligands()

        # Use KLIFS API: Get all bioactivities from these ligand IDs
        if ligands is not None:
            ligand_ids = ligands["ligand.id"].to_list()
            bioactivities = self.from_ligand_ids(ligand_ids)
            return bioactivities

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        # Use KLIFS API: Get all kinase IDs
        ligands_remote = Ligands(self.__client)
        ligands = ligands_remote.from_kinase_ids(kinase_ids)

        # Use KLIFS API: Get all bioactivities from these ligand IDs
        if ligands is not None:
            ligand_ids = ligands["ligand.id"].to_list()
            bioactivities = self.from_ligand_ids(ligand_ids)
            return bioactivities

    def from_ligand_ids(self, ligand_ids):

        bioactivities = self._iterate_over_function(
            self._from_ligand_id, ligand_ids, iterator_isinstance_check=int
        )

        return bioactivities

    def _from_ligand_id(self, ligand_id):
        """
        Get bioactivities by ligand ID.

        Parameters
        ----------
        ligand_id : int
            Ligand ID.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """
        try:
            # Use KLIFS API
            result = (
                self.__client.Ligands.get_bioactivity_list_id(ligand_ID=ligand_id)
                .response()
                .result
            )

            # Convert list of ABC objects to DataFrame and formatting
            bioactivities = self._abc_to_dataframe(result)
            bioactivities = self._format_dataframe(
                bioactivities, REMOTE_COLUMNS_MAPPING["bioactivities"]
            )
            bioactivities["ligand.id (query)"] = ligand_id
            return bioactivities
        except SwaggerMappingError as e:
            _logger.error(f"Ligand ID {ligand_id}: {e}")


class Interactions(InteractionsProvider):
    """
    Extends InteractionsProvider to provide remote kinases requests.
    Refer to InteractionsProvider documentation for more information.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    @property
    def interaction_types(self):

        try:
            # Use KLIFS API
            result = (
                self.__client.Interactions.get_interactions_get_types()
                .response()
                .result
            )

            # Convert list of ABC objects to DataFrame and formatting
            interaction_types = self._abc_to_dataframe(result)
            interaction_types = self._format_dataframe(
                interaction_types, REMOTE_COLUMNS_MAPPING["interactions"]
            )
            return interaction_types
        except SwaggerMappingError as e:
            _logger.error(e)

    def all_interactions(self):

        # Use KLIFS API: Get all structure IDs
        structures_remote = Structures(self.__client)
        structures = structures_remote.all_structures()

        # Use KLIFS API: Get all interactions from these structures IDs
        if structures is not None:
            structure_ids = structures["structure.id"].to_list()
            interactions = self.from_structure_ids(structure_ids)
            return interactions

    def from_structure_ids(self, structure_ids):

        if isinstance(structure_ids, int):
            structure_ids = [structure_ids]

        try:
            # Use KLIFS API
            result = (
                self.__client.Interactions.get_interactions_get_IFP(
                    structure_ID=structure_ids
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame and formatting
            interactions = self._abc_to_dataframe(result)
            interactions = self._format_dataframe(
                interactions, REMOTE_COLUMNS_MAPPING["interactions"]
            )
            return interactions
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(e)

    def from_ligand_ids(self, ligand_ids):

        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Use KLIFS API: Get structure IDs from ligand IDs
        structures_remote = Structures(self.__client)
        structures = structures_remote.from_ligand_ids(ligand_ids)

        # Use KLIFS API: Get interactions from these structure IDs
        if structures is not None:
            structure_ids = structures["structure.id"].to_list()
            interactions = self.from_structure_ids(structure_ids)
            return interactions

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        # Use KLIFS API: Get structure IDs from ligand IDs
        structures_remote = Structures(self.__client)
        structures = structures_remote.from_kinase_ids(kinase_ids)

        # Use KLIFS API: Get interactions from these structure IDs
        if structures is not None:
            structure_ids = structures["structure.id"].to_list()
            interactions = self.from_structure_ids(structure_ids)
            return interactions


class Coordinates(CoordinatesProvider):
    """
    Extends CoordinatesProvider to provide remote coordinates requests, 
    i.e. fetching and saving structural data (coordinates).
    
    Methods
    -------
    fetch(structure_id, entity=, input_format, output_format, compute2d)
        Fetch structural data from KLIFS database in different output formats.
    save(structure_id, output_path, entity, input_format, in_dir)
        Save structural data to file.
    """

    def __init__(self, client):

        super().__init__()
        self.__client = client

    def fetch(
        self,
        structure_id,
        entity="complex",
        input_format="mol2",
        output_format="biopandas",
        compute2d=True,
    ):
        """
        Fetch structural data from KLIFS database in different output formats.

        Parameters
        ----------
        structure_id : str
            KLIFS structure ID.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        input_format : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).
        output_format : str
            Output format: text, biopandas (default), or rdkit (only for entity=ligand).
        compute2d : bool
            For entity=ligand only. Compute 2D coordinates (default) or keep 3D coordinates.
        """

        self.check_parameter_validity(entity, input_format, output_format)

        # Fetch text from KLIFS
        text = self._fetch_text(structure_id, entity, input_format)

        # Return different output formats
        if output_format == "text":
            return text

        elif output_format == "rdkit":
            return self._mol2_text_to_rdkit_mol(text, compute2d)

        elif output_format == "biopandas":
            if input_format == "mol2":
                return self._mol2_text_to_dataframe(text)
            elif input_format == "pdb":
                return self._pdb_text_to_dataframe(text)

    def save(
        self,
        structure_id,
        output_path,
        entity="complex",
        input_format="mol2",
        in_dir=False,
    ):
        """
        Save structural data to file.

        Parameters
        ----------
        structure_id : str
            KLIFS structure ID.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        input_format : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).
        in_dir : bool
            Save file in KLIFS directory structure (default: False).
        """

        self.check_parameter_validity(entity, input_format)
        output_path = Path(output_path)

        # Use KLIFS API: Get structure metadata
        structures_remote = Structures(self.__client)
        metadata = structures_remote.from_structure_ids(structure_id).iloc[0]

        # Set up output path (metadata in the form of directory structure or file name)
        output_path = get_file_path(
            output_path,
            metadata["species.klifs"].upper(),
            metadata["kinase.name"],
            metadata["structure.pdb"],
            metadata["structure.alternate_model"],
            metadata["structure.chain"],
            entity,
            input_format,
            in_dir,
        )

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Get text
        text = self.fetch(structure_id, entity, input_format, output_format="text")

        # Save text to file
        with open(output_path, "w") as f:
            f.write(text)

    def _fetch_text(self, structure_id, entity="complex", input_format="mol2"):
        """
        Get structural data content from KLIFS database as string (text).

        Parameters
        ----------
        structure_id : str
            KLIFS structure ID.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        input_format : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).

        Returns
        -------
        str
            Structural data.
        """

        if entity == "complex" and input_format == "mol2":
            text = (
                self.__client.Structures.get_structure_get_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "complex" and input_format == "pdb":
            text = (
                self.__client.Structures.get_structure_get_pdb_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "ligand" and input_format == "mol2":
            text = (
                self.__client.Structures.get_structure_get_ligand(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "pocket" and input_format == "mol2":
            text = (
                self.__client.Structures.get_structure_get_pocket(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "protein" and input_format == "mol2":
            text = (
                self.__client.Structures.get_structure_get_protein(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        else:
            raise ValueError(f"Entity {entity} is not available remotely.")

        if text is not None:
            return text
        else:
            raise ValueError(f"Data could not be fetched (returned None).")

    @staticmethod
    def _mol2_text_to_dataframe(mol2_text):
        """
        Get structural data from mol2 text.

        Parameters
        ----------
        mol2_text : str
        Mol2 file content from KLIFS database.

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        pmol = PandasMol2()

        try:
            mol2_df = pmol._construct_df(
                mol2_text.splitlines(True),
                col_names=[i[0] for i in MOL2_COLUMNS["n_cols_10"].values()],
                col_types=[i[1] for i in MOL2_COLUMNS["n_cols_10"].values()],
            )
        except ValueError:
            mol2_df = pmol._construct_df(
                mol2_text.splitlines(True),
                col_names=[i[0] for i in MOL2_COLUMNS["n_cols_9"].values()],
                col_types=[i[1] for i in MOL2_COLUMNS["n_cols_9"].values()],
            )

        return mol2_df

    @staticmethod
    def _mol2_text_to_rdkit_mol(mol2_text, compute2d=True):
        """
        Get structural data from mol2 text.

        Parameters
        ----------
        mol2_text : str
            Mol2 file content from KLIFS database.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromMol2Block(mol2_text)

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol

    @staticmethod
    def _pdb_text_to_dataframe(pdb_text):
        """
        Get structural data from pdb text.

        Parameters
        ----------
        pdb_text : str
        Pdb file content from KLIFS database.

        Returns
        -------
        dict of pandas.DataFrame
            Structural data
        """

        ppdb = PandasPdb()

        pdb_dict = ppdb._construct_df(pdb_text.splitlines(True))

        print(f"Structural data keys: {pdb_dict.keys()}")

        return pdb_dict
