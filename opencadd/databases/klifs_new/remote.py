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
from .schema import RENAME_COLUMNS_REMOTE
from .utils import get_file_path, _log_error_empty_query_results


_logger = logging.getLogger(__name__)


class Kinases(KinasesProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_kinase_groups(self):

        results = self.__client.Information.get_kinase_groups().response().result
        kinase_groups = pd.DataFrame(results, columns=["kinase.group"])
        kinase_groups.rename(columns=RENAME_COLUMNS_REMOTE["kinases"], inplace=True)
        return kinase_groups

    def all_kinase_families(self, group=None):

        try:
            results = (
                self.__client.Information.get_kinase_families(kinase_group=group)
                .response()
                .result
            )
            results = pd.DataFrame(results, columns=["kinase.family"])
            results.rename(columns=RENAME_COLUMNS_REMOTE["kinases"], inplace=True)
            return results
        except SwaggerMappingError as e:
            _logger.error(e)

    def all_kinases(self, group=None, family=None, species=None):

        try:
            results = (
                self.__client.Information.get_kinase_names(
                    kinase_group=group, kinase_family=family, species=species
                )
                .response()
                .result
            )
            results = self._abc_to_dataframe(results)
            results.rename(columns=RENAME_COLUMNS_REMOTE["kinases"], inplace=True)
            return results
        except SwaggerMappingError as e:
            _logger.error(e)

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        try:
            results = []
            for kinase_id in kinase_ids:
                result = (
                    self.__client.Information.get_kinase_information(
                        kinase_ID=[kinase_id]
                    )
                    .response()
                    .result
                )
                result_df = self._abc_to_dataframe(result)
                results.append(result_df)
            results = pd.concat(results)
            results.rename(columns=RENAME_COLUMNS_REMOTE["kinases"], inplace=True)
            return results
        except SwaggerMappingError as e:
            _logger.error(e)

    def from_kinase_names(self, kinase_names, species=None):

        if isinstance(kinase_names, str):
            kinase_names = [kinase_names]

        results = []
        for kinase_name in kinase_names:
            try:
                result = (
                    self.__client.Information.get_kinase_ID(
                        kinase_name=kinase_name, species=species
                    )
                    .response()
                    .result
                )
                result_df = self._abc_to_dataframe(result)
                results.append(result_df)
            except SwaggerMappingError as e:
                _logger.error(f"Kinase name {kinase_name}: {e}")

        if len(results) > 0:
            kinases = pd.concat(results)
            # Kinase IDs can occur multiple times if the input kinase names describe the same kinase, thus drop duplicates
            kinases = kinases.drop_duplicates("kinase_ID").reset_index(drop=True)
            kinases.rename(columns=RENAME_COLUMNS_REMOTE["kinases"], inplace=True)
            return kinases


class Ligands(LigandsProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_ligands(self):

        # Get all kinase IDs
        kinases_remote = Kinases(self.__client)
        kinases = kinases_remote.all_kinases()
        # Get ligand results using the swagger API
        ligands_result = (
            self.__client.Ligands.get_ligands_list(
                kinase_ID=kinases["kinase.id"].to_list()
            )
            .response()
            .result
        )
        # Convert list of abstract base classes to DataFrame
        ligands_df = self._abc_to_dataframe(ligands_result)
        # Rename columns
        ligands_df.rename(columns=RENAME_COLUMNS_REMOTE["ligands"], inplace=True)
        return ligands_df

    def from_kinase_ids(self, kinase_ids):

        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        # Get ligands for each kinase ID
        ligand_list = [self._from_kinase_id(kinase_id) for kinase_id in kinase_ids]
        # Remove None values
        ligand_list = [
            ligands_df for ligands_df in ligand_list if ligands_df is not None
        ]

        if len(ligand_list) > 0:
            ligands_df = pd.concat(ligand_list)
            return ligands_df
        else:
            _log_error_empty_query_results()

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
            # Get ligand results using the swagger API
            ligands_result = (
                self.__client.Ligands.get_ligands_list(kinase_ID=[kinase_id])
                .response()
                .result
            )
            # Convert list of abstract base classes to DataFrame
            ligands_df = self._abc_to_dataframe(ligands_result)
            # Rename columns
            ligands_df.rename(columns=RENAME_COLUMNS_REMOTE["ligands"], inplace=True)
            # Add kinase ID to indicate which query was used to retrieve the results
            ligands_df["kinase.id (query)"] = kinase_id
            return ligands_df
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Kinase ID {kinase_id}: {e}")

    def from_kinase_names(self, kinase_names):

        if isinstance(kinase_names, str):
            kinase_names = [kinase_names]
        # Get kinase IDs for input kinase names (remotely)
        # Note: One kinase name can be linked to multiple kinase IDs (due to multiple species)
        kinases_remote = Kinases(self.__client)
        kinases = kinases_remote.from_kinase_names("BMX")[
            ["kinase.id", "kinase.name", "species.klifs"]
        ]
        # Rename columns to indicate columns involved in query
        kinases.rename(
            columns={
                "kinase.id": "kinase.id (query)",
                "kinase.name": "kinase.name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )
        # Get ligands by kinase IDs
        ligands = self.from_kinase_ids(kinases["kinase.id (query)"].to_list())
        # Add kinase name and species details to rationalize kinase IDs
        ligands = ligands.merge(kinases, on="kinase.id (query)", how="left")
        return ligands

    def from_ligand_ids(self, ligand_ids):

        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        ligands_all = self.all_ligands()
        ligands = ligands_all[ligands_all["ligand.id"].isin(ligand_ids)]

        return ligands

    def from_ligand_pdbs(self, ligand_pdbs):

        if isinstance(ligand_pdbs, str):
            ligand_pdbs = [ligand_pdbs]

        ligands_all = self.all_ligands()
        ligands = ligands_all[ligands_all["ligand.pdb"].isin(ligand_pdbs)]

        return ligands


class Structures(StructuresProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_structures(self):
        # Get all kinase IDs
        kinases_remote = Kinases(self.__client)
        kinase_ids = kinases_remote.all_kinases()["kinase.id"].to_list()
        # Get all structures from these kinase IDs
        structures = self.from_kinase_ids(kinase_ids)
        return structures

    def from_structure_ids(self, structure_ids):
        if isinstance(structure_ids, int):
            structure_ids = [structure_ids]

        try:
            structures_result = (
                self.__client.Structures.get_structure_list(structure_ID=structure_ids)
                .response()
                .result
            )
            structures_df = self._abc_to_dataframe(structures_result)
            structures_df.rename(
                columns=RENAME_COLUMNS_REMOTE["structures"], inplace=True
            )
            return structures_df
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Structure ID {structure_ids}: {e}")

    def from_ligand_ids(self, ligand_ids):
        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Get ligand PDB IDs for ligand IDs
        remote_ligands = Ligands(self.__client)
        ligands = remote_ligands.from_ligand_ids(ligand_ids)
        ligand_pdbs = ligands["ligand.pdb"].to_list()

        structures = self.from_ligand_pdbs(ligand_pdbs)

        return structures

    def from_kinase_ids(self, kinase_ids):
        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        try:
            structures_result = (
                self.__client.Structures.get_structures_list(kinase_ID=kinase_ids)
                .response()
                .result
            )
            structures_df = self._abc_to_dataframe(structures_result)
            structures_df.rename(
                columns=RENAME_COLUMNS_REMOTE["structures"], inplace=True
            )
            return structures_df
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Kinase ID {kinase_ids}: {e}")

    def from_structure_pdbs(self, structure_pdbs):
        if isinstance(structure_pdbs, str):
            structure_pdbs = [structure_pdbs]

        try:
            structures_result = (
                self.__client.Structures.get_structures_pdb_list(
                    pdb_codes=structure_pdbs
                )
                .response()
                .result
            )
            structures_df = self._abc_to_dataframe(structures_result)
            structures_df.rename(
                columns=RENAME_COLUMNS_REMOTE["structures"], inplace=True
            )
            return structures_df
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Structure PDB {structure_pdbs}: {e}")

    def from_ligand_pdbs(self, ligand_pdbs):
        if isinstance(ligand_pdbs, str):
            ligand_pdbs = [ligand_pdbs]

        structures_all = self.all_structures()
        structures = structures_all[structures_all["ligand.pdb"].isin(ligand_pdbs)]

        return structures

    def from_kinase_names(self, kinase_names):
        if isinstance(kinase_names, str):
            kinase_names = [kinase_names]

        structures_all = self.all_structures()
        structures = structures_all[structures_all["kinase.name"].isin(kinase_names)]

        return structures


class Bioactivities(BioactivitiesProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_bioactivities(self):
        # Get all kinase IDs
        ligands_remote = Ligands(self.__client)
        ligand_ids = ligands_remote.all_ligands()["ligand.id"].to_list()
        # Get all bioactivities from these ligand IDs
        bioactivities = self.from_ligand_ids(ligand_ids)
        return bioactivities

    def from_kinase_ids(self, kinase_ids):
        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        # Get all kinase IDs
        ligands_remote = Ligands(self.__client)
        ligands = ligands_remote.from_kinase_ids(kinase_ids)
        # Get all bioactivities from these ligand IDs
        if ligands is not None:
            bioactivities = self.from_ligand_ids(ligands["ligand.id"].to_list())
            return bioactivities

    def from_ligand_ids(self, ligand_ids):
        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Get bioactivities for each ligand ID
        bioactivity_list = [self._from_ligand_id(ligand_id) for ligand_id in ligand_ids]
        # Remove None values
        bioactivity_list = [
            bioactivity_df
            for bioactivity_df in bioactivity_list
            if bioactivity_df is not None
        ]

        if len(bioactivity_list) > 0:
            bioactivities = pd.concat(bioactivity_list)
            bioactivities.reset_index(drop=True, inplace=True)
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
            bioactivity_result = (
                self.__client.Ligands.get_bioactivity_list_id(ligand_ID=ligand_id)
                .response()
                .result
            )
            bioactivity_df = self._abc_to_dataframe(bioactivity_result)
            bioactivity_df.rename(
                columns=RENAME_COLUMNS_REMOTE["bioactivities"], inplace=True
            )
            bioactivity_df["ligand.id (query)"] = ligand_id
            return bioactivity_df
        except SwaggerMappingError as e:
            _logger.error(f"Ligand ID {ligand_id}: {e}")


class Interactions(InteractionsProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    @property
    def interaction_types(self):
        interaction_types_result = (
            self.__client.Interactions.get_interactions_get_types().response().result
        )
        interaction_types_df = self._abc_to_dataframe(interaction_types_result)
        interaction_types_df.rename(
            columns=RENAME_COLUMNS_REMOTE["interactions"], inplace=True
        )
        return interaction_types_df

    def all_interactions(self):
        # Get all structure IDs
        structures_remote = Structures(self.__client)
        structure_ids = structures_remote.all_structures()["structure.id"].to_list()
        # Get all interactions from these structures IDs
        interactions = self.from_structure_ids(structure_ids)
        return interactions

    def from_structure_ids(self, structure_ids):
        if isinstance(structure_ids, int):
            structure_ids = [structure_ids]
        try:
            interactions_result = (
                self.__client.Interactions.get_interactions_get_IFP(
                    structure_ID=structure_ids
                )
                .response()
                .result
            )
            interactions_df = self._abc_to_dataframe(interactions_result)
            interactions_df.rename(
                columns=RENAME_COLUMNS_REMOTE["interactions"], inplace=True
            )
            return interactions_df
        except (SwaggerMappingError, ValueError) as e:
            _logger.error(f"Structure ID {structure_ids}: {e}")

    def from_ligand_ids(self, ligand_ids):
        if isinstance(ligand_ids, int):
            ligand_ids = [ligand_ids]

        # Get structure IDs from ligand IDs
        structures_remote = Structures(self.__client)
        structures = structures_remote.from_ligand_ids(ligand_ids)
        # Get interactions from these structure IDs
        if structures is not None:
            interactions = self.from_structure_ids(structures["structure.id"].to_list())
            return interactions

    def from_kinase_ids(self, kinase_ids):
        if isinstance(kinase_ids, int):
            kinase_ids = [kinase_ids]

        # Get structure IDs from ligand IDs
        structures_remote = Structures(self.__client)
        structures = structures_remote.from_kinase_ids(kinase_ids)
        # Get interactions from these structure IDs
        if structures is not None:
            interactions = self.from_structure_ids(structures["structure.id"].to_list())
            return interactions


class Coordinates(CoordinatesProvider):
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
        if not text:  # TODO Ask Albert why no remote water
            raise ValueError(
                f"Entity {entity} is not available remotely but we could ask Albert to add this."
            )

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

        # Get structure metadata
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

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Get text
        text = self.fetch(structure_id, entity, input_format, output_format="text")
        if not text:  # TODO Ask Albert why no remote water
            raise ValueError(
                f"Entity {entity} is not available remotely but we could ask Albert to add this."
            )

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
            return (
                self.__client.Structures.get_structure_get_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "complex" and input_format == "pdb":
            return (
                self.__client.Structures.get_structure_get_pdb_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "ligand" and input_format == "mol2":
            return (
                self.__client.Structures.get_structure_get_ligand(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "pocket" and input_format == "mol2":
            return (
                self.__client.Structures.get_structure_get_pocket(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "protein" and input_format == "mol2":
            return (
                self.__client.Structures.get_structure_get_protein(
                    structure_ID=structure_id
                )
                .response()
                .result
            )

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
                col_names=[
                    "atom_id",
                    "atom_name",
                    "x",
                    "y",
                    "z",
                    "atom_type",
                    "subst_id",
                    "subst_name",
                    "charge",
                    "backbone",
                ],
                col_types=[int, str, float, float, float, str, int, str, float, str],
            )
        except ValueError:
            mol2_df = pmol._construct_df(
                mol2_text.splitlines(True),
                col_names=[
                    "atom_id",
                    "atom_name",
                    "x",
                    "y",
                    "z",
                    "atom_type",
                    "subst_id",
                    "subst_name",
                    "charge",
                ],
                col_types=[int, str, float, float, float, str, int, str, float],
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
