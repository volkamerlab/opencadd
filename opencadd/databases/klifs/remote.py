"""
opencadd.databases.klifs.remote

Defines a remote KLIFS session.
"""

import logging
from copy import deepcopy
from pathlib import Path

import pandas as pd
from bravado.client import SwaggerClient

from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    PocketsProvider,
    CoordinatesProvider,
    DrugsProvider,
    StructureConformationsProvider,
    StructureModifiedResiduesProvider,
)
from .schema import FIELDS
from .utils import metadata_to_filepath, silence_logging
from opencadd.io.dataframe import DataFrame
from opencadd.io.rdkit import Rdkit

_logger = logging.getLogger(__name__)


class SerializableSwaggerClient(SwaggerClient):
    # Since they are using __attributes to mangle the namespace
    # we need to hardcode the parent class name in saved attributes
    # (only also_return_response in this case)
    # Sorry about the hackiness :)

    def __setstate__(self, state):
        self._SwaggerClient__also_return_response = state["also_return_response"]
        self.swagger_spec = state["swagger_spec"]

    def __getstate__(self, *args):
        return {
            "also_return_response": deepcopy(self._SwaggerClient__also_return_response),
            "swagger_spec": deepcopy(self.swagger_spec),
        }


KLIFS_API_DEFINITIONS = "https://dev.klifs.net/swagger_v2/swagger.json"
KLIFS_CLIENT = SerializableSwaggerClient.from_url(
    KLIFS_API_DEFINITIONS, config={"validate_responses": False}
)
# Maximum chunk size for input values for remote client endpoints
N_CHUNK = 10_000


class RemoteInitializer:
    """
    Base class used to define __init__ for all remote classes.

    Attributes
    ----------
    _client : bravado.client.SwaggerClient
        KLIFS client.
    """

    def __init__(self, client, *args, **kwargs):
        self._client = client


class Kinases(RemoteInitializer, KinasesProvider):
    """
    Extends KinasesProvider to provide remote kinases requests.
    Refer to KinasesProvider documentation for more information:
    opencadd.databases.klifs.core.KinasesProvider
    """

    def all_kinase_groups(self):
        # Use KLIFS API
        result = self._client.Information.get_kinase_groups().response().result
        # Convert list to DataFrame (1 column)
        column_name = list(FIELDS.oc_name_to_type("kinase_groups").keys())[0]
        column_dtype = list(FIELDS.oc_name_to_type("kinase_groups").values())[0]
        kinase_groups = pd.DataFrame({column_name: pd.Series(result, dtype=column_dtype)})
        return kinase_groups

    def all_kinase_families(self, groups=None):
        groups = self._ensure_list(groups)
        # Use KLIFS API
        result = (
            self._client.Information.get_kinase_families(kinase_group=groups).response().result
        )
        # Convert list to DataFrame (1 column)
        column_name = list(FIELDS.oc_name_to_type("kinase_families").keys())[0]
        column_dtype = list(FIELDS.oc_name_to_type("kinase_families").values())[0]
        kinase_families = pd.DataFrame({column_name: pd.Series(result, dtype=column_dtype)})
        return kinase_families

    def all_kinases(self, groups=None, families=None, species=None):
        groups = self._ensure_list(groups)
        families = self._ensure_list(families)
        # Use KLIFS API
        result = (
            self._client.Information.get_kinase_names(
                kinase_group=groups, kinase_family=families, species=species
            )
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        kinases = self._abc_to_dataframe(result)
        # Standardize DataFrame
        kinases = self._standardize_dataframe(
            kinases,
            FIELDS.oc_name_to_type("kinases_all"),
            FIELDS.remote_to_oc_names("kinases_all"),
        )
        return kinases

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Use KLIFS API
        result = (
            self._client.Information.get_kinase_information(kinase_ID=kinase_klifs_ids)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        kinases = self._abc_to_dataframe(result)
        # Standardize DataFrame
        kinases = self._standardize_dataframe(
            kinases, FIELDS.oc_name_to_type("kinases"), FIELDS.remote_to_oc_names("kinases")
        )
        return kinases

    def by_kinase_name(self, kinase_names, species=None):
        kinase_names = self._ensure_list(kinase_names)
        # Use KLIFS API
        result = (
            self._client.Information.get_kinase_ID(kinase_name=kinase_names, species=species)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        kinases = self._abc_to_dataframe(result)
        # Standardize DataFrame
        kinases = self._standardize_dataframe(
            kinases, FIELDS.oc_name_to_type("kinases"), FIELDS.remote_to_oc_names("kinases")
        )
        return kinases


class Ligands(RemoteInitializer, LigandsProvider):
    """
    Extends LigandsProvider to provide remote ligands requests.
    Refer to LigandsProvider documentation for more information:
    opencadd.databases.klifs.core.LigandsProvider
    """

    def all_ligands(self):
        # Use KLIFS API: Get all kinase KLIFS IDs
        kinases_remote = Kinases(self._client)
        kinases = kinases_remote.all_kinases()
        # Use KLIFS API: Get ligands
        kinase_klifs_ids = kinases["kinase.klifs_id"].to_list()
        result = (
            self._client.Ligands.get_ligands_list(kinase_ID=kinase_klifs_ids).response().result
        )
        # Convert list of ABC objects to DataFrame
        ligands = self._abc_to_dataframe(result)
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands, FIELDS.oc_name_to_type("ligands"), FIELDS.remote_to_oc_names("ligands")
        )
        return ligands

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        # Use KLIFS API (send requests iteratively)
        ligands = self._multiple_remote_requests(self._by_kinase_klifs_id, kinase_klifs_ids)
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands,
            FIELDS.oc_name_to_type("ligands", {"kinase.klifs_id (query)": "int32"}),
            FIELDS.remote_to_oc_names("ligands"),
        )
        return ligands

    def _by_kinase_klifs_id(self, kinase_klifs_id):
        """
        Get ligands by kinase KLIFS ID.

        Parameters
        ----------
        kinase_klifs_id : int
            Kinase KLIFS ID.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """

        # Use KLIFS API
        result = (
            self._client.Ligands.get_ligands_list(kinase_ID=[kinase_klifs_id]).response().result
        )
        # Convert list of ABC objects to DataFrame
        ligands = self._abc_to_dataframe(result)
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands, FIELDS.oc_name_to_type("ligands"), FIELDS.remote_to_oc_names("ligands")
        )
        # Rename column to indicate query key
        ligands["kinase.klifs_id (query)"] = kinase_klifs_id
        return ligands

    def by_kinase_name(self, kinase_names):
        kinase_names = self._ensure_list(kinase_names)
        # Use KLIFS API: Get kinase KLIFS IDs for input kinase names (remotely)
        # Note: One kinase name can be linked to multiple kinase KLIFS IDs (due to multiple species)
        _logger.info(f"Fetch kinase KLIFS IDs for input kinase names...")
        kinases_remote = Kinases(self._client)
        kinases = kinases_remote.by_kinase_name(kinase_names)
        # Select and rename columns to indicate columns involved in query
        kinases = kinases[
            ["kinase.klifs_id", "kinase.klifs_name", "kinase.gene_name", "species.klifs"]
        ]
        kinases.rename(
            columns={
                "kinase.klifs_id": "kinase.klifs_id (query)",
                "kinase.klifs_name": "kinase.klifs_name (query)",
                "kinase.gene_name": "kinase.gene_name (query)",
                "species.klifs": "species.klifs (query)",
            },
            inplace=True,
        )
        # Use KLIFS API: Get ligands by kinase KLIFS IDs
        _logger.info(f"Fetch ligands based on these KLIFS IDs...")
        kinase_klifs_ids = kinases["kinase.klifs_id (query)"].to_list()
        ligands = self.by_kinase_klifs_id(kinase_klifs_ids)
        # Add kinase name and species details to rationalize kinase KLIFS IDs
        ligands = ligands.merge(kinases, on="kinase.klifs_id (query)", how="left")
        return ligands

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        ligand_klifs_ids = self._ensure_list(ligand_klifs_ids)
        # Use KLIFS API: Get all ligands
        ligands = self.all_ligands()
        # Select ligands by ligand KLIFS IDs
        ligands = ligands[ligands["ligand.klifs_id"].isin(ligand_klifs_ids)]
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands, FIELDS.oc_name_to_type("ligands"), FIELDS.remote_to_oc_names("ligands")
        )
        return ligands

    def by_ligand_expo_id(self, ligand_expo_ids):
        ligand_expo_ids = self._ensure_list(ligand_expo_ids)
        # Use KLIFS API: Get all ligands
        ligands = self.all_ligands()
        # Select ligands by Ligand Expo IDs
        ligands = ligands[ligands["ligand.expo_id"].isin(ligand_expo_ids)]
        # Standardize DataFrame
        ligands = self._standardize_dataframe(
            ligands, FIELDS.oc_name_to_type("ligands"), FIELDS.remote_to_oc_names("ligands")
        )
        return ligands


class Structures(RemoteInitializer, StructuresProvider):
    """
    Extends StructuresProvider to provide remote structures requests.
    Refer to StructuresProvider documentation for more information:
    opencadd.databases.klifs.core.StructuresProvider
    """

    def all_structures(self):
        # Use KLIFS API: Get all kinase KLIFS IDs
        kinases_remote = Kinases(self._client)
        kinases = kinases_remote.all_kinases()
        # Use KLIFS API: Get all structures from these kinase KLIFS IDs
        kinase_klifs_ids = kinases["kinase.klifs_id"].to_list()
        structures = self.by_kinase_klifs_id(kinase_klifs_ids)
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures

    def by_structure_klifs_id(self, structure_klifs_ids):
        structure_klifs_ids = self._ensure_list(structure_klifs_ids)
        # Use KLIFS API
        result = (
            self._client.Structures.get_structure_list(structure_ID=structure_klifs_ids)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        structures = self._abc_to_dataframe(result)
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        ligand_klifs_ids = self._ensure_list(ligand_klifs_ids)
        # Use KLIFS API: Get all structures
        structures = self.all_structures()
        # Select structures by ligand KLIFS IDs
        structures = structures[structures["ligand.klifs_id"].isin(ligand_klifs_ids)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Use KLIFS API
        result = (
            self._client.Structures.get_structures_list(kinase_ID=kinase_klifs_ids)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        structures = self._abc_to_dataframe(result)
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures

    def by_structure_pdb_id(
        self, structure_pdb_ids, structure_alternate_model=None, structure_chain=None
    ):
        structure_pdb_ids = self._ensure_list(structure_pdb_ids)
        # Use KLIFS API
        result = (
            self._client.Structures.get_structures_pdb_list(pdb_codes=structure_pdb_ids)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        structures = self._abc_to_dataframe(result)
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        # If only one structure PDB ID is given, check alternate model and chain filters
        if len(structure_pdb_ids) == 1:
            structures = self._filter_pdb_by_alt_chain(
                structures, structure_alternate_model, structure_chain
            ).reset_index(drop=True)
        return structures

    def by_ligand_expo_id(self, ligand_expo_ids):
        ligand_expo_ids = self._ensure_list(ligand_expo_ids)
        # Use KLIFS API: Get all structures
        structures = self.all_structures()
        # Select structures by Ligand Expo IDs
        structures = structures[structures["ligand.expo_id"].isin(ligand_expo_ids)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures

    def by_kinase_name(self, kinase_names):
        kinase_names = self._ensure_list(kinase_names)
        # Use KLIFS API: Get all structures
        structures = self.all_structures()
        # Select structures by kinase names
        structures = structures[structures["kinase.klifs_name"].isin(kinase_names)]
        # Standardize DataFrame
        structures = self._standardize_dataframe(
            structures,
            FIELDS.oc_name_to_type("structures"),
            FIELDS.remote_to_oc_names("structures"),
        )
        return structures


class Bioactivities(RemoteInitializer, BioactivitiesProvider):
    """
    Extends BioactivitiesProvider to provide remote bioactivities requests.
    Refer to BioactivitiesProvider documentation for more information:
    opencadd.databases.klifs.core.BioactivitiesProvider
    """

    def all_bioactivities(self, _top_n=None):
        # Use KLIFS API: Get all kinase KLIFS IDs
        ligands_remote = Ligands(self._client)
        ligands = ligands_remote.all_ligands()
        # Optional: Select top n ligands for bioactivity query!
        if _top_n:
            ligands = ligands[:_top_n]
        # Use KLIFS API: Get all bioactivities from these ligand KLIFS IDs
        ligand_klifs_ids = ligands["ligand.klifs_id"].to_list()
        # Many ligands do not have bioactivities in ChEMBL,
        # Thus, disable logging messages for this query
        with silence_logging():
            bioactivities = self.by_ligand_klifs_id(ligand_klifs_ids)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities"),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        return bioactivities

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Use KLIFS API: Get all kinase KLIFS IDs
        ligands_remote = Ligands(self._client)
        ligands = ligands_remote.by_kinase_klifs_id(kinase_klifs_ids)
        # Use KLIFS API: Get all bioactivities from these ligand KLIFS IDs
        ligand_klifs_ids = ligands["ligand.klifs_id"].to_list()
        bioactivities = self.by_ligand_klifs_id(ligand_klifs_ids)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities"),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        return bioactivities

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        # Use KLIFS API (send requests iteratively)
        bioactivities = self._multiple_remote_requests(self._by_ligand_klifs_id, ligand_klifs_ids)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities", {"ligand.klifs_id (query)": "int32"}),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        return bioactivities

    def by_ligand_expo_id(self, ligand_expo_id):
        # Use KLIFS API (send requests iteratively)
        bioactivities = self._multiple_remote_requests(self._by_ligand_expo_id, ligand_expo_id)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities", {"ligand.expo_id (query)": "string"}),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        return bioactivities

    def _by_ligand_klifs_id(self, ligand_klifs_id):
        """
        Get bioactivities by ligand KLIFS ID.

        Parameters
        ----------
        ligand_klifs_id : int
            Ligand KLIFS ID.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """

        # Use KLIFS API
        result = (
            self._client.Ligands.get_bioactivity_list_id(ligand_ID=ligand_klifs_id)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        bioactivities = self._abc_to_dataframe(result)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities"),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        # Rename column to indicate query key
        bioactivities["ligand.klifs_id (query)"] = ligand_klifs_id
        return bioactivities

    def _by_ligand_expo_id(self, ligand_expo_id):
        """
        Get bioactivities by ligand Expo ID.

        Parameters
        ----------
        ligand_expo_id : int
            Ligand Expo ID.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """

        # Use KLIFS API
        result = (
            self._client.Ligands.get_bioactivity_list_pdb(ligand_PDB=ligand_expo_id)
            .response()
            .result
        )
        # Convert list of ABC objects to DataFrame
        bioactivities = self._abc_to_dataframe(result)
        # Standardize DataFrame
        bioactivities = self._standardize_dataframe(
            bioactivities,
            FIELDS.oc_name_to_type("bioactivities"),
            FIELDS.remote_to_oc_names("bioactivities"),
        )
        # Rename column to indicate query key
        bioactivities["ligand.expo_id (query)"] = ligand_expo_id
        return bioactivities


class Interactions(RemoteInitializer, InteractionsProvider):
    """
    Extends InteractionsProvider to provide remote kinases requests.
    Refer to InteractionsProvider documentation for more information:
    opencadd.databases.klifs.core.InteractionsProvider
    """

    @property
    def interaction_types(self):
        # Use KLIFS API
        result = self._client.Interactions.get_interactions_get_types().response().result
        # Convert list of ABC objects to DataFrame
        interaction_types = self._abc_to_dataframe(result)
        # Standardize DataFrame
        interaction_types = self._standardize_dataframe(
            interaction_types,
            FIELDS.oc_name_to_type("interaction_types"),
            FIELDS.remote_to_oc_names("interaction_types"),
        )
        return interaction_types

    def all_interactions(self):
        # Use KLIFS API: Get all structure KLIFS IDs
        structures_remote = Structures(self._client)
        structures = structures_remote.all_structures()
        # Use KLIFS API: Get all interactions from these structures KLIFS IDs
        structure_klifs_ids = structures["structure.klifs_id"].to_list()
        interactions = self.by_structure_klifs_id(structure_klifs_ids)
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions,
            FIELDS.oc_name_to_type("interactions"),
            FIELDS.remote_to_oc_names("interactions"),
        )
        return interactions

    def by_structure_klifs_id(self, structure_klifs_ids):
        structure_klifs_ids = self._ensure_list(structure_klifs_ids)
        interactions = []
        # NOTE: Chunk IDs to avoid HTTPURITooLong
        for i in range(0, len(structure_klifs_ids), N_CHUNK):
            _structure_klifs_ids = structure_klifs_ids[i : i + N_CHUNK]
            # Use KLIFS API
            result = (
                self._client.Interactions.get_interactions_get_IFP(
                    structure_ID=_structure_klifs_ids
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame
            _interactions = self._abc_to_dataframe(result)
            # Standardize DataFrame
            _interactions = self._standardize_dataframe(
                _interactions,
                FIELDS.oc_name_to_type("interactions"),
                FIELDS.remote_to_oc_names("interactions"),
            )
            interactions.append(_interactions)
        interactions = pd.concat(interactions)
        return interactions

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        ligand_klifs_ids = self._ensure_list(ligand_klifs_ids)
        # Use KLIFS API: Get structure KLIFS IDs from ligand KLIFS IDs
        structures_remote = Structures(self._client)
        structures = structures_remote.by_ligand_klifs_id(ligand_klifs_ids)
        # Use KLIFS API: Get interactions from these structure KLIFS IDs
        structure_klifs_ids = structures["structure.klifs_id"].to_list()
        interactions = self.by_structure_klifs_id(structure_klifs_ids)
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions,
            FIELDS.oc_name_to_type("interactions"),
            FIELDS.remote_to_oc_names("interactions"),
        )
        return interactions

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        kinase_klifs_ids = self._ensure_list(kinase_klifs_ids)
        # Use KLIFS API: Get structure KLIFS IDs from ligand KLIFS IDs
        structures_remote = Structures(self._client)
        structures = structures_remote.by_kinase_klifs_id(kinase_klifs_ids)
        # Use KLIFS API: Get interactions from these structure KLIFS IDs
        structure_klifs_ids = structures["structure.klifs_id"].to_list()
        interactions = self.by_structure_klifs_id(structure_klifs_ids)
        # Standardize DataFrame
        interactions = self._standardize_dataframe(
            interactions,
            FIELDS.oc_name_to_type("interactions"),
            FIELDS.remote_to_oc_names("interactions"),
        )
        return interactions


class Pockets(RemoteInitializer, PocketsProvider):
    """
    Extends PocketsProvider to provide remote pocket requests.
    Refer to PocketsProvider documentation for more information:
    opencadd.databases.klifs.core.PocketsProvider
    """

    def by_structure_klifs_id(self, structure_klifs_id):
        # Use KLIFS API
        result = (
            self._client.Interactions.get_interactions_match_residues(
                structure_ID=structure_klifs_id
            )
            .response()
            .result
        )
        # Convert to DataFrame and formatting
        pocket = pd.DataFrame(result)
        # Standardize DataFrame
        pocket = self._standardize_dataframe(
            pocket, FIELDS.oc_name_to_type("pockets"), FIELDS.remote_to_oc_names("pockets")
        )
        # Add KLIFS region and color  TODO not so nice to have this after standardization
        pocket = self._add_klifs_region_details(pocket)
        return pocket


class Coordinates(RemoteInitializer, CoordinatesProvider):
    """
    Extends CoordinatesProvider to provide remote coordinates requests,
    i.e. fetching and saving structural data (coordinates).
    Refer to CoordinatesProvider documentation for more information:
    opencadd.databases.klifs.core.CoordinatesProvider
    """

    def to_text(self, structure_klifs_id, entity="complex", extension="mol2"):
        self._raise_invalid_extension(extension)

        if entity == "complex" and extension == "mol2":
            text = (
                self._client.Structures.get_structure_get_complex(structure_ID=structure_klifs_id)
                .response()
                .result
            )
        elif entity == "complex" and extension == "pdb":
            text = (
                self._client.Structures.get_structure_get_pdb_complex(
                    structure_ID=structure_klifs_id
                )
                .response()
                .result
            )
        elif entity == "ligand" and extension == "mol2":
            text = (
                self._client.Structures.get_structure_get_ligand(structure_ID=structure_klifs_id)
                .response()
                .result
            )
        elif entity == "pocket" and extension == "mol2":
            text = (
                self._client.Structures.get_structure_get_pocket(structure_ID=structure_klifs_id)
                .response()
                .result
            )
        elif entity == "protein" and extension == "mol2":
            text = (
                self._client.Structures.get_structure_get_protein(structure_ID=structure_klifs_id)
                .response()
                .result
            )
        else:
            raise ValueError(f"Entity {entity} is not available or not available remotely.")

        if text:
            return text
        else:
            raise ValueError(f"Data could not be fetched.")

    def to_dataframe(self, structure_klifs_id, entity="complex", extension="mol2"):
        text = self.to_text(structure_klifs_id, entity, extension)
        dataframe = DataFrame.from_text(text, extension)
        dataframe = self._add_residue_klifs_ids(dataframe, structure_klifs_id)
        return dataframe

    def to_rdkit(self, structure_klifs_id, entity="complex", extension="mol2", compute2d=True):
        text = self.to_text(structure_klifs_id, entity, extension)
        rdkit_mol = Rdkit.from_text(text, extension, compute2d)
        return rdkit_mol

    def to_pdb(self, structure_klifs_id, output_path, entity="complex", in_dir=False):
        """
        Save structural data as pdb file.

        Parameters
        ----------
        structure_klifs_id : str
            Structure KLIFS ID.
        output_path : pathlib.Path or str
            Path to output folder.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        in_dir : bool
            Save file in KLIFS directory structure (default: False).

        Returns
        -------
        pathlib.Path
            Path to file.

        Raises
        ------
        ValueError
            If input yields not result.
        """

        return self._to_file(structure_klifs_id, output_path, entity, "pdb", in_dir)

    def to_mol2(self, structure_klifs_id, output_path, entity="complex", in_dir=False):
        """
        Save structural data as mol2 file.

        Parameters
        ----------
        structure_klifs_id : str
            Structure KLIFS ID.
        output_path : pathlib.Path or str
            Path to output folder.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        in_dir : bool
            Save file in KLIFS directory structure (default: False).

        Returns
        -------
        pathlib.Path
            Path to file.

        Raises
        ------
        ValueError
            If input yields not result.
        """

        return self._to_file(structure_klifs_id, output_path, entity, "mol2", in_dir)

    def _to_file(self, structure_klifs_id, output_path, entity, extension, in_dir=False):
        """
        Save structural data to file.

        Parameters
        ----------
        structure_klifs_id : str
            Structure KLIFS ID.
        output_path : pathlib.Path or str
            Path to output folder.
        entity : str
            Structural entity: complex, ligand, pocket, or protein.
        extension : str
            Input file format (fetched from KLIFS): mol2 or pdb (only for entity=complex).
        in_dir : bool
            Save file in KLIFS directory structure (default: False).

        Returns
        -------
        pathlib.Path
            Path to file.

        Raises
        ------
        ValueError
            If input yields not result.
        """

        # self._check_parameter_validity(entity, extension)
        output_path = Path(output_path)

        # Use KLIFS API: Get structure metadata
        structures_remote = Structures(self._client)
        metadata = structures_remote.by_structure_klifs_id(structure_klifs_id).iloc[0]

        # Set up output path (metadata in the form of directory structure or file name)
        output_path = metadata_to_filepath(
            output_path,
            metadata["species.klifs"].upper(),
            metadata["kinase.klifs_name"],
            metadata["structure.pdb_id"],
            metadata["structure.alternate_model"],
            metadata["structure.chain"],
            entity,
            extension,
            in_dir,
        )

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Get text
        text = self.to_text(structure_klifs_id, entity, extension)

        # Save text to file
        with open(output_path, "w") as f:
            f.write(text)

        return output_path

    def _add_residue_klifs_ids(self, dataframe, structure_klifs_id):
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

        # Get pocket residue details: KLIFS and PDB residue IDs
        pockets_remote = Pockets(self._client)
        pocket = pockets_remote.by_structure_klifs_id(structure_klifs_id)
        # Merge DataFrames
        dataframe = dataframe.merge(pocket, on="residue.id", how="left")
        dataframe = dataframe.astype({"residue.klifs_id": "Int64"})

        return dataframe


class Drugs(RemoteInitializer, DrugsProvider):
    """
    Extends DrugsProvider to provide remote drug requests.
    Refer to DrugsProvider documentation for more information:
    opencadd.databases.klifs.core.DrugsProvider
    """

    def all_drugs(self):
        # Use KLIFS API
        result = self._client.Ligands.get_drug_list().response().result
        # Convert list of ABC objects to DataFrame
        drugs = self._abc_to_dataframe(result)
        # Standardize DataFrame
        drugs = self._standardize_dataframe(
            drugs, FIELDS.oc_name_to_type("drugs"), FIELDS.remote_to_oc_names("drugs")
        )
        return drugs


class StructureConformations(RemoteInitializer, StructureConformationsProvider):
    """
    Extends StructureConformationsProvider to provide remote conformations requests.
    Refer to StructureConformationsProvider documentation for more information:
    opencadd.databases.klifs.core.StructureConformationsProvider
    """

    def all_conformations(self):
        # Use KLIFS API: Get all structure KLIFS IDs
        structures_remote = Structures(self._client)
        structures = structures_remote.all_structures()
        # Use KLIFS API: Get all interactions from these structures KLIFS IDs
        structure_klifs_ids = structures["structure.klifs_id"].to_list()
        conformations = self.by_structure_klifs_id(structure_klifs_ids)
        # Standardize DataFrame
        conformations = self._standardize_dataframe(
            conformations,
            FIELDS.oc_name_to_type("structure_conformations"),
            FIELDS.remote_to_oc_names("structure_conformations"),
        )
        return conformations

    def by_structure_klifs_id(self, structure_klifs_ids):
        structure_klifs_ids = self._ensure_list(structure_klifs_ids)
        conformations = []
        # NOTE: Chunk IDs to avoid HTTPURITooLong
        for i in range(0, len(structure_klifs_ids), N_CHUNK):
            _structure_klifs_ids = structure_klifs_ids[i : i + N_CHUNK]
            # Use KLIFS API
            result = (
                self._client.Structures.get_structure_conformation(
                    structure_ID=_structure_klifs_ids
                )
                .response()
                .result
            )
            # Convert list of ABC objects to DataFrame
            _conformations = self._abc_to_dataframe(result)
            # Standardize DataFrame
            _conformations = self._standardize_dataframe(
                _conformations,
                FIELDS.oc_name_to_type("structure_conformations"),
                FIELDS.remote_to_oc_names("structure_conformations"),
            )
            conformations.append(_conformations)
        conformations = pd.concat(conformations)
        return conformations


class StructureModifiedResidues(RemoteInitializer, StructureModifiedResiduesProvider):
    """
    Extends StructureModifiedResiduesProvider to provide remote modified residues requests.
    Refer to StructureModifiedResiduesProvider documentation for more information:
    opencadd.databases.klifs.core.StructureModifiedResiduesProvider
    """

    def by_structure_klifs_id(self, structure_klifs_id):
        # Use KLIFS API
        result = (
            self._client.Structures.get_structure_modified_residues(
                structure_ID=structure_klifs_id
            )
            .response()
            .result
        )
        if len(result) == 0:
            # Usually we return errors if the DataFrame is empty because it means that the input
            # data was incorrect or KLIFS is missing data that should be there
            # In this case here, however, most structures won't have modifications and therefore
            # will have nothing to return
            # To make the transition between 0 and >0 modifications smoother, let's always return
            # DataFrames will all columns and 0 or more rows.
            fields = FIELDS.df
            column_names = fields[fields["field_type"] == "structure_modified_residues"][
                "opencadd.df_name"
            ].to_list()
            return pd.DataFrame([], columns=column_names)
        else:
            # Convert to DataFrame and formatting
            modifications = pd.DataFrame(result)
            # Standardize DataFrame
            modifications = self._standardize_dataframe(
                modifications,
                FIELDS.oc_name_to_type("structure_modified_residues"),
                FIELDS.remote_to_oc_names("structure_modified_residues"),
            )
            # Add KLIFS region and color  TODO not so nice to have this after standardization
            modifications = self._add_klifs_region_details(modifications)
            return modifications
