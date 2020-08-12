"""
remote.py

Defines remote KLIFS session.
"""

import logging

from bravado_core.exception import SwaggerMappingError
import pandas as pd

from .utils import KLIFS_CLIENT
from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    CoordinatesProvider,
)
from .utils import _abc_idlist_to_dataframe, _log_error_empty_query_results
from .utils import RENAME_COLUMNS_REMOTE_KINASE, RENAME_COLUMNS_REMOTE_LIGAND

_logger = logging.getLogger(__name__)


class Kinases(KinasesProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_kinase_groups(self):

        results = self.__client.Information.get_kinase_groups().response().result
        kinase_groups = pd.DataFrame(results, columns=["kinase.group"])
        kinase_groups.rename(columns=RENAME_COLUMNS_REMOTE_KINASE, inplace=True)
        return kinase_groups

    def all_kinase_families(self, group=None):

        try:
            results = (
                self.__client.Information.get_kinase_families(kinase_group=group)
                .response()
                .result
            )
            results = pd.DataFrame(results, columns=["kinase.family"])
            results.rename(columns=RENAME_COLUMNS_REMOTE_KINASE, inplace=True)
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
            results = _abc_idlist_to_dataframe(results)
            results.rename(columns=RENAME_COLUMNS_REMOTE_KINASE, inplace=True)
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
                result_df = _abc_idlist_to_dataframe(result)
                results.append(result_df)
            results = pd.concat(results)
            results.rename(columns=RENAME_COLUMNS_REMOTE_KINASE, inplace=True)
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
                result_df = _abc_idlist_to_dataframe(result)
                results.append(result_df)
            except SwaggerMappingError as e:
                _logger.error(f"Kinase name {kinase_name}: {e}")

        if len(results) > 0:
            kinases = pd.concat(results)
            # Kinase IDs can occur multiple times if the input kinase names describe the same kinase, thus drop duplicates
            kinases = kinases.drop_duplicates("kinase_ID").reset_index(drop=True)
            kinases.rename(columns=RENAME_COLUMNS_REMOTE_KINASE, inplace=True)
            return kinases


class Ligands(LigandsProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client

    def all_ligands(self):

        # Get all kinase IDs
        kinases_remote = Kinases(KLIFS_CLIENT)
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
        ligands_df = _abc_idlist_to_dataframe(ligands_result)
        # Rename columns
        ligands_df.rename(columns=RENAME_COLUMNS_REMOTE_LIGAND, inplace=True)
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
            ligands_df = _abc_idlist_to_dataframe(ligands_result)
            # Rename columns
            ligands_df.rename(columns=RENAME_COLUMNS_REMOTE_LIGAND, inplace=True)
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
        kinases_remote = Kinases(KLIFS_CLIENT)
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

        if ligands.shape[0] > 0:
            return ligands

    def from_ligand_pdbs(self, ligand_pdbs):

        if isinstance(ligand_pdbs, str):
            ligand_pdbs = [ligand_pdbs]

        ligands_all = self.all_ligands()
        ligands = ligands_all[ligands_all["ligand.pdb"].isin(ligand_pdbs)]

        if ligands.shape[0] > 0:
            return ligands


class Bioactivities(BioactivitiesProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client


class Structures(StructuresProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client


class Interactions(InteractionsProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client


class Coordinates(CoordinatesProvider):
    def __init__(self, client):
        super().__init__()
        self.__client = client
