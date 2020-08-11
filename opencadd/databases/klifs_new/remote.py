"""
remote.py

Defines remote KLIFS session.
"""

import logging

from bravado_core.exception import SwaggerMappingError
import pandas as pd

from .core import (
    KinasesProvider,
    LigandsProvider,
    StructuresProvider,
    BioactivitiesProvider,
    InteractionsProvider,
    CoordinatesProvider,
)
from .utils import _abc_idlist_to_dataframe
from .utils import RENAME_COLUMNS_REMOTE_KINASE

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
                _logger.error(f"Kinase {kinase_name}: {e}")

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
