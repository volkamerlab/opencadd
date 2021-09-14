"""
opencadd.databases.klifs.core

Defines core classes and functions.
"""

import logging
import html

from bravado_core.exception import SwaggerMappingError
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .schema import POCKET_KLIFS_REGION_COLORS

_logger = logging.getLogger(__name__)


class BaseProvider:
    """
    Base class for KLIFS requests (local and remote).
    """

    @staticmethod
    def _ensure_list(value):
        """
        Cast value to list if it is not a list, else return as-is.
        """

        # TODO python iterator protocol
        # check items (if they are str or int)
        # https://stackoverflow.com/questions/16301253/what-exactly-is-pythons-iterator-protocol
        # test for behaviour (value.__iter__()) instead of type
        # can be removed once input lists are unpacked *kinases_ids
        if value is None:
            return value
        elif not isinstance(value, list):
            return [value]
        else:
            return value

    @staticmethod
    def _abc_to_dataframe(abc_object):
        """
        Transform ABC object into a DataFrame (needed for KLIFS API results).
        Note: These ABC objects are bravado wrappers for KLIFS responses.

        Parameters
        ----------
        abc_object : list of abc.IDList/KinaseInformation/ligandDetails/structureDetails
            List of labeled list objects from abstract base classes module.

        Returns
        -------
        pandas.DataFrame
            Table with list labels as column names.
        """

        result = abc_object

        keys = list(result[0])

        results_dict = {key: [] for key in keys}

        for result in abc_object:
            for key in keys:
                results_dict[key].append(result[key])

        return pd.DataFrame(results_dict).apply(html.unescape)

    @staticmethod
    def _map_old_to_new_column_names(dataframe, columns_mapping):
        """
        Rename column names that are returned from local or remote queries (old column names)
        into standardized column names to be used in this module (new column names).

        Parameters
        ----------
        dataframe : pandas.DataFrame
            KLIFS data.
        columns_mapping : dict
            Mapping of old to new column names. If None, no changes are made.

        Notes
        -----
        This method will only change remote data, since column name standardization for
        local data is already performed upon session initialization.
        """

        dataframe.rename(columns=columns_mapping, inplace=True)
        return dataframe

    @staticmethod
    def _add_missing_columns(dataframe, column_names):
        """
        Add missing columns in local or remote result that are available in result pendant
        from remote or local module, respectively.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Result from query in local or remote session.
        column_names : list of str
            Column names.

        Returns
        -------
        pandas.DataFrame
            Input DataFrame with additional columns representing missing data.
        """

        for column_name in column_names:
            if column_name not in dataframe.columns:
                dataframe[column_name] = None

        return dataframe

    @staticmethod
    def _standardize_column_values(dataframe):
        """
        Standardize column values to be consistent across local and remote return/response values.
        Mainly, this concerns the notation of missing values (0, "-", "", ...).

        Parameters
        ----------
        dataframe : pandas.DataFrame
            KLIFS data.

        Returns
        -------
        pandas.DataFrame
            KLIFS data.

        Notes
        -----
        This method will only change remote data, since column value standardization for
        local data is already performed upon session initialization.
        """

        if "structure.alternate_model" in dataframe.columns:
            # Remote
            dataframe["structure.alternate_model"].replace("", "-", inplace=True)
        if "ligand.expo_id" in dataframe.columns:
            # Remote
            dataframe["ligand.expo_id"].replace(0, "-", inplace=True)
        if "ligand_allosteric.expo_id" in dataframe.columns:
            # Remote
            dataframe["ligand_allosteric.expo_id"].replace(0, "-", inplace=True)
        if "structure.resolution" in dataframe.columns:
            # Remote
            dataframe["structure.resolution"].replace(0, np.nan, inplace=True)

        return dataframe

    def _standardize_dataframe(self, dataframe, columns, columns_mapping=None):
        """
        Standardize a DataFrame across local and remote query results.
        - Map old to new column names (applied to remote data only, local data is standardized upon
          session initialization already).
        - Add missing columns (if a column is missing in remote data, add empty column to local
          data - and vice versa).
        - Standardize column values (applies mostly to missing value notations).
        - Select and sort column names.
        - Drop duplicate rows.
        - Reset dataframe index.
        - Handle empty dataframes: raise ValueError.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Remote query result.
        columns : dict
            Column names and dtypes (in the order of interest for output).
        columns_mapping : dict or None
            Mapping of old to new column names. If None, no changes are made.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """

        # Standardize column names
        if columns_mapping:
            dataframe = self._map_old_to_new_column_names(dataframe, columns_mapping)

        # Add missing columns (None values)
        column_names = list(columns.keys())
        dataframe = self._add_missing_columns(dataframe, column_names)

        # Standardize column values
        dataframe = self._standardize_column_values(dataframe)

        # Standardize dtypes
        dataframe = dataframe.astype(columns)

        # Select and sort columns
        dataframe = dataframe[column_names].copy()

        # Drop duplicate rows
        dataframe.drop_duplicates(inplace=True)

        # Reset index
        dataframe.reset_index(drop=True, inplace=True)

        # Handle empty results
        if dataframe.shape[0] > 0:
            return dataframe
        else:
            raise ValueError(f"Input values yield no results.")

    def _multiple_remote_requests(self, function, iterator, *args, **kwargs):
        """
        Wrap remote requests using multiple inputs, where KLIFS API only allows a single input
        or where for some reason - in this package - single input requests are preferred over
        using the KLIFS API for multiple input requests.

        Parameters
        ----------
        function : function
            Function that returns result for a single input.
        iterator : int/str or list of int/str
            Iterator, here IDs or names.
            If single input is given (thus no iterator), it will be cast to a list.
        """

        iterator = self._ensure_list(iterator)

        # If single request fails, error will be raised. Catch errors in list.
        errors = []

        # Get result for each iterator element
        # If more than 10 requests, print progress bar
        progressbar = tqdm(iterator, desc="Processing...")
        result_list = []

        for i in progressbar:
            progressbar.set_description(f"Processing {i}...")

            try:
                result = function(i, *args, **kwargs)
                result_list.append(result)
            except (SwaggerMappingError, ValueError) as e:
                errors.append(f"Error for {i}: {e}")

        # Remove None values
        result_list = [result_df for result_df in result_list if result_df is not None]

        # Log failed requests
        if len(errors) > 0:
            _logger.error(
                f"There was (were) {len(errors)}/{len(iterator)} failed request(s).\n"
                f"Show error messages (up to 5 messages only):"
            )
            _logger.error("\n".join([e for e in errors]))

        # If request returned any results, return as DataFrame, else raise SwaggerMappingError
        if len(result_list) > 0:
            result_df = pd.concat(result_list)
            result_df.reset_index(drop=True, inplace=True)
            return result_df

        else:
            raise SwaggerMappingError(f"Input values yield no results.")


class KinasesProvider(BaseProvider):
    """
    Class for kinases requests.

    Methods
    -------
    all_kinase_groups()
        Get all available kinase groups.
    all_kinase_families(group=None)
        Get all available kinase groups.
    all_kinases(groups=None, families=None, species=None)
        Get all available kinase names (optional: select kinase group, family and/or species).
    by_kinase_klifs_id(kinase_klifs_ids)
        Get kinases by one or more kinase KLIFS IDs.
    by_kinase_name(kinase_names)
        Get kinases by one or more kinase names (KLIFS or HGNC name).

    Notes
    -----
    Class methods all return a pandas.DataFrame of kinases (rows) with the (or a subset of the)
    following attributes (columns):

    Remote only:

        kinase.klifs_name : str
            Kinase name according to KLIFS.
        kinase.subfamily : str
            Kinase class.
            Available remotely only.
        kinase.full_name : str
            Full kinase name.
            Available remotely only.
        kinase.uniprot : str
            UniProt ID.
            Available remotely only.
        kinase.iuphar : int
            IUPHAR ID.
            Available remotely only.

    Local only:

        -

    Both local and remote:

        kinase.klifs_id : int
            Kinase KLIFS ID.
        kinase.gene_name : str
            Kinase name according to the HUGO Gene Nomenclature Committee.
            Available remotely only.
        kinase.family : str
            Kinase family.
        kinase.group : str
            Kinase group.
        species.klifs : str
            Species (KLIFS notation).
        kinase.pocket : str
            One-letter amino acid sequence for the kinase's 85 residue KLIFS pocket (gaps "-").
    """

    def all_kinase_groups(self):
        """
        Get all available kinase groups.

        Returns
        -------
        pandas.DataFrame
            Kinase groups (rows) with the following column: "kinase.group". Check class docstring
            for more information on columns.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_kinase_families(self, group=None):
        """
        Get all available kinase families (optional: select kinase group).

        Parameters
        ----------
        group : None or str
            Kinase group name (default is None, i.e. all kinase groups are selected).

        Returns
        -------
        pandas.DataFrame
            Kinase families (rows) with the following column: "kinase.family". Check class
            docstring for more information on columns.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If group does not exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_kinases(self, group=None, family=None, species=None):
        """
        Get all available kinase names (optional: select kinase group, family and/or species).

        Parameters
        ----------
        group : None or str
            Kinase group name (default is None, i.e. all kinase groups are selected).
        family : None or str
            Kinase family name (default is None, i.e. all kinase families are selected).
        species : None or str
            Species name (default is None, i.e. all species are selected).

        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with the following columns: "kinase.klifs_id", "kinase.hgnc_name",
            "kinase.full_name", "species.klifs". Check class docstring for more information on
            columns.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If group or family or species do not exist
            # TODO in the future: use ValueError instead but keep the original message
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_klifs_id(self, kinase_klifs_ids):  # TODO *kinases_ids
        """
        Get kinases by one or more kinase KLIFS IDs.

        Parameters
        ----------
        kinase_klifs_ids : int or list of int
            Kinase KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_name(self, kinase_names, species=None):
        """
        Get kinases by one or more kinase names (KLIFS or HGNC name).

        Parameters
        ----------
        kinase_names : str or list of str
            Kinase names (remote: KLIFS or HGNC name; local: any of the given kinase names).
        species : None or str
            Species name (default is None, i.e. all species are selected).


        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase names exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class LigandsProvider(BaseProvider):
    """
    Class for ligands requests.

    Methods
    -------
    all_ligands()
        Get all available ligands.
    by_kinase_klifs_id(kinase_klifs_ids)
        Get ligands by one or more kinase KLIFS IDs.
    by_kinase_name(kinase_names)
        Get ligands by one or more kinase names (KLIFS or HGNC name).
    by_ligand_klifs_id(ligand_klifs_ids)
        Get ligands by one or more ligand KLIFS IDs.
    by_ligand_expo_id(ligand_expo_ids)
        Get ligands by one or more Ligand Expo IDs (3-letter codes), i.e. the chemical component
        identifiers as defined by Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.

    Notes
    -----
    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the)
    following attributes (columns):

    Remote only:

        ligand.klifs_id : int
            Ligand KLIFS ID.
            Available remotely only.
        ligand.smiles : str
            Ligand SMILES.
            Available remotely only.
        ligand.inchikey : str
            Ligand InChI key.
            Available remotely only.

    Local only:

        -

    Both local and remote:

        ligand.expo_id : str
            Ligand PDB name.
        ligand.name : str
            Ligand name.
    """

    def all_ligands(self):
        """
        Get all available ligands.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        """
        Get ligands by one or more kinase KLIFS IDs.

        Parameters
        ----------
        kinase_klifs_ids : int or list of int
            Kinase KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_name(self, kinase_names):
        """
        Get ligands by one or more kinase names (KLIFS or HGNC name).

        Parameters
        ----------
        kinase_names : str or list of str
            Kinase names (remote: KLIFS or HGNC name; local: any of the given kinase names).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase names exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        """
        Get ligands by one or more ligand KLIFS IDs.

        Parameters
        ----------
        ligand_klifs_ids : int or list of int
            Ligand KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_expo_id(self, ligand_expo_ids):
        """
        Get ligands by one or more Ligand Expo IDs (3-letter codes), i.e. the chemical component
        identifiers as defined by Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.

        Parameters
        ----------
        ligand_expo_id : str or list of str
            Ligand Expo ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class StructuresProvider(BaseProvider):
    """
    Class for structures requests.

    Methods
    -------
    all_structures()
        Get all available structures.
    by_structure_klifs_id(structure_klifs_ids)
        Get structures by one or more structure KLIFS IDs.
    by_ligand_klifs_id(ligand_klifs_ids)
        Get structures by one or more ligand KLIFS IDs.
    by_kinase_klifs_id(kinase_klifs_ids)
        Get structures by one or more kinase KLIFS IDs.
    by_structure_pdb_id(structure_pdb_ids)
        Get structures by one or more structure PDB IDs.
    by_ligand_expo_id(ligand_expo_ids)
        Get structures by one or more Ligand Expo IDs (3-letter codes), i.e. the chemical component
        identifiers as defined by Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.
    by_kinase_name(kinase_names)
        Get structures by one or more kinase names (KLIFS or HGNC name).

    Notes
    -----
    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the)
    following attributes (columns):

    structure.klifs_id : int
        Structure KLIFS ID.
    structure.pdb_id : str
        Structure PDB ID.
    structure.alternate_model : str
        Alternate model. "-" if no alternate model.
    structure.chain : str
        Chain.
    species.klifs : str
        Species (KLIFS notation).
    kinase.klifs_id : int
        Kinase KLIFS ID.
    kinase.klifs_name : str
        Kinase name according to KLIFS.
    kinase.name_full : str
        Full kinase name.
        Available remotely when querying all kinases only.
    kinase.family : str
        Kinase family.
        Available locally only.
    kinase.group : str
        Kinase group.
        Available locally only.
    structure.pocket : str
        One-letter amino acid sequence for the structure's 85 residue KLIFS pocket (gaps "-").
    ligand.expo_id : str or int (0)
        Orthosteric ligand Ligand Expo ID. None if no ligand.
        Ligand Expo IDs (3-letter codes) are chemical component identifiers as defined by
        Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.
    ligand_allosteric.expo_id : str or
        Allosteric ligand Ligand Expo ID. None if no ligand.
        Ligand Expo IDs (3-letter codes) are chemical component identifiers as defined by
        Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.
    ligand.name : str
        Orthosteric ligand name. None if no ligand.
        Available locally only.
    ligand_allosteric.name : str
        Allosteric ligand name. None if no ligand.
        Available locally only.
    structure.dfg : str
        DFG conformation (in, out, out-like, na).
    structure.ac_helix : str
        aC helix conformation (in, out, out-like, na).
    structure.resolution : float
        Structure resolution in Angström.
    structure.qualityscore : float
        KLIFS quality score (considering missing residues/atoms and RMSD1/RMSD2).
    structure.missing_residues : int
        Number of missing residues.
    structure.missing_atoms : int
        Number of missing atoms.
    structure.rmsd1 : float
        RMSD between structure and reference structures based on full kinase domain.
    structure.rmsd2 : float
        RMSD between structure and reference structures based on kinase pocket residues.
    structure.front : bool
        Orthosteric ligand occupies KLIFS front cleft.
        Available remotely only.
    structure.gate : bool
        Orthosteric ligand occupies gate area.
        Available remotely only.
    structure.back : bool
        Orthosteric ligand occupies KLIFS back cleft.
        Available remotely only.
    structure.fp_i : bool
        Orthosteric ligand occupies KLIFS front pocket I (FP-I).
    structure.fp_ii : bool
        Orthosteric ligand occupies KLIFS front pocket II (FP-II).
    structure.bp_i_a : bool
        Orthosteric ligand occupies KLIFS back pocket I-A (BP-I-A).
    structure.bp_i_b : bool
        Orthosteric ligand occupies KLIFS back pocket I-B (BP-I-B).
    structure.bp_ii_in : bool
        Orthosteric ligand occupies KLIFS back pocket II DFG-in (BP-II-in).
    structure.bp_ii_a_in : bool
        Orthosteric ligand occupies KLIFS back pocket II-A DFG-in (BP-II-A-in).
    structure.bp_ii_b_in : bool
        Orthosteric ligand occupies KLIFS back pocket II-B DFG-in (BP-II-B-in).
    structure.bp_ii_out : bool
        Orthosteric ligand occupies KLIFS back pocket II DFG-out (BP-II-out).
    structure.bp_ii_b : bool
        Orthosteric ligand occupies KLIFS back pocket II-B (BP-II-B).
    structure.bp_iii : bool
        Orthosteric ligand occupies KLIFS back pocket III (BP-III).
    structure.bp_iv : bool
        Orthosteric ligand occupies KLIFS back pocket IV (BP-IV).
    structure.bp_v : bool
        Orthosteric ligand occupies KLIFS back pocket V (BP-V).
    structure.grich_distance : float
        G-rich loop conformation: Distance between the catalytic loop and the G-rich loop.
        Available remotely only.
    structure.grich_angle : float
        G-rich loop conformation: Angle of loop compared to hinge reagion and catalytic loop.
        Available remotely only.
    structure.grich_rotation : float
        G-rich loop conformation: Rotation of the G-rich loop compared to the catalytic loop.
        Available remotely only.
    structure.filepath : str
        Path to folder with the structure's coordinate files.
        Available locally only.
    """

    def all_structures(self):
        """
        Get all available structures.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_structure_klifs_id(self, structure_klifs_ids):
        """
        Get structures by one or more structure KLIFS IDs.

        Parameters
        ----------
        structure_klifs_ids : int or list of int
            Structure KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the structure KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        """
        Get structures by one or more ligand KLIFS IDs.

        Parameters
        ----------
        ligand_klifs_ids : int or list of int
            Ligand KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the ligand KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        """
        Get structures by one or more kinase KLIFS IDs.

        Parameters
        ----------
        kinase_klifs_ids : int or list of int
            Kinase KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_structure_pdb_id(
        self, structure_pdb_ids, structure_alternate_model=None, structure_chain=None
    ):
        """
        Get structures by one or more structure PDB IDs.

        Parameters
        ----------
        structure_pdb_ids : str or list of str
            Structure PDB ID(s).
        structure_alternate_model : None or str
            Alternate model (only used if only one structure PDB ID is given).
        structure_chain : None or str
            Chain (only used if only one structure PDB ID is given).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the structure PDB IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    @staticmethod
    def _filter_pdb_by_alt_chain(structures, structure_alternate_model=None, structure_chain=None):
        """
        If only one structure PDB ID given, filter structures by alternate model and/or chain.

        Parameters
        ----------
        structures : pandas.DataFrame
            Structures for one structures PDB ID (rows) with columns as described in the class docstring.
        structure_alternate_model : None or str
            Alternate model (only used if only one structure PDB ID is given).
        structure_chain : None or str
            Chain (only used if only one structure PDB ID is given).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """

        # Filter by alternate model if given
        if structure_alternate_model:
            structures = structures[
                structures["structure.alternate_model"] == structure_alternate_model
            ]
        # Filter by chain if given
        if structure_chain:
            structures = structures[structures["structure.chain"] == structure_chain]
        return structures

    def by_ligand_expo_id(self, ligand_expo_ids):
        """
        Get structures by one or more Ligand Expo IDs (3-letter codes), i.e. the chemical component
        identifiers as defined by Ligand Expo (http://ligand-expo.rcsb.org/) and used in the PDB.

        Parameters
        ----------
        ligand_expo_ids : str or list of str
            Ligand Expo ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_name(self, kinase_names):
        """
        Get structures by one or more kinase names (KLIFS or HGNC name).

        Parameters
        ----------
        kinase_names : str or list of str
            Kinase names (remote: KLIFS or HGNC name; local: any of the given kinase names).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class BioactivitiesProvider(BaseProvider):
    """
    Class for bioactivities requests (ChEMBL data linked to KLIFS data).

    Methods
    -------
    all_bioactivities()
        Get all available bioactivities.
    by_kinase_klifs_id(kinase_klifs_ids)
        Get bioactivities by one or more kinase KLIFS IDs.
    by_ligand_klifs_id(ligand_klifs_ids)
        Get bioactivities by one or more ligand KLIFS IDs.

    Notes
    -----
    Class methods all return a pandas.DataFrame of bioactivities (rows) with the (or a subset of
    the) following attributes (columns):
    kinase.pref_name : str
        Target name (ChEMBL notation).
    kinase.uniprot : str
        UniProt ID.
    species.chembl : str
        Species (ChEMBL notation).
    ligand.bioactivity_standard_value : str
        Standard bioactivity type (observation model) from ChEMBL.
    ligand.bioactivity_standard_relation : str
        Standard bioactivity relation from ChEMBL.
    ligand.bioactivity_standard_units : str
        Standard bioactivity unit from ChEMBL.
    ligand.bioactivity_pchembl_value : str
        pChEMBL value from ChEMBL: -Log(molar IC50, XC50, EC50, AC50, Ki, Kd or Potency).
    """

    def all_bioactivities(self, _top_n=None):
        """
        Get all available bioactivities.

        Parameters
        ----------
        _top_n: None or int
            If None, bioactivities for all ligands are returned (takes a few minutes).
            If set to n, bioactivities for the top n ligands are returned. This parameter
            is used for testing.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        """
        Get bioactivities by one or more kinase KLIFS IDs.

        Parameters
        ----------
        kinase_klifs_ids : int or list of int
            Kinase KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        """
        Get bioactivities by one or more ligand KLIFS IDs.

        Parameters
        ----------
        ligand_klifs_ids : int or list of int
            Ligand KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the ligand KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_expo_id(self, ligand_expo_id):
        """
        Get bioactivities by one or more ligand Expo IDs.

        Parameters
        ----------
        ligand_expo_id : str or list of str
            Ligand Expo ID(s).

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the ligand Expo IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class InteractionsProvider(BaseProvider):
    """
    Class for interactions requests.

    Methods
    -------
    interaction_types()
        Get all available interaction types.
    all_interactions()
        Get all available interaction fingerprints.
    by_structure_klifs_id(structure_klifs_ids)
        Get interactions by one or more structure KLIFS IDs.
    by_ligand_klifs_id(ligand_klifs_ids)
        Get interactions by one or more ligand KLIFS IDs.
    by_kinase_klifs_id(kinase_klifs_ids)
        Get interactions by one or more kinase KLIFS IDs.

    Notes
    -----
    Class methods all return a pandas.DataFrame of interactions (rows) with the (or a subset of
    the) following attributes (columns):

    structure.klifs_id : int
        Structure KLIFS ID.
    interaction.fingerprint : str
        Interaction fingerprint, a string of 595 zeros and ones for
        85 (pocket residues) * 7 (interaction features).
    interaction.id : int
        Interaction ID (1-7).
    interaction.name : str
        Interaction name (KLIFS notation).
    """

    @property
    def interaction_types(self):
        """
        Get all available interaction types.

        Returns
        -------
        pandas.DataFrame
            7 interaction types (rows) with the following columns:
            "interaction.id", "interaction.name".
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_interactions(self):
        """
        Get all available interaction fingerprints.

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns:
            "structure.klifs_id", "interaction.fingerprint".
            Check class docstring for more information on columns.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_structure_klifs_id(self, structure_klifs_ids):
        """
        Get interactions by one or more structure KLIFS IDs.

        Parameters
        ----------
        structure_klifs_ids : int or list of int
            Structure KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns:
            "structure.klifs_id", "interaction.fingerprint".
            Check class docstring for more information on columns.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the structure KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_ligand_klifs_id(self, ligand_klifs_ids):
        """
        Get interactions by one or more ligand KLIFS IDs.

        Parameters
        ----------
        ligand_klifs_ids : int or list of int
            Ligand KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns:
            "structure.klifs_id", "interaction.fingerprint".
            Check class docstring for more information on columns.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def by_kinase_klifs_id(self, kinase_klifs_ids):
        """
        Get interactions by one or more kinase KLIFS IDs.

        Parameters
        ----------
        kinase_klifs_ids : int or list of int
            Kinase KLIFS ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns:
            "structure.klifs_id", "interaction.fingerprint".
            Check class docstring for more information on columns.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase KLIFS IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class PocketsProvider(BaseProvider):
    """
    Class for pocket requests.
    Get residue ID and the corresponding residue KLIFS ID (plus kinase region label for each residue).

    Methods
    -------
    by_structure_klifs_id()
        Get a structure's residue ID and KLIFS ID by structure ID
        (plus kinase region label for each residue).

    Notes
    -----
    Class methods all return a pandas.DataFrame of interactions (rows) with the (or a subset of
    the) following attributes (columns):

    residue.klifs_id : int
        Residue KLIFS ID.
    residue.id : str
        Residue ID.
    residue.klifs_region_id : str
        KLIFS region assigned to pocket residue.
    residue.klifs_region : str
        KLIFS region name assigned to pocket residue.
    residue.klifs_color : str
        KLIFS color assigned to pocket residue.
    """

    def by_structure_klifs_id(self, structure_klifs_id):
        """
        Get a structure's residue ID and KLIFS ID by structure ID
        (plus kinase region label for each residue).

        Parameters
        ----------
        structure_klifs_id : str
            Structure KLIFS ID.

        Returns
        -------
        pandas.DataFrame
            Residue numbering details.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: Structure KLIFS ID does not exist.
        """
        raise NotImplementedError("Implement in your subclass!")

    @staticmethod
    def _add_klifs_region_details(pocket):
        """
        Add KLIFS region name and color as additional columns to the pocket DataFrame.

        Parameters
        ----------
        pandas.DataFrame
            Pocket data.

        Returns
        -------
        pandas.DataFrame
            Pocket data including KLIFS region names and colors.
        """

        pocket["residue.klifs_region"] = pocket["residue.klifs_region_id"].apply(
            lambda x: ".".join(x.split(".")[:-1])
        )
        pocket["residue.klifs_color"] = pocket["residue.klifs_region"].apply(
            lambda x: POCKET_KLIFS_REGION_COLORS[x]
        )

        return pocket


class CoordinatesProvider(BaseProvider):
    """
    Class for coordinates requests:
    Get coordinates for different entities (complex, ligand, pocket, protein, water) and
    input formats (mol2, pdb)
    in the form of different output formats (text, biopandas, rdkit).

    Attributes
    ----------
    options : dict of str: list of str
        Allowed input parameter options.

    Methods
    -------
    to_text(structure_klifs_id, entity, extension)
        Load coordinates as text.
    to_dataframe(structure_klifs_id, entity, extension)
        Load coordinates as DataFrame.
    to_rdkit(structure_klifs_id, entity, extension, compute2d)
        Load coordinates as RDKit molecule.

    Notes
    -----
    Class methods all return a pandas.DataFrame of atom coordinates (rows) with the (or a subset of
    the) following attributes (columns):

    atom.id : int
        Atom ID.
    atom.name : str
        Atom name.
    atom.x : float
        Atom x coordinate.
    atom.y : float
        Atom y coordinate.
    atom.z : float
        Atom z coordinate.
    atom.type : float
        Atom type. TODO from where?
    residue.subst_id : str
        Residue's substructure ID. TODO from where?
    residue.subst_name : str
        Residue's substructure name. TODO from where?
    atom.charge : float
        Atom charge.
    residue.name : float
        Residue name.
    residue.id : float
        Residue ID.
    residue.klifs_id : int
        Residue KLIFS ID (pocket residues only, other NaN).
    residue.klifs_region : float
        KLIFS region that the residue belongs to (pocket residues only, other NaN).
    """

    def __init__(self, *args, **kwargs):

        self.options = {
            "entities": ["complex", "ligand", "pocket", "protein", "water"],
            "extensions": ["mol2", "pdb"],
            "output_formats": ["text", "biopandas", "rdkit"],
        }

    def to_text(self, structure_klifs_id, entity="complex", extension="mol2"):
        """
        Get structural data content from KLIFS database as string (text).

        Parameters
        ----------
        structure_klifs_id : int
            Structure KLIFS ID.
            In the local module, this parameter is called "structure_klifs_id_or_filepath".
        structure_klifs_id_or_filepath : int / str or pathlib.Path
            Structure KLIFS ID or path to file.
            In the remote module, this parameter is called "structure_klifs_id".
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.
        extension : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.

        Returns
        -------
        str
            Structural data.
        """
        raise NotImplementedError("Implement in your subclass!")

    def to_dataframe(self, structure_klifs_id, entity="complex", extension="mol2"):
        """
        Get structural data content as DataFrame.

        Parameters
        ----------
        structure_klifs_id : int
            Structure KLIFS ID.
            In the local module, this parameter is called "structure_klifs_id_or_filepath".
        structure_klifs_id_or_filepath : int / str or pathlib.Path
            Structure KLIFS ID or path to file.
            In the remote module, this parameter is called "structure_klifs_id".
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.
        extension : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: Structure KLIFS ID does not exist.
        ValueError
            If input yields not result.
        """
        raise NotImplementedError("Implement in your subclass!")

    def to_rdkit(self, structure_klifs_id, entity="complex", extension="mol2", compute2d=True):
        """
        Get structural data content as RDKit molecule.

        Parameters
        ----------
        structure_klifs_id : int
            Structure KLIFS ID.
            In the local module, this parameter is called "structure_klifs_id_or_filepath".
        structure_klifs_id_or_filepath : int / str or pathlib.Path
            Structure KLIFS ID or path to file.
            In the remote module, this parameter is called "structure_klifs_id".
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.
        extension : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb.
            In the local module, when a filepath is passed to "structure_klifs_id_or_filepath", this
            parameter will be ignored and inferred from the filepath instead.
        compute2d : bool
            For entity=ligand only. Compute 2D coordinates (default) or keep 3D coordinates.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: Structure KLIFS ID does not exist.
        ValueError
            If input yields not result.
        """
        raise NotImplementedError("Implement in your subclass!")

    @staticmethod
    def _raise_invalid_extension(extension):
        """
        Check if extension is valid.

        Parameters
        ----------
        extension : str
            Extension name.

        Returns
        -------
        bool
            Extension is valid or invalid.
        """

        extensions = ["pdb", "mol2"]
        if extension not in extensions:
            raise ValueError(f"Invalid extension. Select from: {', '.join(extensions)}")
