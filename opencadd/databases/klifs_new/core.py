"""
core.py

Defines core classes and functions.
"""

import logging

from bravado_core.exception import SwaggerMappingError
import pandas as pd
from tqdm import tqdm

_logger = logging.getLogger(__name__)

COORDINATES_PARAMETERS = {
    "entities": ["complex", "ligand", "pocket", "protein", "water"],
    "input_formats": ["mol2", "pdb"],
    "output_formats": ["text", "biopandas", "rdkit"],
}


class BaseProvider:
    """
    Base class for KLIFS requests (local and remote).
    """

    def __init__(self):
        """Empty init."""

    @staticmethod
    def _cast_to_list(value):
        """
        Cast value to list if it is not a list, else return as-is.
        """
        if not isinstance(value, list):
            return [value]
        else:
            return value

    @staticmethod
    def _abc_to_dataframe(abc_object):
        """
        Transform ABC object into a DataFrame (needed for KLIFS API results).

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

        return pd.DataFrame(results_dict)

    @staticmethod
    def _rename_dataframe_columns(dataframe, columns_mapping):
        """
        Rename DataFrame columns.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Remote query result.
        columns_mapping : dict
            Mapping of old to new column names.
        """

        # Rename columns
        dataframe.rename(columns=columns_mapping, inplace=True)

        return dataframe

    @staticmethod
    def _format_dataframe(dataframe, column_names):
        """
        Format a DataFrame: Sort columns, drop duplicate rows, and reset index.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Remote query result.
        column_names : list of str
            Column names in the order of interest for output.

        Raises
        ------
        ValueError
            If DataFrame is empty.
        """

        # Select and sort columns
        dataframe = dataframe[column_names].copy()

        # Drop duplicate rows
        dataframe.drop_duplicates(inplace=True)

        # Reset index
        dataframe.reset_index(drop=True, inplace=True)

        # Handy empty results
        if dataframe.shape[0] > 0:
            return dataframe
        else:
            raise ValueError(f"Input values yield no results.")

    @staticmethod
    def _multiple_remote_requests(
        function, iterator, additional_parameters=None,
    ):
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
        additional_parameters : None or list
            List of additional parameters (not iterated) that are needed for function.
            None if no additional parameters.
        """

        if not isinstance(iterator, list):
            iterator = [iterator]

        # If single request fails, error will be raised. Catch errors in list.
        errors = []

        # Get result for each iterator element
        # If more than 10 requests, print progress bar
        progressbar = tqdm(iterator, desc="Processing...")
        result_list = []

        for i in progressbar:
            progressbar.set_description(f"Processing {i}...")

            if additional_parameters is not None:
                try:
                    result = function(i, *additional_parameters)
                    result_list.append(result)
                except Exception as e:
                    errors.append(f"Error for {i}: {e}")

            else:
                try:
                    result = function(i)
                    result_list.append(result)
                except Exception as e:
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
            raise SwaggerMappingError(
                f"None of the input values exist, thus no results are returned."
            )


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
    from_kinase_ids(kinase_ids)
        Get kinases by one or more kinase IDs.
    from_kinase_names(kinase_names)
        Get kinases by one or more kinase names (KLIFS or HGNC name).

    Notes
    -----
    Class methods all return a pandas.DataFrame of kinases (rows) with the (or a subset of the) 
    following attributes (columns):

    Remote only:

        kinase.id : int
            Kinase ID.
        kinase.hgnc : str
            Kinase name according to the HUGO Gene Nomenclature Committee.
        kinase.class : str
            Kinase class.
        kinase.name_full : str
            Full kinase name.
        kinase.uniprot : str
            UniProt ID.
        kinase.iuphar : int
            IUPHAR ID.

    Local only:

        -

    Both local and remote:
        
        kinase.name : str
            Kinase name according to KLIFS.
        kinase.family : str
            Kinase family.
        kinase.group : str
            Kinase group.
        species.klifs : str
            Species (KLIFS notation).
        kinase.pocket : str
            One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps notated with -).
    """

    def __init__(self):

        super().__init__()

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
            Kinases (rows) with the following columns: "kinase.id", "kinase.name", "kinase.name_", 
            "species.klifs". Check class docstring for more information on columns.
        
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If group or family or species do not exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get kinases by one or more kinase IDs.

        Parameters
        ----------
        kinase_ids : int or list of int
            KLIFS kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with columns as described in the class docstring.
        
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
        """
        Get kinases by one or more kinase names (KLIFS or HGNC name).

        Parameters
        ----------
        kinase_names : str or list of str
            Kinase names (remote: KLIFS or HGNC name; local: any of the given kinase names).

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
    from_kinase_ids(kinase_ids)
        Get ligands by one or more kinase IDs.
    from_kinase_names(kinase_names)
        Get ligands by one or more kinase names (KLIFS or HGNC name).
    from_ligand_ids(ligand_ids)
        Get ligands by one or more ligand IDs.
    from_ligand_pdbs(ligand_pdbs)
        Get ligands by one or more ligand PDB IDs.

    Notes
    -----
    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) 
    following attributes (columns):

    Remote only:

        ligand.id : int
            Ligand ID.
        ligand.smiles : str
            Ligand SMILES.
        ligand.inchikey : str
            Ligand InChI key.

    Local only:

        -

    Both local and remote:
    
        ligand.pdb : str
            Ligand PDB name.
        ligand.name : str
            Ligand name.
    """

    def __init__(self):

        super().__init__()

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

    def from_kinase_ids(self, kinase_ids):
        """
        Get ligands by one or more kinase IDs.

        Parameters
        ----------
        kinase_ids : int or list of int
            KLIFS kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
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

    def from_ligand_ids(self, ligand_ids):
        """
        Get ligands by one or more ligand IDs.

        Parameters
        ----------
        ligand_ids : int or list of int
            KLIFS ligand ID(s).

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

    def from_ligand_pdbs(self, ligand_pdbs):
        """
        Get ligands by one or more ligand PDB IDs.

        Parameters
        ----------
        ligand.pdb : str or list of str
            Ligand PDB ID(s).

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
    from_structure_ids(structure_ids)
        Get structures by one or more structure IDs.
    from_ligand_ids(ligand_ids)
        Get structures by one or more ligand IDs.
    from_kinase_ids(kinase_ids)
        Get structures by one or more kinase IDs.
    from_structure_pdbs(structure_pdbs)
        Get structures by one or more structure PDB IDs.
    from_ligand_pdbs(ligand_pdbs)
        Get structures by one or more ligand PDB IDs.
    from_kinase_names(kinase_names)
        Get structures by one or more kinase names (KLIFS or HGNC name).

    Notes
    -----
    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) 
    following attributes (columns):

    Remote only:

        structure.id : int
            Structure ID.
        kinase.id : int
            Kinase ID.
        structure.front : bool
            Orthosteric ligand occupies KLIFS front cleft.
        structure.gate : bool
            Orthosteric ligand occupies gate area.
        structure.back : bool
            Orthosteric ligand occupies KLIFS back cleft.
        structure.grich_distance : float
            G-rich loop conformation: Distance between the catalytic loop and the G-rich loop.
        structure.grich_angle : float
            G-rich loop conformation: Angle of loop compared to hinge reagion and catalytic loop.
        structure.grich_rotation : float
            G-rich loop conformation: Rotation of the G-rich loop compared to the catalytic loop

    Local only:

        structure.filepath : str
            Path to folder with the structure's coordinate files.
        kinase.family : str
            Kinase family.
        kinase.group : str
            Kinase group.
        kinase.name_full : str
            Full kinase name.
        ligand.name : str
            Orthosteric ligand name. None if no ligand.
        ligand.name_allosteric : str
            Allosteric ligand name. None if no ligand.
    
    Both remote and local:

        kinase.name : str
            Kinase name according to KLIFS.
        species.klifs : str
            Species (KLIFS notation).
        structure.pdb : str
            Structure PDB ID.
        structure.alternate_model : str
            Alternate model. None if no alternate model.
        structure.chain : str
            Chain.
        structure.rmsd1 : float
            RMSD between structure and reference structures based on full kinase domain.
        structure.rmsd2 : float
            RMSD between structure and reference structures based on kinase pocket residues.
        kinase.pocket : str
            One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps "-").
        structure.resolution : float
            Structure resolution in Angström.
        structure.qualityscore : float
            KLIFS quality score (considering missing residues/atoms and RMSD1/RMSD2).
        structure.missing_residues : int
            Number of missing residues.
        structure.missing_atoms : int
            Number of missing atoms.
        ligand.pdb : str or int (0)
            Orthosteric ligand PDB ID. None if no ligand.
        ligand.pdb_allosteric : str or 
            Allosteric ligand PDB ID.  None if no ligand.
        structure.DFG : str
            DFG conformation (in, out, out-like, na).
        structure.ac_helix : str
            aC helix conformation (in, out, out-like, na).
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
    """

    def __init__(self):
        super().__init__()

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

    def from_structure_ids(self, structure_ids):
        """
        Get structures by one or more structure IDs.

        Parameters
        ----------
        structure_ids : int or list of int
            KLIFS structure ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the structure IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get structures by one or more ligand IDs.

        Parameters
        ----------
        ligand_ids : int or list of int
            KLIFS ligand ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the ligand IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get structures by one or more kinase IDs.

        Parameters
        ----------
        kinase_ids : int or list of int
            KLIFS kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_structure_pdbs(
        self, structure_pdbs, structure_alternate_model=None, structure_chain=None
    ):
        """
        Get structures by one or more structure PDB IDs.

        Parameters
        ----------
        structure_pdbs : str or list of str
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
    def _filter_pdb_by_alt_chain(
        structures, structure_alternate_model=None, structure_chain=None
    ):
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

    def from_ligand_pdbs(self, ligand_pdbs):
        """
        Get structures by one or more ligand PDB IDs.

        Parameters
        ----------
        ligand_pdbs : str or list of str
            Ligand PDB ID(s).

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

    def from_kinase_names(self, kinase_names):
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
    from_kinase_ids(kinase_ids)
        Get bioactivities by one or more kinase IDs.
    from_ligand_ids(ligand_ids)
        Get bioactivities by one or more ligand IDs.

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

    def __init__(self):
        super().__init__()

    def all_bioactivities(self, n=None):
        """
        Get all available bioactivities.

        Parameters
        ----------
        n : None or int
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

    def from_kinase_ids(self, kinase_ids):
        """
        Get bioactivities by one or more kinase IDs.
        
        Parameters
        ----------
        kinase_ids : int or list of int
            KLIFS kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get bioactivities by one or more ligand IDs.

        Parameters
        ----------
        ligand_ids : int or list of int
            KLIFS ligand ID(s).

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the ligand IDs exist.
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
    from_structure_ids(structure_ids)
        Get interactions by one or more structure IDs.
    from_ligand_ids(ligand_ids)
        Get interactions by one or more ligand IDs.
    from_kinase_ids(kinase_ids)
        Get interactions by one or more kinase IDs.
    
    Notes
    -----
    Class methods all return a pandas.DataFrame of interactions (rows) with the (or a subset of 
    the) following attributes (columns):

    structure.id : int
        Structure ID.
    interaction.fingerprint : str
        Interaction fingerprint, a string of 595 zeros and ones for 
        85 (pocket residues) * 7 (interaction features).
    interaction.id : int
        Interaction ID (1-7).
    interaction.name : str
        Interaction name (KLIFS notation).
    """

    def __init__(self):
        super().__init__()

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
            "structure.id", "interaction.fingerprint". 
            Check class docstring for more information on columns.
            
        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_structure_ids(self, structure_ids):
        """
        Get interactions by one or more structure IDs.

        Parameters
        ----------
        structure_ids : int or list of int
            KLIFS structure ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: 
            "structure.id", "interaction.fingerprint". 
            Check class docstring for more information on columns.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the structure IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get interactions by one or more ligand IDs.

        Parameters
        ----------
        ligand_ids : int or list of int
            KLIFS ligand ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: 
            "structure.id", "interaction.fingerprint". 
            Check class docstring for more information on columns.
            
        Raises
        ------
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get interactions by one or more kinase IDs.
        
        Parameters
        ----------
        kinase_ids : int or list of int
            KLIFS kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: 
            "structure.id", "interaction.fingerprint". 
            Check class docstring for more information on columns.
            
        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: If none of the kinase IDs exist.
        ValueError
            If DataFrame is empty.
        """
        raise NotImplementedError("Implement in your subclass!")


class PocketsProvider(BaseProvider):
    """
    Class for pocket requests.
    Get PDB and KLIFS pocket residues numbering (plus kinase region label for each residue).

    Methods
    -------
    from_structure_id()
        Get a structure's residue numbering in PDB and KLIFS by structure ID
        (plus kinase region label for each residue).
    """

    def __init__(self):
        super().__init__()

    def from_structure_id(self, structure_id):
        """
        Get a structure's residue numbering in PDB and KLIFS by structure ID
        (plus kinase region label for each residue).

        Parameters
        ----------
        structure_id : str
            KLIFS structure ID.

        Returns
        -------
        pandas.DataFrame
            Residue numbering details.

        Raises
        ------
        bravado_core.exception.SwaggerMappingError
            Remote module: Structure ID does not exist.
        """
        raise NotImplementedError("Implement in your subclass!")


class CoordinatesProvider(BaseProvider):
    """
    Class for coordinates requests: 
    Get coordinates for different entities (complex, ligand, pocket, protein, water) and 
    input formats (mol2, pdb)
    in the form of different output formats (text, biopandas, rdkit).
    """

    def __init__(self):
        super().__init__()

    @staticmethod
    def _fetch_text(structure_id, entity="complex", input_format="mol2"):
        """
        Get structural data content from KLIFS database as string (text).

        Parameters
        ----------
        structure_id : str
            KLIFS structure ID.
        entity : str
            Structural entity: complex (default), ligand, pocket, or protein.
        input_format : str
            Input file format (fetched from KLIFS): mol2 (default) or pdb (only if entity=complex).

        Returns
        -------
        str
            Structural data.
        """

        if entity == "complex" and input_format == "mol2":
            return (
                KLIFS_CLIENT.Structures.get_structure_get_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "complex" and input_format == "pdb":
            return (
                KLIFS_CLIENT.Structures.get_structure_get_pdb_complex(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "ligand" and input_format == "mol2":
            return (
                KLIFS_CLIENT.Structures.get_structure_get_ligand(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "pocket" and input_format == "mol2":
            return (
                KLIFS_CLIENT.Structures.get_structure_get_pocket(
                    structure_ID=structure_id
                )
                .response()
                .result
            )
        elif entity == "protein" and input_format == "mol2":
            return (
                KLIFS_CLIENT.Structures.get_structure_get_protein(
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
            Structural data.
        """

        ppdb = PandasPdb()

        pdb_dict = ppdb._construct_df(pdb_text.splitlines(True))

        print(f"Structural data keys: {pdb_dict.keys()}")

        return pdb_dict

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
            mol2_df = pmol.read_mol2(
                str(mol2_file),
                columns={
                    0: ("atom_id", int),
                    1: ("atom_name", str),
                    2: ("x", float),
                    3: ("y", float),
                    4: ("z", float),
                    5: ("atom_type", str),
                    6: ("subst_id", int),
                    7: ("subst_name", str),
                    8: ("charge", float),
                    9: ("backbone", str),
                },
            )

        except ValueError:
            mol2_df = pmol.read_mol2(str(mol2_file))

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

    @staticmethod
    def check_parameter_validity(entity, input_format, output_format=None):
        """
        Check if entity and input/output format (and their combinations) are valid.

        Parameters
        ----------
        entity : str
            Structural entity: complex, ligand, pocket, protein, or water (only in local module).
        input_format : str
            File input format: mol2 or pdb (only for entity=complex).
        output_format : None or str
            Output format: text (only in remote module), biopandas, or rdkit (only for entity=ligand).
        """

        # Check if parameters are valid
        if entity not in COORDINATES_PARAMETERS["entities"]:
            raise ValueError(
                f"Invalid entity. Select from {', '.join(COORDINATES_PARAMETERS['entities'])}."
            )
        if input_format not in COORDINATES_PARAMETERS["input_formats"]:
            raise ValueError(
                f"Invalid input format. Select from {', '.join(COORDINATES_PARAMETERS['input_formats'])}."
            )
        if output_format:
            if output_format not in COORDINATES_PARAMETERS["output_formats"]:
                raise ValueError(
                    f"Invalid output format. Select from {', '.join(COORDINATES_PARAMETERS['output_formats'])}."
                )

        # Check if parameter combination is valid
        if input_format == "pdb" and entity != "complex":
            raise ValueError(f"Entity {entity} is only available in mol2 format.")
        if output_format:
            if output_format == "rdkit" and entity != "ligand":
                raise ValueError(
                    f"Only entity ligand can be fetched as rdkit molecule."
                )
