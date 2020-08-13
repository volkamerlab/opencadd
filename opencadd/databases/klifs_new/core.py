"""
core.py

Defines core classes and functions.
"""

import logging

import pandas as pd

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
    def _abc_to_dataframe(abc_object):
        """
        Transform ABC object into a DataFrame (needed for KLIFS API results).

        Parameters
        ----------
        abc_object : list of abc.IDList or abc.KinaseInformation or abc.ligandDetails or abc.structureDetails
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


class KinasesProvider(BaseProvider):
    """
    Class for kinases requests.

    Class methods all return a pandas.DataFrame of kinases (rows) with the (or a subset of the) following attributes (columns):
    kinase.id : int
        Kinase ID.
    kinase.name : str
        Kinase name according to KLIFS.
    kinase.hgnc : str
        Kinase name according to the HUGO Gene Nomenclature Committee.
    kinase.family : str
        Kinase family.
    kinase.group : str
        Kinase group.
    kinase.class : str
        Kinase class. TODO: Where does this come from?
    species : str
        Species.
    kinase.name_full : str
        Full kinase name.
    kinase.uniprot : str
        UniProt ID.
    kinase.iuphar : int
        IUPHAR ID.
    kinase.pocket : str
        One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps notated with -).
    """

    def __init__(self):
        super().__init__()

    def all_kinase_groups(self):
        """
        Get all kinase groups available.
        
        Returns
        -------
        pandas.DataFrame
            Kinase groups (rows) with the following column: group. 
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_kinase_families(self, group=None):
        """
        Get all kinase families available.

        Parameters
        ----------
        group : None or str
            Kinase group name (default is None, i.e. all kinase groups are selected).
        
        Returns
        -------
        pandas.DataFrame
            Kinase families (rows) with the following column: family.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_kinases(self, groups=None, families=None, species=None):
        """
        Get all kinase names available. 
        TODO Format input strings: species capitalized, groups (Other?) and families upper/lower?
        
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
            Kinases (rows) with the following columns: kinase_ID, name, full_name, species.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get kinases by one or more kinase IDs.

        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
        """
        Get kinases by one or more kinase names (KLIFS or HGNC name).

        Returns
        -------
        pandas.DataFrame
            Kinases (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")


class LigandsProvider(BaseProvider):
    """
    Class for ligands requests.

    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) following attributes (columns):
    ligand.id : int
        Ligand ID.
    ligand.pdb : str
        Ligand PDB name.
    ligand.name : str
        Ligand name.
    ligand.smiles : str
        Ligand SMILES.
    ligand.inchikey : str
        InChI key.
    """

    def __init__(self):
        super().__init__()

    def all_ligands(self):
        """
        Get all ligands available.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get ligands by one or more kinase IDs.

        Parameters
        ----------
        kinase_ids : int or list of int
            Kinase ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
        """
        Get ligands by one or more kinase names (KLIFS or HGNC name).

        Parameters
        ----------
        kinase_names : str or list of str
            Kinase names(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get ligands by one or more ligand IDs.

        Parameters
        ----------
        ligand_ids : int or list of int
            Ligand ID(s).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
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
        """
        raise NotImplementedError("Implement in your subclass!")


class StructuresProvider(BaseProvider):
    """
    Class for structures requests.

    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) following attributes (columns):
    structure.id : int
        Structure ID.
    kinase.name : str
        Kinase name according to KLIFS.
    species : str
        Species.
    kinase.id : int
        Kinase ID.
    structure.pdb : str
        Structure PDB ID.
    structure.alternate_model : str
        Alternate model. If no alternate model, empty string.
    structure.chain : str
        Chain.
    structure.rmsd1 : float
        RMSD
    structure.rmsd2 : float
        RMSD
    kinase.pocket : str
        One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps notated with -).
    structure.resolution : float
        Structure resolution in Angström.
    structure.qualityscore : float
        KLIFS quality score.
    structure.missing_residues : int
        Number of missing residues.
    structure.missing_atoms : int
        Number of missing atoms.
    ligand.pdb : str or int (0)
        Orthosteric ligand PDB ID. If no ligand, 0.
    ligand.pdb_allosteric : str or 
        Allosteric ligand PDB ID. If no ligand, 0.
    structure.DFG : str
        DFG conformation (in, out, out-like, na).
    structure.ac_helix : str
        aC helix conformation (in, out, out-like, na).
    structure.grich_distance : float
        G-rich loop conformation annotation: Distance between the catalytic loop and the G-rich loop.
    structure.grich_angle : float
        G-rich loop conformation annotation: Angle of loop compared to hinge reagion and catalytic loop.
    structure.grich_rotation : float
        G-rich loop conformation annotation: Rotation of the G-rich loop compared to the catalytic loop
    structure.front : bool
        Orthosteric ligand occupies KLIFS front cleft.
    structure.gate : bool
        Orthosteric ligand occupies gate area.
    structure.back : bool
        Orthosteric ligand occupies KLIFS back cleft.
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
        Get all structures available.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_structure_ids(self, structure_ids):
        """
        Get structures by one or more structure IDs.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get structures by one or more ligand IDs.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get structures by one or more kinase IDs.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_structure_pdbs(self, structure_pdbs):
        """
        Get structures by one or more structure PDB IDs.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_pdbs(self, ligand_pdbs):
        """
        Get structures by one or more ligand PDB IDs.

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
        """
        Get structures by one or more kinase names (KLIFS or HGNC name).

        Returns
        -------
        pandas.DataFrame
            Structures (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")


class BioactivitiesProvider(BaseProvider):
    """
    Class for bioactivities requests (ChEMBL data linked to KLIFS data).
    
    Class methods all return a pandas.DataFrame of bioactivities (rows) with the (or a subset of the) following attributes (columns):
    pref_name : str
        Target name from ChEMBL.
    accession : str
        UniProt ID.
    organism : str
        Organism from ChEMBL.
    standard_type : str
        Standard bioactivity type (observation model) from ChEMBL.
    standard_relation : str
        Standard bioactivity relation from ChEMBL.
    standard_units : str
        Standard bioactivity unit from ChEMBL.
    pchembl_value : str
        pChEMBL value from ChEMBL, defined as -Log(molar IC50, XC50, EC50, AC50, Ki, Kd or Potency).
    """

    def __init__(self):
        super().__init__()

    def all_bioactivities(self):
        """
        Get all bioactivities available.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get bioactivities by one or more kinase IDs.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get structures by one or more ligand IDs.

        Returns
        -------
        pandas.DataFrame
            Bioactivities (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")


class InteractionsProvider(BaseProvider):
    """
    Class for interactions requests.

    Class methods all return a pandas.DataFrame of interactions (rows) with the (or a subset of the) following attributes (columns):
    structure_ID : int
        Structure ID.
    IFP : str
        Interaction fingerpint, a string of 595 zeros and ones for 85 (pocket residues) * 7 (interaction features).
    """

    def __init__(self):
        super().__init__()

    @property
    def interaction_types(self):
        """
        Get all interaction types available in KLIFS.

        Returns
        -------
        pandas.DataFrame
            7 interaction types (rows) with the following columns:
            position : int TODO
                Position 1-7.
            name : str
                Interaction type name.
        """
        raise NotImplementedError("Implement in your subclass!")

    def all_interactions(self):
        """
        Get all interactions available.

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: structure_ID, IFP.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_structure_ids(self, structure_ids):
        """
        Get interactions by one or more structure IDs.

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: structure_ID, IFP.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get interactions by one or more ligand IDs.

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: structure_ID, IFP.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_ids(self, kinase_ids):
        """
        Get interactions by one or more kinase IDs.

        Returns
        -------
        pandas.DataFrame
            Interactions (rows) with the following columns: structure_ID, IFP.
            Check class docstring for more information on columns.
        """
        raise NotImplementedError("Implement in your subclass!")


class CoordinatesProvider(BaseProvider):
    """
    Class for coordinates requests: Get coordinates for different entities for different input formats in the form of different output formats.
    """

    def __init__(self):
        super().__init__()

    def from_structure_id(self, structure_id, entity, input_format, output_format):
        raise NotImplementedError("Implement in your subclass!")

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
            Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).

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
            Structural data
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
        output_format : str or None
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

