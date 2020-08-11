"""
core.py

Defines core classes and functions.
"""

import logging

_logger = logging.getLogger(__name__)


class KinasesProvider:
    """
    Class for kinases requests.

    Class methods all return a pandas.DataFrame of kinases (rows) with the (or a subset of the) following attributes (columns):
    kinase_ID : int
        Kinase ID.
    name : str
        Kinase name according to KLIFS.
    HGNC : str
        Kinase name according to the HUGO Gene Nomenclature Committee.
    family : str
        Kinase family.
    group : str
        Kinase group.
    kinase_class : str
        Kinase class. TODO: Where does this come from?
    species : str
        Species.
    full_name : str
        Full kinase name.
    uniprot : str
        UniProt ID.
    iuphar : int
        IUPHAR ID.
    pocket : str
        One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps notated with -).
    """

    def __init__(self):
        """Empty init."""

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

    def _rename_remote_columns(kinases):
        """
        
        """


class LigandsProvider:
    """
    Class for ligands requests.

    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) following attributes (columns):
    ligand_ID : int
        Ligand ID.
    PDB-code : str
        Ligand PDB name.
    Name : str
        Ligand name.
    SMILES : str
        Ligand SMILES.
    InChIKey : str
        InChI key.
    """

    def __init__(self):
        """Empty init."""

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

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_ids(self, ligand_ids):
        """
        Get ligands by one or more ligand IDs.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_kinase_names(self, kinase_names):
        """
        Get ligands by one or more kinase names (KLIFS or HGNC name).

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")

    def from_ligand_pdbs(self, ligand_pdbs):
        """
        Get ligands by one or more ligand PDB IDs.

        Returns
        -------
        pandas.DataFrame
            Ligands (rows) with columns as described in the class docstring.
        """
        raise NotImplementedError("Implement in your subclass!")


class StructuresProvider:
    """
    Class for structures requests.

    Class methods all return a pandas.DataFrame of ligands (rows) with the (or a subset of the) following attributes (columns):
    structure_ID : int
        Structure ID.
    kinase : str
        Kinase name according to KLIFS.
    species : str
        Species.
    kinase_ID : int
        Kinase ID.
    pdb : str
        Complex PDB ID.
    alt : str
        Alternate model. If no alternate model, empty string.
    chain : str
        Chain.
    rmsd1 : float
        RMSD
    rmsd2 : float
        RMSD
    pocket : str
        One-letter amino acid sequence for the 85 residue KLIFS pocket (gaps notated with -).
    resolution : float
        Structure resolution in Angström.
    quality_score : float
        KLIFS quality score.
    missing_residues : int
        Number of missing residues.
    missing_atoms : int
        Number of missing atoms.
    ligand : str or int (0)
        Orthosteric ligand PDB ID. If no ligand, 0.
    allosteric_ligand : str or 
        Allosteric ligand PDB ID. If no ligand, 0.
    DFG : str
        DFG conformation (in, out, out-like, na).
    aC_helix : str
        aC helix conformation (in, out, out-like, na).
    Grich_distance : float
        G-rich loop conformation annotation: Distance between the catalytic loop and the G-rich loop.
    Grich_angle : float
        G-rich loop conformation annotation: Angle of loop compared to hinge reagion and catalytic loop.
    Grich_rotation : float
        G-rich loop conformation annotation: Rotation of the G-rich loop compared to the catalytic loop
    front : bool
        Orthosteric ligand occupies KLIFS front cleft.
    gate : bool
        Orthosteric ligand occupies gate area.
    back : bool
        Orthosteric ligand occupies KLIFS back cleft.
    fp_I : bool
        Orthosteric ligand occupies KLIFS front pocket I (FP-I).
    fp_II : bool
        Orthosteric ligand occupies KLIFS front pocket II (FP-II).
    bp_I_A : bool
        Orthosteric ligand occupies KLIFS back pocket I-A (BP-I-A).
    bp_I_B : bool
        Orthosteric ligand occupies KLIFS back pocket I-B (BP-I-B).
    bp_II_in : bool
        Orthosteric ligand occupies KLIFS back pocket II DFG-in (BP-II-in).
    bp_II_A_in : bool
        Orthosteric ligand occupies KLIFS back pocket II-A DFG-in (BP-II-A-in).
    bp_II_B_in : bool
        Orthosteric ligand occupies KLIFS back pocket II-B DFG-in (BP-II-B-in).
    bp_II_out : bool
        Orthosteric ligand occupies KLIFS back pocket II DFG-out (BP-II-out).
    bp_II_B : bool
        Orthosteric ligand occupies KLIFS back pocket II-B (BP-II-B).
    bp_III : bool
        Orthosteric ligand occupies KLIFS back pocket III (BP-III).
    bp_IV : bool
        Orthosteric ligand occupies KLIFS back pocket IV (BP-IV).
    bp_V : bool
        Orthosteric ligand occupies KLIFS back pocket V (BP-V).
    """

    def __init__(self):
        """Empty init."""

    @property
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


class BioactivitiesProvider:
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
        """Empty init."""

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


class InteractionsProvider:
    """
    Class for interactions requests.

    Class methods all return a pandas.DataFrame of interactions (rows) with the (or a subset of the) following attributes (columns):
    structure_ID : int
        Structure ID.
    IFP : str
        Interaction fingerpint, a string of 595 zeros and ones for 85 (pocket residues) * 7 (interaction features).
    """

    def __init__(self):
        """Empty init."""

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


class CoordinatesProvider:
    """
    Class for coordinates requests: Get coordinates for different entities for different input formats in the form of different output formats.
    """

    def __init__(self):
        """Empty init."""

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
