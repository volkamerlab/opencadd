"""
This function converts:
- SMILES
- InChI
- SDF
"""

from rdkit import Chem

__all__ = [
    "convert_smiles_to_mol",
    "convert_inchi_to_mol",
    "convert_mol_to_smiles",
    "convert_mol_to_inchi",
    "convert_sdf_to_mol_array",
    "convert_mol_to_sdf",
]


def _check_transform_file_ending(fn, f_format=".sdf"):
    """Checks if there is sdf(or optional other) file ending. If not it
    adds a file ending.
    """
    if fn.endswith(f_format):
        return fn
    else:
        return fn + f_format


def convert_smiles_to_mol(smiles, *args, **kwargs):
    """Converts SMILES to mol.

    Parameters
    ---------
    smiles: str
        The SMILES string that has to be converted.

    Returns
    -------
    mol: rdkit.Chem.Mol
        The molecule genrated from the SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles, *args, **kwargs)
    return mol


def convert_inchi_to_mol(inchi, *args, **kwargs):
    """Converts InChI to mol.

    Parameters
    ---------
    inchi: str
        The InChI string that has to be converted.

    Returns
    -------
    mol: rdkit.Chem.Mol
        The molecule genrated from the InChI string.
    """
    mol = Chem.MolFromInchi(inchi, *args, **kwargs)
    return mol


def convert_mol_to_smiles(mol, *args, **kwargs):
    """Converts mol to SMILES.

    Parameters
    ---------
    mol: rdkit.Chem.Mol
        The mol string that has to be converted.

    Returns
    -------
    smiles: str
        The SMILES string genrated from the mol.
    """
    smiles = Chem.MolToSmiles(mol, *args, **kwargs)
    return smiles


def convert_mol_to_inchi(mol, *args, **kwargs):
    """Converts mol to InChI.

    Parameters
    ---------
    mol: rdkit.Chem.Mol
        The mol string that has to be converted.

    Returns
    -------
    inchi: str
        The InchI string genrated from the mol.
    """
    inchi = Chem.MolToInchi(mol, *args, **kwargs)
    return inchi


def convert_sdf_to_mol_array(fn):
    """Converts molecules stored in an file to mol.

    Parameters
    ---------
    fn: str
        The filename of the sdf-file containing the molecules.

    Returns
    -------
    mol_array: list of rdkit.Chem.Mol
        A list of all molecules contained in the input sdf-file.

    Notes
    -----
    Files are load from the data folder.
    """
    # #This version doesn't work
    # suppl = Chem.SDMolSupplier(data_path(_check_transform_file_ending(fn)))
    # mol_array = []
    # for mol in suppl:
    #     mol_array.append(mol)
    #     return mol_array
    suppl = Chem.SDMolSupplier(_check_transform_file_ending(fn))
    mol_array = []
    for mol in suppl:
        mol_array.append(mol)
        return mol_array


def convert_mol_to_sdf(mol_array, fn="unnamed_mol_file"):
    """Generates an sdf-file containing molecules.

    Parameters
    ---------
    mol_array: list of rdkit.Chem.Mol
        A list of mol strings.
    fn: str, optional
        The name of the to be generated file.

    Returns
    -------
    fn.sdf: sdf-file
        Returns a sdf-file with the dafault-name "unnamed_mol_file.sdf"
        if not customized

    Notes
    -----
    Files are saved to the data folder.
    """
    # #This version doesn't work
    # w = Chem.SDWriter(data_path(_check_transform_file_ending(fn)))
    # for mol in mol_array:
    #     w.write(mol)
    w = Chem.SDWriter(_check_transform_file_ending(fn))
    for mol in mol_array:
        w.write(mol)
