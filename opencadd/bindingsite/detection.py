from . import pdb
import re  # for filtering floats from a list of strings
from pathlib import Path  # for handling local paths



def select_best_pocket(binding_site_df, selection_method, selection_criteria, ascending=False):
    """
    Select the best binding site from the table of all detected binding sites,
    either by sorting the binding sites based on a set of properties in the table,
    or by applying a function on the property values.

    Parameters
    ----------
    binding_site_df : pandas.DataFrame
        Binding site data retrieved from the DoGSiteScorer webserver.
    selection_method : str
        Selection method for selecting the best binding site.
        Either 'sorting' or 'function'.
    selection_criteria : str or list
        If 'selection_method' is 'sorting':
            List of one or more property names.
        If 'selection_method' is 'function':
            Any valid python syntax that generates a list-like object
            with the same length as the number of detected binding sites.
    ascending : bool
        Optional; default: False.
        If set to True, the binding site with the lowest value will be selected,
        otherwise, the binding site with the highest value is selected.

    Returns
    -------
    str
        Name of the selected binding site.
    """
    df = binding_site_df

    if selection_method == "sorting":
        sorted_df = df.sort_values(by=selection_criteria, ascending=ascending)
    elif selection_method == "function":
        df["function_score"] = eval(selection_criteria)
        sorted_df = df.sort_values(by="function_score", ascending=ascending)
    else:
        raise ValueError(f"Binding site selection method unknown: {selection_method}")

    selected_pocket_name = sorted_df.iloc[0].name
    return selected_pocket_name


def calculate_pocket_coordinates_from_pocket_pdb_file(filepath):
    """
    Calculate the coordinates of a binding site using the binding site's PDB file
    downloaded from DoGSiteScorer.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Local filepath of the binding site's PDB file.

    Returns
    -------
    dict of list of int
        Binding site coordinates in format:
        `{'center': [x, y, z], 'size': [x, y, z]}`
    """
    with open(Path(filepath).with_suffix(".pdb")) as f:
        pdb_file_text_content = f.read()
    pdb_file_df = pdb.load_pdb_file_as_dataframe(pdb_file_text_content)
    pocket_coordinates_data = pdb_file_df["OTHERS"].loc[5, "entry"]
    coordinates_data_as_list = pocket_coordinates_data.split()
    # select strings representing floats from a list of strings
    coordinates = [float(element) for element in coordinates_data_as_list if
                   re.compile(r'\d+(?:\.\d*)').match(element)]
    pocket_coordinates = {
        "center": coordinates[:3],
        "size": [coordinates[-1] * 2 for dim in range(3)],
    }
    return pocket_coordinates


def get_pocket_residues(pocket_pdb_filepath):
    """
    Get residue-IDs and names of a specified pocket.

    Parameters
    ----------
    pocket_pdb_filepath : str or pathlib.Path
        Path of pocket's PDB file.

    Returns
    -------
    pandas.DataFrame
        Table of residues names and IDs for the selected binding site.
    """

    with open(Path(pocket_pdb_filepath).with_suffix(".pdb")) as f:
        pdb_content = f.read()
    atom_info = pdb.load_pdb_file_as_dataframe(pdb_content)["ATOM"]
    # Drop duplicates, since the PDB file contains one entry per atom,
    # but we only need one entry per residue
    atom_info.sort_values("residue_number", inplace=True)
    atom_info.drop_duplicates(subset="residue_number", keep="first", inplace=True)
    atom_info.reset_index(drop=True, inplace=True)
    atom_info.index += 1
    return atom_info[["residue_number", "residue_name"]]
