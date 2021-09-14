"""
opencadd.databases.klifs.fields

Defines the fields available in KLIFS (remote and local) and their mapping to the names used
in opencadd.
"""

import pandas as pd


class Fields:
    """
    Class for KLIFS fields (remote and local) and their mapping to the names used in opencadd.

    Attributes
    ----------
    df : pandas.DataFrame
        Contains all field names in KLIFS (remote and local) and the names/dtypes used in opencadd.
        Column names:
        - field_types: KLIFS models as named in opencadd:
          - kinase_groups
          - kinase_families
          - kinases_all
          - kinases
          - ligands
          - structures
          - bioactivities
          - interactions
          - interaction_types
          - pockets
          - coordinates
        - opencadd.df_name: Field name as used in opencadd
        - opencadd.df_type: Field dtype as used in opencadd
        - klifs.remote: Field names as used in KLIFS remote
        - klifs.local_export: Field names as used in the KLIFS download file `KLIFS_export.csv`
        - klifs.local_overview: Field names as use din the KLIFS download file `overview.csv`

    Notes
    -----
    - "kinase_names": Kinase names: "kinase.gene_name (kinase.klifs_name)"
    - "kinase.klifs_name": Depending on availability: Manning name or UniProt gene name
    - "kinase.full_name": Depending on availability: HGNC gene name or Manning name or UniProt gene name
    - "kinase.gene_name": HGNC or MGI name
    - "kinase.uniprot": UniProt accession
    """

    def __init__(self, fields_path):
        self.df = pd.read_csv(fields_path)

    def _select_field_type(self, field_type):
        """
        Subset DataFrame by field type.

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.

        Returns
        -------
        pandas.DataFrame
            Subset DataFrame with selected fields.
        """
        df = self.df.copy()
        return df.groupby("field_type").get_group(field_type)

    def _to_dict(self, field_type, key_column_name, value_column_name):
        """
        Select two DataFrame columns as keys and values to generate a dictionary.

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.
        key_column_name : str
            Column name for the column to be used as keys.
        value_column_name : str
            Column name for the column to be used as values.

        Returns
        -------
        dict
            Selected columns formatted as dict.
        """
        df = self._select_field_type(field_type)

        # Select 2 columns and drop all rows with any NaN values
        df = df[[key_column_name, value_column_name]].dropna(how="any")

        # Cast DataFrame to dict
        if len(df) == 1:
            # If DataFrame has only one row, use ugly hack
            return {df.values[0][0]: df.values[0][1]}
        else:
            return df.set_index(key_column_name).squeeze().to_dict()

    def remote_to_oc_names(self, field_type):
        """
        Get a KLIFS remote to opencadd names mapping as dict.

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.

        Returns
        -------
        dict
            KLIFS remote to opencadd names mapping.
        """
        dict_ = self._to_dict(field_type, "klifs.remote", "opencadd.df_name")
        return dict_

    def local_export_to_oc_name(self, field_type="structures"):
        """
        Get a KLIFS local (`KLIFS_export.csv` download) to opencadd names mapping as dict.

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.

        Returns
        -------
        dict
            KLIFS local (`KLIFS_export.csv` download) to opencadd names mapping.
        """
        dict_ = self._to_dict(field_type, "klifs.local_export", "opencadd.df_name")
        return dict_

    def local_overview_to_oc_name(self, field_type="structures"):
        """
        Get a KLIFS local (`overview.csv` download) to opencadd names mapping as dict.

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.

        Returns
        -------
        dict
            KLIFS local (`overview.csv` download) to opencadd names mapping.
        """
        dict_ = self._to_dict(field_type, "klifs.local_overview", "opencadd.df_name")
        return dict_

    def oc_name_to_type(self, field_type, additional_dict=None):
        """
        Get an opencadd name to dtype mapping as dict. Used to standardize the opencadd output
        DataFrames!

        Parameters
        ----------
        field_type : str
            Field type, check class docstring for possible types.
        additional_dict : dict or None.
            If set, add this dictionary to the default dictionary.

        Returns
        -------
        dict
            opencadd name to dtype mapping.
        """
        dict_ = self._to_dict(field_type, "opencadd.df_name", "opencadd.df_type")
        if additional_dict is not None:
            for key, value in additional_dict.items():
                dict_[key] = value
        return dict_
