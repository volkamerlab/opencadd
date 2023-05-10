from typing import Dict
import polars as pl
import numpy as np

from .._filestruct import DDL2CIFFile
from .._df_tools import CIFItemCategoryDataFrame


mapping = {
    "data": {
        'category_group_list',
        'datablock',
        'dictionary',
        'dictionary_history',
        'item_type_list',
        'item_units_conversion',
        'item_units_list',
        'pdbx_comparison_operator_list',
        'pdbx_conditional_context_list',
        'pdbx_dictionary_component',
        'pdbx_dictionary_component_history',
        'pdbx_item_linked_group',
        'pdbx_item_linked_group_list',
        'sub_category'
    },
    "def_cat": {
        'category',
        'category_examples',
        'category_group',
        'category_key',
        'ndb_category_examples',
        'pdbx_category_conditional_context',
        'pdbx_category_context',
        'pdbx_category_description'
    },
    "def_key": {
        'item',
        'item_aliases',
        'item_default',
        'item_dependent',
        'item_description',
        'item_enumeration',
        'item_examples',
        'item_linked',
        'item_range',
        'item_related',
        'item_sub_category',
        'item_type',
        'item_type_conditions',
        'item_units',
        'pdbx_item',
        'pdbx_item_conditional_context',
        'pdbx_item_context',
        'pdbx_item_description',
        'pdbx_item_enumeration',
        'pdbx_item_enumeration_details',
        'pdbx_item_examples',
        'pdbx_item_range',
        'pdbx_item_type'
    }
}

process_func = {
    "data": None,
    "def_cat": None,
    "def_key": None
}


def process_def_key(df: pl.DataFrame):
    if np.isin(["category", "keyword"], df.columns).any():
        raise ValueError
    df = df.rename(
        {"frame_code_category": "category", "frame_code_keyword": "keyword"}
    )
    return CIFItemCategoryDataFrame(
        df=df,
        col_name_block_code=None,
        col_name_frame_code_category="category",
        col_name_frame_code_keyword="keyword",
    )



def start(dic: DDL2CIFFile):
    if dic.df.select(pl.col("data_name_keyword").is_in(["category", "keyword"]).any())[0, 0]:
        raise ValueError()
    if dic.count_data_blocks != 1:
        raise ValueError()
    dfs_dict = dic.df_per_category(part="all", reduce=True)
    for data_type, dfs_dict_type in dfs_dict.keys():
        if not isinstance(dfs_dict_type, dict):
            continue
        for cat_name, cat_df in dfs_dict_type.items():
            df_struct = process_func[data_type](cat_df)
            dfs_dict_type[cat_name] = mapping[data_type][cat_name](df_struct)
    return dfs_dict

