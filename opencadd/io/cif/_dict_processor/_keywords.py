import warnings
import polars as pl
import numpy as np

from .. import _df_tools as tools


CAT_KEY = ["frame_code_category", "frame_code_keyword"]
COL_CAT = pl.col("frame_code_category")
COL_KEY = pl.col("frame_code_keyword")



def check_item_names(df: pl.DataFrame, col_name: str = "name"):
    if df.df.select(["category", "keyword"]).is_duplicated().any():
        pass # TODO
    if tools.match_names():
        # TODO
        pass


def item_default(df: tools.CIFItemCategoryDataFrame) -> pl.DataFrame:

    if not df.has_keyword("value"):
        raise ValueError("DataFrame doesn't contain necessary columns.")
    if "name" in df.keywords:
        if not df.match_names(col_name="name", reduce=True):
            raise ValueError()
    duplicated = df.is_duplicated()
    if duplicated.any():
        warnings.warn(
            f"Found {duplicated.sum()} duplicated rows:\n"
            f"{df.filter(duplicated)}\n"
            f"Removing duplicates."
        )
        df = df.unique()
    return rename_frame_code(df).select(["category", "keyword", "value"]).sort(
        by=["category", "keyword"])



def item_enumeration(cifdf: tools.CIFItemCategoryDataFrame):
    col_names = cifdf.df.columns
    if np.isin(["value", "detail"], col_names, invert=True).any():
        raise ValueError("Required columns not found.")
    if "name" in col_names and not cifdf.match_names(col_name="name"):
        warnings.warn("Column 'name' doesn't match implicit names. Ignoring 'name'.")
    return cifdf.df.select(["category", "keyword", "value", "detail"])





def process_df_item(df: pl.DataFrame):
    if collections.Counter(df.columns) != collections.Counter(
        [
            'frame_code_category',
            'frame_code_keyword',
            'name',
            'category_id',
            'mandatory_code',
        ]
    ):
        raise ValueError("Incorrect set of column names.")
    with pl.StringCache():
        df_1 = df.with_columns(pl.col("mandatory_code").str.to_uppercase().cast(pl.Categorical))
        if not df_1["mandatory_code"].unique().is_in(["YES", "NO", "IMPLICIT", "IMPLICIT-ORDINAL"]).all():
            raise ValueError("Unexpected value for mandatory_code.")
    if match_names(df_1, col_name="name", reduce=True):
        df_2 = rename_frame_code(df_1)
    else:
        df_2 = break_down_data_name_column(
            df=df_1, col_name="name", cat_col_name="category", key_col_name="keyword"
        )
        if not match_cat_key(df_2, col_name_category="category", col_name_keyword="keyword"):
            raise ValueError()
    if not column_is_subset(df=df_2, main_col_name="category", sub_col_name="category_id", reduce=True):
        raise ValueError()
    df_3 = df_2.select(["category", "keyword", "mandatory_code"])
    duplicated = df_3.is_duplicated()
    if duplicated.any():
        warnings.warn(
            f"Found {duplicated.sum()} duplicated rows:\n"
            f"{df_3.filter(duplicated)}\n"
            f"Removing duplicates."
        )
    return df_3.unique().sort(by=["category", "keyword"])







def process_df_item_type(df: pl.DataFrame) -> pl.DataFrame:
    if "name" in df.columns:
        if not match_names(df, col_name="name", reduce=True):
            raise ValueError()
    return


def process_df_item_dependent(df: pl.DataFrame) -> pl.DataFrame:
    return