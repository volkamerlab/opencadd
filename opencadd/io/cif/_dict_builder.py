import collections
import warnings
import polars as pl


def break_down_data_name_column(
        df: pl.DataFrame,
        col_name: str = "name",
        cat_col_name: str = "category",
        key_col_name: str = "keyword",
):
    df_new = df.with_columns(
        pl.col(col_name).str.extract(r'_([^.]+)\.(\S+)', 1).alias(cat_col_name),
        pl.col(col_name).str.extract(r'_([^.]+)\.(\S+)', 2).alias(key_col_name)
    )
    if not ("_" + df_new["category"] + "." + df_new["keyword"] == df_new["name"]).all():
        raise ValueError()
    return df_new.select(pl.exclude("name"))


def match_names(df: pl.DataFrame, col_name: str = "name", reduce: bool = True):
    """
    Whether the column with given name `col_name` is equal to '_{frame_code_category}.{frame_code_keyword}'.
    This will ignore null values in the `col_name` column.
    """
    df_ = df.with_columns(
        ("_" + pl.col("frame_code_category") + "." + pl.col("frame_code_keyword")).alias("__name__")
    )
    return column_is_subset(df_, main_col_name="__name__", sub_col_name=col_name, reduce=reduce)


def column_is_subset(df: pl.DataFrame, main_col_name: str, sub_col_name: str, reduce: bool = True):
    per_row = df.select(
        (
            (pl.col(sub_col_name).is_null()) |
            (pl.col(sub_col_name) == pl.col(main_col_name))
        ).alias("mask")
    )["mask"]
    return per_row.all() if reduce else per_row


def match_cat_key(df, col_name_category: str = "category", col_name_keyword: str = "keyword"):

    unique_nonmatching_frame_codes = df.filter(
        (pl.col(col_name_category) != pl.col("frame_code_category")) |
        (pl.col(col_name_keyword) != pl.col("frame_code_keyword"))
    ).select(
        pl.col(["frame_code_category", "frame_code_keyword"])
    ).unique()

    unique_matching_frame_codes = df.join(
        unique_nonmatching_frame_codes,
        on=["frame_code_category", "frame_code_keyword"],
        how="inner"
    ).filter(
        (pl.col(col_name_category) == pl.col("frame_code_category")) &
        (pl.col(col_name_keyword) == pl.col("frame_code_keyword"))
    ).select(
        pl.col(["frame_code_category", "frame_code_keyword"])
    ).unique()

    return unique_matching_frame_codes.sort(["frame_code_category", "frame_code_keyword"]).frame_equal(
        unique_nonmatching_frame_codes.sort(["frame_code_category", "frame_code_keyword"])
    )


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


def process_df_item_default(df: pl.DataFrame) -> pl.DataFrame:
    if not {"frame_code_category", "frame_code_keyword", "value"}.issubset(df.columns):
        raise ValueError("DataFrame doesn't contain necessary columns.")
    if "name" in df.columns:
        if not match_names(df, col_name="name", reduce=True):
            raise ValueError()
    duplicated = df.is_duplicated()
    if duplicated.any():
        warnings.warn(
            f"Found {duplicated.sum()} duplicated rows:\n"
            f"{df.filter(duplicated)}\n"
            f"Removing duplicates."
        )
        df = df.unique()
    return rename_frame_code(df).select(["category", "keyword", "value"]).sort(by=["category", "keyword"])




def process_df_item_type(df: pl.DataFrame) -> pl.DataFrame:
    if "name" in df.columns:
        if not match_names(df, col_name="name", reduce=True):
            raise ValueError()


def rename_frame_code(df: pl.DataFrame) -> pl.DataFrame:
    return df.rename({"frame_code_category": "category", "frame_code_keyword": "keyword"})