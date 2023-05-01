from typing import Optional, Union, Sequence

import collections
import warnings

import numpy as np
import polars as pl

from opencadd._df_tools import PolarsDataFrame


class CIFDictItemDataFrame(PolarsDataFrame):

    def __init__(
            self,
            df: pl.DataFrame,
            col_name_cat: Optional[str] = "frame_code_category",
            col_name_key: Optional[str] = "frame_code_keyword",
    ):
        super().__init__(df)
        self._col_name_cat = col_name_cat
        self._col_name_key = col_name_key
        self._col_cat = pl.col(self._col_name_cat)
        self._col_key = pl.col(self._col_name_key)
        return

    def match_names(
            self,
            col_name: str = "name",
            reduce: bool = True
    ):
        """
        Whether the column with given name `col_name` is equal to '_{frame_code_category}.{frame_code_keyword}'.
        This will ignore null values in the `col_name` column.
        """
        df_ = self.df.with_columns(
            ("_" + self._col_cat + "." + self._col_key).alias("__name__")
        )
        col_is_sub = PolarsDataFrame(df_).column_is_subset(col_1="__name__", col_2=col_name)
        return (col_is_sub[1] is True or col_is_sub[1] == "__name__") if reduce else col_is_sub[2]["2_is_sub_1"]

    def rename_id_cols(self, cat_name, key_name):
        if self.has_keyword((cat_name, key_name)):
            raise ValueError()
        self._df = self._df.rename({self._col_name_cat: cat_name, self._col_name_key: key_name})
        self._col_name_cat = cat_name
        self._col_name_key = key_name
        self._col_cat = pl.col(self._col_name_cat)
        self._col_key = pl.col(self._col_name_key)
        return

    def break_down_data_name_column(
            self,
            col_name_data_name: str = "name",
            col_name_cat_new: str = "name_category",
            col_name_key_new: str = "name_keyword"
    ):
        df_new = self.df.with_columns(
            pl.col(col_name_data_name).str.extract(r'_([^.]+)\.(\S+)', 1).alias(col_name_cat_new),
            pl.col(col_name_data_name).str.extract(r'_([^.]+)\.(\S+)', 2).alias(col_name_key_new)
        )
        if not (
            "_" + df_new[col_name_key_new] + "." + df_new[col_name_key_new] == df_new[col_name_data_name]
        ).all():
            raise ValueError()
        self._df = df_new.select(pl.exclude(col_name_data_name))
        return

    def match_cat_key(self, col_name_category: str = "category", col_name_keyword: str = "keyword"):

        unique_nonmatching_frame_codes = self.df.filter(
            (pl.col(col_name_category) != pl.col("frame_code_category")) |
            (pl.col(col_name_keyword) != pl.col("frame_code_keyword"))
        ).select(
            pl.col(["frame_code_category", "frame_code_keyword"])
        ).unique()

        unique_matching_frame_codes = self.df.join(
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


class CIFDictCatDataFrame(PolarsDataFrame):

    def __init__(
            self,
            df: pl.DataFrame,
            col_name_cat: Optional[str] = "frame_code_category",
    ):
        super().__init__(df)
        self._col_name_cat = col_name_cat
        self._col_cat = pl.col(self._col_name_cat)
        return


class CIFDataDataFrame(PolarsDataFrame):

    def __init__(
            self,
            df: pl.DataFrame,
            col_name_block: Optional[str] = "block_code",
    ):
        super().__init__(df)
        self._col_name_block = col_name_block
        self._col_block = pl.col(self._col_name_block)
        return
