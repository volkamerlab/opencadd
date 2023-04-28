from typing import Dict, Union, Optional, Sequence, Literal
import numpy as np
import polars as pl

import opencadd._exceptions as _exception

__author__ = "Armin Ariamajd"


class DDL2CIFFile:

    def __init__(self, df: pl.DataFrame):
        _exception.raise_for_type(
            f"{self.__class__.__name__}.__init__",
            ("df", df, pl.DataFrame)
        )
        column_names = df.columns
        if len(column_names) != 6:
            raise ValueError()
        for column_name in (
                'block_code',
                'frame_code_category',
                'frame_code_keyword',
                'data_name_category',
                'data_name_keyword',
                'data_value',
        ):
            if column_name not in column_names:
                raise ValueError()
        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    @property
    def block_codes(self) -> pl.Series:
        return self._df["block_code"].unique()

    @property
    def count_data_blocks(self) -> int:
        return self.block_codes.shape[0]

    def has_save_frame(self, per_data_block: bool = False) -> Union[pl.DataFrame, bool]:
        df = self._df.groupby(
            "block_code"
        ).agg(
            ~pl.col("frame_code_category").is_null().all().alias("has_save_frame")
        )
        if per_data_block:
            return df
        return df["has_save_frame"].any()

    def data_block(self, block_code_or_idx: Optional[Union[str, int]] = 0):
        if isinstance(block_code_or_idx, int):
            block_code = self.block_codes[block_code_or_idx]
        else:
            block_code = block_code_or_idx
        return DDL2CIFBlock(
            block_code=block_code,
            df=self.df.filter(pl.col("block_code") == block_code).select(pl.exclude("block_code"))
        )

    def extract(
            self,
            part: Literal["data", "def", "def_cat", "def_key", "all"] = "data",
            reduce: bool = True
    ):
        _exception.raise_literal(
            "_extract",
            ("part", part, ["data", "def", "def_cat", "def_key", "all"])
        )

        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).extract(part=part)
        if part == "all":
            parts = ("data", "def_cat", "def_key")
            sub_dfs = [_extract(df=self._df, part=part) for part in parts]
            return {
                part: datastruct(df=sub_df) if sub_df.shape[0] != 0 else None
                for part, datastruct, sub_df in zip(
                    parts, (DDL2CIFDataFile, DDL2CIFCatDefFile, DDL2CIFDefFile), sub_dfs
                )
            }
        sub_df = _extract(df=self._df, part=part)
        if sub_df.shape[0] == 0:
            return
        if part == "data":
            return DDL2CIFDataFile(df=sub_df)
        if part in ("def", "def_key"):
            return DDL2CIFDefFile(df=sub_df)
        return DDL2CIFCatDefFile(df=sub_df)

    def df_per_category(
            self,
            part: Optional[Literal["data", "def", "def_cat", "def_key", "all"]] = None,
            reduce: bool = True
    ):
        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).df_per_category(part=part)
        if part is None:
            return _dataframe_per_table(self._df)
        if part == "all":
            return {
                part: datastruct.df_per_category() if datastruct is not None else None
                for part, datastruct in self.extract(part=part, reduce=False).items()
            }
        return self.extract(part=part, reduce=False).df_per_category()


class DDL2CIFDataFile:

    def __init__(self, df: pl.DataFrame):

        column_names = df.columns
        if len(column_names) != 4:
            raise ValueError()
        for column_name in (
                'block_code',
                'data_name_category',
                'data_name_keyword',
                'data_value',
        ):
            if column_name not in column_names:
                raise ValueError()

        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    @property
    def block_codes(self) -> pl.Series:
        return self._df["block_code"].unique()

    @property
    def count_data_blocks(self) -> int:
        return self.block_codes.shape[0]

    def df_per_category(self, reduce: bool = True):
        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).df_per_category()
        return _dataframe_per_table(self._df, col_name__other_ids=("block_code", ))

    def data_block(self, block_code_or_idx: Optional[Union[str, int]] = 0):
        if isinstance(block_code_or_idx, int):
            block_code = self.block_codes[block_code_or_idx]
        else:
            block_code = block_code_or_idx
        return DDL2CIFDataBlock(
            block_code=block_code,
            df=self.df.filter(pl.col("block_code") == block_code).select(pl.exclude("block_code"))
        )


class DDL2CIFDefFile:

    def __init__(self, df: pl.DataFrame):
        column_names = df.columns
        if len(column_names) != 6:
            raise ValueError()
        for column_name in (
                'block_code',
                'frame_code_category',
                'frame_code_keyword',
                'data_name_category',
                'data_name_keyword',
                'data_value',
        ):
            if column_name not in column_names:
                raise ValueError()
        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    @property
    def block_codes(self) -> pl.Series:
        return self._df["block_code"].unique()

    @property
    def count_data_blocks(self) -> int:
        return self.block_codes.shape[0]

    def data_block(self, block_code_or_idx: Optional[Union[str, int]] = 0):
        if isinstance(block_code_or_idx, int):
            block_code = self.block_codes[block_code_or_idx]
        else:
            block_code = block_code_or_idx
        return DDL2CIFDefBlock(
            block_code=block_code,
            df=self.df.filter(pl.col("block_code") == block_code).select(pl.exclude("block_code"))
        )

    def extract(self, part: Literal["def_cat", "def_key", "all"] = "all", reduce: bool = True):
        _exception.raise_literal(
            "_extract",
            ("part", part, ["def_cat", "def_key", "all"])
        )
        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).extract(part=part)
        if part == "all":
            parts = ("def_cat", "def_key")
            sub_dfs = [_extract(df=self._df, part=part) for part in parts]
            return {
                part: datastruct(df=sub_df) if sub_df.shape[0] != 0 else None
                for part, datastruct, sub_df in zip(
                    parts, (DDL2CIFCatDefFile, DDL2CIFDefFile), sub_dfs
                )
            }
        sub_df = _extract(df=self._df, part=part)
        if sub_df.shape[0] == 0:
            return
        if part == "def_cat":
            return DDL2CIFCatDefFile(df=sub_df)
        return DDL2CIFDefFile(df=sub_df)

    def df_per_category(
            self,
            reduce: bool = True,
            part: Optional[Literal["def_cat", "def_key", "all"]] = None,
    ):
        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).df_per_category()
        if part is None:
            return _dataframe_per_table(self._df)
        if part == "all":
            return {
                part: datastruct.df_per_category() if datastruct is not None else None
                for part, datastruct in self.extract(part=part, reduce=False)
            }
        return self.extract(part=part, reduce=False).df_per_category()


class DDL2CIFCatDefFile:

    def __init__(self, df: pl.DataFrame):
        column_names = df.columns
        if len(column_names) != 5:
            raise ValueError()
        for column_name in (
                'block_code',
                'frame_code_category',
                'data_name_category',
                'data_name_keyword',
                'data_value',
        ):
            if column_name not in column_names:
                raise ValueError()
        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    @property
    def block_codes(self) -> pl.Series:
        return self._df["block_code"].unique()

    @property
    def count_data_blocks(self) -> int:
        return self.block_codes.shape[0]

    def data_block(self, block_code_or_idx: Optional[Union[str, int]] = 0):
        if isinstance(block_code_or_idx, int):
            block_code = self.block_codes[block_code_or_idx]
        else:
            block_code = block_code_or_idx
        return DDL2CIFCatDefBlock(
            block_code=block_code,
            df=self.df.filter(pl.col("block_code") == block_code).select(pl.exclude("block_code"))
        )

    def df_per_category(self, reduce: bool = True):
        if reduce and self.count_data_blocks == 1:
            return self.data_block(block_code_or_idx=0).df_per_category()
        return _dataframe_per_table(self._df, col_name__other_ids=("block_code", "frame_code_category"))


class DDL2CIFBlock:

    def __init__(self, block_code: str, df):
        self.block_code: str = block_code
        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    def df_per_category(self, part: Optional[Literal["data", "def", "def_cat", "def_key", "all"]] = "all"):
        if part is None:
            return _dataframe_per_table(
                self.df, col_name__other_ids=("frame_code_category", "frame_code_keyword")
            )
        if part == "all":
            return {
                part: datastruct.df_per_category() if datastruct is not None else None
                for part, datastruct in self.extract(part="all").items()
            }
        return self.extract(part=part).df_per_category()

    def extract(
        self,
        part: Literal["data", "def", "def_cat", "def_key", "all"] = "all",
    ):
        _exception.raise_literal(
            "_extract",
            ("part", part, ["data", "def", "def_cat", "def_key", "all"])
        )
        if part == "all":
            parts = ("data", "def_cat", "def_key")
            sub_dfs = [_extract(df=self._df, part=part) for part in parts]
            return {
                part: datastruct(block_code=self.block_code, df=sub_df) if sub_df.shape[0] != 0 else None
                for part, datastruct, sub_df in zip(
                    parts, (DDL2CIFDataBlock, DDL2CIFCatDefBlock, DDL2CIFDefBlock), sub_dfs
                )
            }
        sub_df = _extract(df=self._df, part=part)
        if sub_df.shape[0] == 0:
            return
        if part == "data":
            return DDL2CIFDataBlock(block_code=self.block_code, df=sub_df)
        if part in ("def", "def_key"):
            return DDL2CIFDefBlock(block_code=self.block_code, df=sub_df)
        return DDL2CIFCatDefBlock(block_code=self.block_code, df=sub_df)


class DDL2CIFDataBlock:

    def __init__(self, block_code: str, df: pl.DataFrame):
        self._block_code = block_code
        self._df = df
        return

    @property
    def block_code(self):
        return self._block_code

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    def df_per_category(self):
        return _dataframe_per_table(self.df, col_name__other_ids=tuple())


class DDL2CIFDefBlock:

    def __init__(self, block_code: str, df: pl.DataFrame):
        self._block_code = block_code
        self._df = df
        return

    @property
    def block_code(self):
        return self._block_code

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    def extract(
        self,
        part: Literal["def_cat", "def_key", "all"] = "all"
    ):
        _exception.raise_literal(
            "_extract",
            ("part", part, ["def_cat", "def_key", "all"])
        )
        if part == "all":
            parts = ("def_cat", "def_key")
            sub_dfs = [_extract(df=self._df, part=part) for part in parts]
            return {
                part: datastruct(df=sub_df) if sub_df.shape[0] != 0 else None
                for part, datastruct, sub_df in zip(
                    parts, (DDL2CIFCatDefFile, DDL2CIFDefFile), sub_dfs
                )
            }
        sub_df = _extract(df=self._df, part=part)
        if sub_df.shape[0] == 0:
            return
        if part == "def_cat":
            return DDL2CIFCatDefBlock(block_code=self.block_code, df=sub_df)
        return DDL2CIFDefBlock(block_code=self.block_code, df=sub_df)

    def df_per_category(self):
        return _dataframe_per_table(self.df, col_name__other_ids=("frame_code_category", "frame_code_keyword"))


class DDL2CIFCatDefBlock:

    def __init__(self, block_code: str, df: pl.DataFrame):
        self._block_code = block_code
        self._df = df
        return

    @property
    def block_code(self):
        return self._block_code

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    def df_per_category(self):
        return _dataframe_per_table(self.df, col_name__other_ids=("frame_code_category", ))



def _dataframe_per_table(
        df: pl.DataFrame,
        col_name__table_id: str = "data_name_category",
        col_name__col_id: str = "data_name_keyword",
        col_name__values: str = "data_value",
        col_name__other_ids: Sequence[str] = ("block_code", "frame_code_category", "frame_code_keyword")
) -> Dict[str, pl.DataFrame]:
    table_dfs = {
        table_id: df.filter(
            pl.col(col_name__table_id) == table_id
        ).with_columns(
            idx_data=pl.arange(0, pl.col(col_name__values).arr.lengths())
        ).explode(
            [col_name__values, 'idx_data']
        ).pivot(
            values=col_name__values,
            index=[*col_name__other_ids, 'idx_data'],
            columns=col_name__col_id,
            aggregate_function='first'
        ).drop(
            'idx_data'
        ) for table_id in df.get_column(col_name__table_id).unique()
    }
    return table_dfs


def _extract(df, part: Literal["data", "def", "def_cat", "def_key"] = "data"):
    col_frame_cat = pl.col("frame_code_category")
    col_frame_key = pl.col("frame_code_keyword")
    condition = col_frame_cat.is_null() if part == "data" else col_frame_cat.is_not_null()
    if part == "def_cat":
        condition &= col_frame_key.is_null()
        final_columns = pl.exclude(["frame_code_keyword"])
    elif part == "def_key":
        condition &= col_frame_key.is_not_null()
        final_columns = pl.all()
    elif part == "data":
        final_columns = pl.exclude(["frame_code_category", "frame_code_keyword"])
    else:
        final_columns = pl.all()
    return df.filter(condition).select(final_columns)
