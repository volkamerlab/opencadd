from typing import Dict, Union, Optional, Sequence
import numpy as np
import polars as pl
import pandas as pd


class CIFFile:
    def __init__(self, dict_hor: dict, df: pl.DataFrame):
        self._dict = dict_hor
        self._df = df
        return

    @property
    def block_codes(self) -> pl.Series:
        return self._df["block_code"].unique()

    @property
    def count_data_blocks(self) -> int:
        return self.block_codes.shape[0]

    @property
    def has_null_block_code(self):
        return self.block_codes.is_null().any()

    @property
    def has_null_data_name(self):
        return self._df["data_name_category"].is_null().any()

    def has_save_frame(self, per_data_block: bool = True) -> Union[pl.DataFrame, bool]:
        df = self._df.groupby("block_code").agg(~pl.col("frame_code_category").is_null().all())
        if per_data_block:
            return df
        return df["frame_code_category"].any()

    @property
    def all_data_names_are_ddl1_conformant(self) -> bool:
        return self._df.select(
            (pl.all(pl.col(["frame_code_keyword", "data_name_keyword"]).is_null()).all())
        )[0, 0]

    def data_names_are_ddl2_conformant(self):
        return self._df.select(~pl.col("data_name_keyword").is_null().any())[0, 0]

    def columns_have_same_length(self, per_loop: bool = False):
        df_per_loop = self._df.with_columns(
            pl.col("data_value").arr.lengths().alias("list_lengths")
        ).groupby(
            "loop"
        ).agg(
            (pl.col("list_lengths").n_unique() == 1).alias("same_length")
        )
        if per_loop:
            return df_per_loop
        return df_per_loop["same_length"].all()

    def loops_have_same_category(self, per_loop: bool = False):
        df_per_loop = self._df.filter(
            pl.col("loop_id") > 0
        ).groupby(
            "loop_id"
        ).agg(
            (pl.col("data_name_category").n_unique() == 1).alias("same_category")
        )
        if per_loop:
            return df_per_loop
        return df_per_loop["same_category"].all()

    def categories_have_same_length_columns(self):
        df_per_address = self._df.with_columns(
            pl.col("data_value").arr.lengths().alias("list_lengths")
        ).groupby(
            ["block_code", "frame_code_category", "frame_code_keyword", "data_name_category"]
        ).agg(
            (pl.col("list_lengths").n_unique() == 1).alias("same_length")
        )
        return df_per_address["same_length"].all()

    @property
    def has_duplicated_address(self):
        return self._df.select(pl.exclude(["data_value", "loop"])).is_duplicated().any()

    @staticmethod
    def dataframe_per_table(
            df: pl.DataFrame,
            col_name__table_id: str = "data_name_category",
            col_name__col_id: str = "data_name_keyword",
            col_name__values: str = "data_value",
            col_name__other_ids: Sequence[str] = ("block_code", "frame_code_category", "frame_code_keyword")
    ) -> Dict[str, pl.DataFrame]:

        col_name__other_ids = list(col_name__other_ids)
        table_dfs = {}

        for (table_name, *rest), table in df.groupby([col_name__table_id] + col_name__other_ids):
            new_table = table.select(
                pl.col(col_name__other_ids + [col_name__col_id, col_name__values])
            ).pivot(
                index=col_name__other_ids,
                columns=col_name__col_id,
                values=col_name__values,
                aggregate_function=None,
            ).explode(
                columns=table[col_name__col_id].unique().to_list()
            )

            table_dfs[table_name] = pl.concat(
                [table_dfs.setdefault(table_name, pl.DataFrame()), new_table],
                how="diagonal"
            )
        return table_dfs

    def data_block(self, block_code_or_idx: Optional[Union[str, int]] = None):
        if isinstance(block_code_or_idx, int):
            block_code = self.block_codes[block_code_or_idx]
        else:
            block_code = block_code_or_idx
        return CIFDataBlock(
            block_code=block_code,
            dic=self._dict[block_code],
            df=self._df[self._df.block_code == block_code].drop("block_code", axis=1)
        )

    def dataframe_per_data_category(self):
        category_dfs = dict()

        for block_code, data_block_dic in self._dict.items():
            for frame_code_category, save_frame_category_dic in data_block_dic.items():
                for frame_code_keyword, save_frame_keyword_dic in save_frame_category_dic.items():
                    for data_name_category, data_category_dic in save_frame_keyword_dic.items():

                        category_df = category_dfs.setdefault(
                            data_name_category,
                            dict(block_code=[], frame_code_category=[], frame_code_keyword=[])
                        )
                        lens = np.array(
                            [
                                len(data_value) if not isinstance(data_value, str) else 0
                                for data_value in data_category_dic.values()
                            ]
                        )
                        unique_lens = np.unique(lens)
                        if unique_lens.size != 1:
                            raise ValueError
                        if unique_lens[0] > 0:
                            cat_type = "tab"
                            num_vals = unique_lens[0]
                        else:
                            cat_type = "single"
                            num_vals = 1

                        for data_name_keyword, data_value in data_category_dic.items():
                            col = category_df.setdefault(data_name_keyword, [])
                            if len(col) != len(category_df["block_code"]):
                                col.extend([None] * (len(category_df["block_code"]) - len(col)))
                            if cat_type == "tab":
                                col.extend(data_value)
                            else:
                                col.append(data_value)

                        category_df["block_code"].extend([block_code] * num_vals)
                        category_df["frame_code_category"].extend([frame_code_category] * num_vals)
                        category_df["frame_code_keyword"].extend([frame_code_keyword] * num_vals)

        category_dfss = dict()
        for cat_name, category_df in category_dfs.items():
            num_entries = len(category_df["block_code"])
            for col_name, col in category_df.items():
                if len(col) != num_entries:
                    col.extend([None] * (num_entries - len(col)))
            category_dfss[cat_name] = pl.DataFrame(category_df)
        return category_dfss


    @staticmethod
    def _dic_to_df(dic):
        block_codes = []
        frame_code_categories = []
        frame_code_keywords = []
        data_name_categories = []
        data_name_keywords = []
        data_values = []
        value_dimensions = []

        for block_code, data_block_dic in dic.items():
            for frame_code_category, save_frame_category_dic in data_block_dic.items():
                for frame_code_keyword, save_frame_keyword_dic in save_frame_category_dic.items():
                    for data_name_category, data_category_dic in save_frame_keyword_dic.items():
                        for data_name_keyword, data_value in data_category_dic.items():
                            block_codes.append(block_code)
                            frame_code_categories.append(frame_code_category)
                            frame_code_keywords.append(frame_code_keyword)
                            data_name_categories.append(data_name_category)
                            data_name_keywords.append(data_name_keyword)
                            data_values.append(data_value)
                            value_dimensions.append(value_dimension)

        data = dict(
            block_code=block_codes,
            frame_code_category=frame_code_categories,
            frame_code_keyword=frame_code_keywords,
            data_name_category=data_name_categories,
            data_name_keyword=data_name_keywords,
            data_value=data_values,
            value_dimension=value_dimensions,
        )
        return pd.DataFrame(data)


class CIFDataBlock:
    def __init__(self, block_code: str, dic, df):
        self.block_code: str = block_code

        frame_code_categories = list(dic.keys())

        if frame_code_categories == [None]:
            frame_code_keywords = list(dic[None].keys())
            if frame_code_keywords == [None]:
                dic = dic[None][None]
                df = df.drop(["frame_code_category", "frame_code_keyword"], axis=1)


        self.dic = dic
        self.df = df
        return

def f2(
    df: pl.DataFrame,
    col_name__table_id: str = "table_id",
    col_name__col_id: str = "column_id",
    col_name__values: str = "data",
    col_name__other_ids: Sequence[str] = ("db_id",)
) -> Dict[str, pl.DataFrame]:
    x = {
            tab: df.filter(
                pl.col(col_name__table_id) == tab
            ).with_columns(
                datai=pl.arange(0, pl.col(col_name__values).arr.lengths())
            ).explode(
                [col_name__values, 'datai']
            ).pivot(
                values=col_name__values,
                index=[*col_name__other_ids, 'datai'],
                columns=col_name__col_id,
                aggregate_function='first'
            ).drop(
                'datai'
            ) for tab in df.get_column(col_name__table_id).unique()
    }
    return x


def f3(
        df: pl.DataFrame,
        col_name__table_id: str = "table_id",
        col_name__col_id: str = "column_id",
        col_name__values: str = "data",
        col_name__other_ids: Sequence[str] = ("db_id",)
) -> Dict[str, pl.DataFrame]:
    df_long = df.with_row_count().explode(
        col_name__values
    ).with_columns(
        len=pl.count().over(col_name__table_id, col_name__col_id)
    ).with_columns(
        max_len=pl.max("len").over(col_name__table_id)
    )

    table_dfs = df_long.with_columns(
        pl.col("row_nr") + 1
    ).join_asof(
        df_long, by=col_name__table_id, on="row_nr"
    ).with_columns(
        diff=pl.col("max_len") - pl.col("len_right")
    ).with_columns(
        pl.when(pl.col("len") != pl.col("max_len")).then(pl.col("diff")).fill_null(0)
    ).with_columns(
        pl.col("row_nr").cumcount().over(col_name__table_id, col_name__col_id) + pl.col("diff")
    ).pivot(
        index=["row_nr", *col_name__other_ids, col_name__table_id],
        columns=col_name__col_id,
        values=col_name__values,
        aggregate_function=None,
    ).groupby(
        "row_nr", col_name__table_id, maintain_order=True
    ).agg(
        pl.all().drop_nulls().first()
    ).partition_by(
        col_name__table_id
    )

    return table_dfs
