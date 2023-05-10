from typing import Union, Sequence
import functools
import operator

import polars as pl
import numpy as np


class PolarsDataFrame:

    def __init__(self, df: pl.DataFrame):
        self._df = df
        return

    @property
    def df(self) -> pl.DataFrame:
        return self._df

    @property
    def keywords(self) -> list:
        return self.df.columns

    def has_keyword(self, keywords: Union[str, Sequence[str]], reduce: bool = True):
        is_in = np.isin(keywords, self.keywords)
        return is_in.all() if reduce else is_in

    def filter(self, *args):
        predicate = functools.reduce(operator.and_, [pl.col(arg[0]) == arg[1] for arg in args])
        return self.df.filter(predicate)

    def column_is_subset(self, col_1: str, col_2: str):
        per_row = self.df.select(
            (
                    (pl.col(col_2).is_null()) |
                    (pl.col(col_2) == pl.col(col_1))
            ).alias("2_is_sub_1"),
            (
                    (pl.col(col_1).is_null()) |
                    (pl.col(col_2) == pl.col(col_1))
            ).alias("1_is_sub_2")
        )
        reduced = per_row.select(pl.all().all())
        if reduced.select(pl.all(pl.all()).alias("both"))[0, "both"]:
            return True, True, per_row
        if reduced[0, "1_is_sub_2"]:
            return True, col_2, per_row
        if reduced[0, "2_is_sub_1"]:
            return True, col_1, per_row
        return False, False, per_row
