
from typing import Dict

import polars as pl


class Item:

    def __init__(self, dfs, col_name_cat: str = "category", col_name_key: str = "keyword"):
        self._dfs: Dict[str, pl.DataFrame] = dfs
        self._col_cat = pl.col(col_name_cat)
        self._col_key = pl.col(col_name_key)
        self._col_names_id = (col_name_cat, col_name_key)

    def description(self, category: str, keyword: str) -> str:
        df = self._select(df_name="item_description", category=category, keyword=keyword)
        return df[0, "description"]

    def type(self, category: str, keyword: str):
        return self._select(df_name="item_type", category=category, keyword=keyword)[0, "code"]

    def enumeration(self, category: str, keyword: str) -> pl.DataFrame:
        return self._select(
            df_name="item_enumeration", category=category, keyword=keyword
        ).select(["value", "detail"])

    def range(self, category: str, keyword: str):
        return self._select(
            df_name="item_range", category=category, keyword=keyword
        ).select(["minimum", "maximum"])

    def mandatory(self):
        pass

    def _select(self, df_name, category, keyword):
        return self._dfs[df_name].filter(
            (self._col_cat == category) & (self._col_key == keyword)
        ).select(pl.exclude(self._col_names_id))

class CIFDict:
    pass
