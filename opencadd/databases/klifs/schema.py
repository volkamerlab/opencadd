"""
opencadd.databases.klifs.schema

Defines schema used across the klifs module.
"""

from opencadd.io.schema import DATAFRAME_COLUMNS
from pathlib import Path

import pandas as pd

# PATH_DATA = Path(__name__).parent / "opencadd/data"
PATH_DATA = Path(__file__).parent / "../../data"

FIELDS = pd.read_csv(PATH_DATA / "klifs_fields.csv")

# "kinase_names": Kinase names: "kinase.gene_name (kinase.klifs_name)"
# "kinase.klifs_name": Depending on availability: Manning name or UniProt gene name
# "kinase.full_name": Depending on availability: HGNC gene name or Manning name or UniProt gene name
# "kinase.gene_name": HGNC or MGI name
# "kinase.uniprot": UniProt accession
LOCAL_COLUMNS_MAPPING = {
    "klifs_export": FIELDS[["local.klifs_export", "df.name"]]
    .dropna(how="any")
    .set_index("local.klifs_export")
    .squeeze()
    .to_dict(),
    "klifs_overview": FIELDS[["local.klifs_overview", "df.name"]]
    .dropna(how="any")
    .set_index("local.klifs_overview")
    .squeeze()
    .to_dict(),
}
REMOTE_COLUMNS_MAPPING = {
    field_type: fields[["remote.fields", "df.name"]]
    .dropna(how="any")
    .set_index("remote.fields")
    .squeeze()
    .to_dict()
    for field_type, fields in FIELDS.groupby("field_type")
}
DATAFRAME_COLUMNS = {}
for field_type, fields in FIELDS.groupby("field_type"):
    fields = fields[["df.name", "df.type"]].dropna(how="any").values.tolist()
    fields = [tuple(i) for i in fields]
    DATAFRAME_COLUMNS[field_type] = fields

POCKET_KLIFS_REGIONS = [
    (1, "I"),
    (2, "I"),
    (3, "I"),
    (4, "g.l"),
    (5, "g.l"),
    (6, "g.l"),
    (7, "g.l"),
    (8, "g.l"),
    (9, "g.l"),
    (10, "II"),
    (11, "II"),
    (12, "II"),
    (13, "II"),
    (14, "III"),
    (15, "III"),
    (16, "III"),
    (17, "III"),
    (18, "III"),
    (19, "III"),
    (20, "αC"),
    (21, "αC"),
    (22, "αC"),
    (23, "αC"),
    (24, "αC"),
    (25, "αC"),
    (26, "αC"),
    (27, "αC"),
    (28, "αC"),
    (29, "αC"),
    (30, "αC"),
    (31, "b.l"),
    (32, "b.l"),
    (33, "b.l"),
    (34, "b.l"),
    (35, "b.l"),
    (36, "b.l"),
    (37, "b.l"),
    (38, "IV"),
    (39, "IV"),
    (40, "IV"),
    (41, "IV"),
    (42, "V"),
    (43, "V"),
    (44, "V"),
    (45, "GK"),
    (46, "hinge"),
    (47, "hinge"),
    (48, "hinge"),
    (49, "linker"),
    (50, "linker"),
    (51, "linker"),
    (52, "linker"),
    (53, "αD"),
    (54, "αD"),
    (55, "αD"),
    (56, "αD"),
    (57, "αD"),
    (58, "αD"),
    (59, "αD"),
    (60, "αE"),
    (61, "αE"),
    (62, "αE"),
    (63, "αE"),
    (64, "αE"),
    (65, "VI"),
    (66, "VI"),
    (67, "VI"),
    (68, "c.l"),
    (69, "c.l"),
    (70, "c.l"),
    (71, "c.l"),
    (72, "c.l"),
    (73, "c.l"),
    (74, "c.l"),
    (75, "c.l"),
    (76, "VII"),
    (77, "VII"),
    (78, "VII"),
    (79, "VIII"),
    (80, "xDFG"),
    (81, "xDFG"),
    (82, "xDFG"),
    (83, "xDFG"),
    (84, "a.l"),
    (85, "a.l"),
]
POCKET_KLIFS_REGIONS = pd.DataFrame(
    [
        (klifs_id, klifs_region, ".".join([klifs_region, str(klifs_id)]))
        for (klifs_id, klifs_region) in POCKET_KLIFS_REGIONS
    ],
    columns=["residue.klifs_id", "residue.klifs_region", "residue.klifs_region_id"],
)

POCKET_KLIFS_REGION_COLORS = {
    "I": "khaki",
    "g.l": "green",
    "II": "khaki",
    "III": "khaki",
    "αC": "red",
    "b.l": "green",
    "IV": "khaki",
    "V": "khaki",
    "GK": "orange",
    "hinge": "magenta",
    "linker": "cyan",
    "αD": "red",
    "αE": "red",
    "VI": "khaki",
    "c.l": "darkorange",
    "VII": "khaki",
    "VIII": "khaki",
    "xDFG": "cornflowerblue",
    "a.l": "cornflowerblue",
}
