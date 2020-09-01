"""
opencadd.io.schema

Defines schema used accross the io module.
"""

DATAFRAME_COLUMNS = [
    ("record.name", "object"),
    ("atom.id", "Int64"),
    ("atom.name", "object"),
    ("atom.type", "object"),
    ("atom.symbol", "object"),
    ("atom.x", "float64"),
    ("atom.y", "float64"),
    ("atom.z", "float64"),
    ("atom.charge", "float64"),
    ("atom.backbone", "object"),
    ("atom.alterante_model", "object"),
    ("atom.occupancy", "float64"),
    ("atom.bfactor", "float64"),
    ("residue.pdb_id", "object"),
    ("residue.name", "object"),
    ("residue.subst_id", "Int64"),
    ("residue.subst_name", "object"),
    ("structure.chain", "object"),
]

PDB_COLUMNS = {
    0: ("record.name", str),
    1: ("atom.id", int),
    3: ("atom.name", str),
    4: ("atom.alternate_model", str),
    5: ("residue.name", str),
    7: ("structure.chain", str),
    8: ("residue.id", str),
    9: ("residue.insertion", str),
    11: ("atom.x", float),
    12: ("atom.y", float),
    13: ("atom.z", float),
    14: ("atom.occupancy", float),
    15: ("atom.bfactor", float),
    17: ("segment.id", str),
    18: ("atom.symbol", str),
    19: ("atom.charge", float),
}

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom.id", int),
        1: ("atom.name", str),
        2: ("atom.x", float),
        3: ("atom.y", float),
        4: ("atom.z", float),
        5: ("atom.type", str),
        6: ("residue.subst_id", int),
        7: ("residue.subst_name", str),
        8: ("atom.charge", float),
        9: ("atom.backbone", str),
    },
    "n_cols_9": {
        0: ("atom.id", int),
        1: ("atom.name", str),
        2: ("atom.x", float),
        3: ("atom.y", float),
        4: ("atom.z", float),
        5: ("atom.type", str),
        6: ("residue.subst_id", int),
        7: ("residue.subst_name", str),
        8: ("atom.charge", float),
    },
}
