"""
opencadd.io.schema

Defines schema used across the io module.
"""

DATAFRAME_COLUMNS = {
    "default": [
        ("atom.id", "int64"),
        ("atom.name", "object"),
        ("atom.x", "float64"),
        ("atom.y", "float64"),
        ("atom.z", "float64"),
        ("residue.pdb_id", "object"),
        ("residue.name", "object"),
    ],
    "verbose": [
        ("atom.type", "object"),
        ("residue.subst_id", "Int64"),
        ("residue.subst_name", "object"),
        ("record.name", "object"),
        ("atom.symbol", "object"),
        ("atom.charge", "float64"),
        ("atom.status_bit", "object"),
        ("atom.occupancy", "float64"),
        ("atom.bfactor", "float64"),
        ("atom.alternative_model", "object"),
        ("structure.chain", "object"),
    ],
}


PDB_COLUMNS = {
    0: ("record.name", "object"),
    1: ("atom.id", "Int64"),
    3: ("atom.name", "object"),
    4: ("atom.alternative_model", "object"),
    5: ("residue.name", "object"),
    7: ("structure.chain", "object"),
    8: ("residue.pdb_id", "object"),
    9: ("residue.insertion", "object"),
    11: ("atom.x", "float64"),
    12: ("atom.y", "float64"),
    13: ("atom.z", "float64"),
    14: ("atom.occupancy", "float64"),
    15: ("atom.bfactor", "float64"),
    17: ("segment.id", "object"),
    18: ("atom.symbol", "object"),
    19: ("atom.charge", "float64"),
}

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom.id", "int64"),
        1: ("atom.name", "object"),
        2: ("atom.x", "float64"),
        3: ("atom.y", "float64"),
        4: ("atom.z", "float64"),
        5: ("atom.type", "object"),
        6: ("residue.subst_id", "int64"),
        7: ("residue.subst_name", "object"),
        8: ("atom.charge", "float64"),
        9: ("atom.status_bit", "object"),
    },
    "n_cols_9": {
        0: ("atom.id", "int64"),
        1: ("atom.name", "object"),
        2: ("atom.x", "float64"),
        3: ("atom.y", "float64"),
        4: ("atom.z", "float64"),
        5: ("atom.type", "object"),
        6: ("residue.subst_id", "int64"),
        7: ("residue.subst_name", "object"),
        8: ("atom.charge", "float64"),
    },
}
