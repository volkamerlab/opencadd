"""
opencadd.io.schema

Defines schema used across the io module.
"""

DATAFRAME_COLUMNS = {
    "default": [
        ("atom.id", "int32"),
        ("atom.name", "object"),
        ("atom.x", "float32"),
        ("atom.y", "float32"),
        ("atom.z", "float32"),
        ("residue.id", "object"),
        ("residue.name", "object"),
    ],
    "verbose": [
        ("atom.type", "object"),
        ("residue.subst_id", "Int64"),
        ("residue.subst_name", "object"),
        ("record.name", "object"),
        ("atom.symbol", "object"),
        ("atom.charge", "float32"),
        ("atom.status_bit", "object"),
        ("atom.occupancy", "float32"),
        ("atom.bfactor", "float32"),
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
    8: ("residue.id", "object"),
    9: ("residue.insertion", "object"),
    11: ("atom.x", "float32"),
    12: ("atom.y", "float32"),
    13: ("atom.z", "float32"),
    14: ("atom.occupancy", "float32"),
    15: ("atom.bfactor", "float32"),
    17: ("segment.id", "object"),
    18: ("atom.symbol", "object"),
    19: ("atom.charge", "float32"),
}

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom.id", "int32"),
        1: ("atom.name", "object"),
        2: ("atom.x", "float32"),
        3: ("atom.y", "float32"),
        4: ("atom.z", "float32"),
        5: ("atom.type", "object"),
        6: ("residue.subst_id", "int32"),
        7: ("residue.subst_name", "object"),
        8: ("atom.charge", "float32"),
        9: ("atom.status_bit", "object"),
    },
    "n_cols_9": {
        0: ("atom.id", "int32"),
        1: ("atom.name", "object"),
        2: ("atom.x", "float32"),
        3: ("atom.y", "float32"),
        4: ("atom.z", "float32"),
        5: ("atom.type", "object"),
        6: ("residue.subst_id", "int32"),
        7: ("residue.subst_name", "object"),
        8: ("atom.charge", "float32"),
    },
}
