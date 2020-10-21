"""
opencadd.io.schema

Defines schema used across the io module.
"""

DATAFRAME_COLUMNS = {
    "default": [
        ("atom.id", "int32"),
        ("atom.name", "string"),
        ("atom.x", "float32"),
        ("atom.y", "float32"),
        ("atom.z", "float32"),
        ("residue.id", "string"),
        ("residue.name", "string"),
    ],
    "verbose": [
        ("atom.type", "string"),
        ("residue.subst_id", "Int64"),
        ("residue.subst_name", "string"),
        ("record.name", "string"),
        ("atom.symbol", "string"),
        ("atom.charge", "float32"),
        ("atom.status_bit", "string"),
        ("atom.occupancy", "float32"),
        ("atom.bfactor", "float32"),
        ("atom.alternative_model", "string"),
        ("structure.chain", "string"),
    ],
}


PDB_COLUMNS = {
    0: ("record.name", "string"),
    1: ("atom.id", "Int64"),
    3: ("atom.name", "string"),
    4: ("atom.alternative_model", "string"),
    5: ("residue.name", "string"),
    7: ("structure.chain", "string"),
    8: ("residue.id", "string"),
    9: ("residue.insertion", "string"),
    11: ("atom.x", "float32"),
    12: ("atom.y", "float32"),
    13: ("atom.z", "float32"),
    14: ("atom.occupancy", "float32"),
    15: ("atom.bfactor", "float32"),
    17: ("segment.id", "string"),
    18: ("atom.symbol", "string"),
    19: ("atom.charge", "float32"),
}

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom.id", "int32"),
        1: ("atom.name", "string"),
        2: ("atom.x", "float32"),
        3: ("atom.y", "float32"),
        4: ("atom.z", "float32"),
        5: ("atom.type", "string"),
        6: ("residue.subst_id", "int32"),
        7: ("residue.subst_name", "string"),
        8: ("atom.charge", "float32"),
        9: ("atom.status_bit", "string"),
    },
    "n_cols_9": {
        0: ("atom.id", "int32"),
        1: ("atom.name", "string"),
        2: ("atom.x", "float32"),
        3: ("atom.y", "float32"),
        4: ("atom.z", "float32"),
        5: ("atom.type", "string"),
        6: ("residue.subst_id", "int32"),
        7: ("residue.subst_name", "string"),
        8: ("atom.charge", "float32"),
    },
}
