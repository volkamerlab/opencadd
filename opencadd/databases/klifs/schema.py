"""
schema.py

Defines schema used accross the klifs module.
"""


LOCAL_COLUMNS_MAPPING = {
    "klifs_export": {
        "NAME": "kinase.name",
        "FAMILY": "kinase.family",
        "GROUPS": "kinase.group",
        "PDB": "structure.pdb",
        "CHAIN": "structure.chain",
        "ALTERNATE_MODEL": "structure.alternate_model",
        "SPECIES": "species.klifs",
        "LIGAND": "ligand.name",
        "PDB_IDENTIFIER": "ligand.pdb",
        "ALLOSTERIC_NAME": "ligand.name_allosteric",
        "ALLOSTERIC_PDB": "ligand.pdb_allosteric",
        "DFG": "structure.dfg",
        "AC_HELIX": "structure.ac_helix",
    },
    "klifs_overview": {
        "species": "species.klifs",
        "kinase": "kinase.name",
        "pdb": "structure.pdb",
        "alt": "structure.alternate_model",
        "chain": "structure.chain",
        "orthosteric_PDB": "ligand.pdb",
        "allosteric_PDB": "ligand.pdb_allosteric",
        "rmsd1": "structure.rmsd1",
        "rmsd2": "structure.rmsd2",
        "qualityscore": "structure.qualityscore",
        "pocket": "kinase.pocket",
        "resolution": "structure.resolution",
        "missing_residues": "structure.missing_residues",
        "missing_atoms": "structure.missing_atoms",
        "full_ifp": "interaction.fingerprint",
        "fp_i": "structure.fp_i",
        "fp_ii": "structure.fp_ii",
        "bp_i_a": "structure.bp_i_a",
        "bp_i_b": "structure.bp_i_b",
        "bp_ii_in": "structure.bp_ii_in",
        "bp_ii_a_in": "structure.bp_ii_a_in",
        "bp_ii_b_in": "structure.bp_ii_b_in",
        "bp_ii_out": "structure.bp_ii_out",
        "bp_ii_b": "structure.bp_ii_b",
        "bp_iii": "structure.bp_iii",
        "bp_iv": "structure.bp_iv",
        "bp_v": "structure.bp_v",
    },
}

REMOTE_COLUMNS_MAPPING = {
    # Information.get_kinase_information()
    "kinases": {
        "kinase_ID": "kinase.id",
        "name": "kinase.name",
        "HGNC": "kinase.hgnc",
        "family": "kinase.family",
        "group": "kinase.group",
        "kinase_class": "kinase.class",
        "species": "species.klifs",
        "full_name": "kinase.name_full",
        "uniprot": "kinase.uniprot",
        "iuphar": "kinase.iuphar",
        "pocket": "kinase.pocket",
    },
    # Ligands.get_ligands_list
    "ligands": {
        "ligand_ID": "ligand.id",
        "PDB-code": "ligand.pdb",
        "Name": "ligand.name",
        "SMILES": "ligand.smiles",
        "InChIKey": "ligand.inchikey",
    },
    # Structures.get_structure_list()
    # Structures.get_structure_lists()
    "structures": {
        "structure_ID": "structure.id",
        "kinase": "kinase.name",
        "species": "species.klifs",
        "kinase_ID": "kinase.id",
        "pdb": "structure.pdb",
        "alt": "structure.alternate_model",
        "chain": "structure.chain",
        "rmsd1": "structure.rmsd1",
        "rmsd2": "structure.rmsd2",
        "pocket": "kinase.pocket",
        "resolution": "structure.resolution",
        "quality_score": "structure.qualityscore",
        "missing_residues": "structure.missing_residues",
        "missing_atoms": "structure.missing_atoms",
        "ligand": "ligand.pdb",
        "allosteric_ligand": "ligand.pdb_allosteric",
        "DFG": "structure.dfg",
        "aC_helix": "structure.ac_helix",
        "Grich_distance": "structure.grich_distance",
        "Grich_angle": "structure.grich_angle",
        "Grich_rotation": "structure.grich_rotation",
        "front": "structure.front",
        "gate": "structure.gate",
        "back": "structure.back",
        "fp_I": "structure.fp_i",
        "fp_II": "structure.fp_ii",
        "bp_I_A": "structure.bp_i_a",
        "bp_I_B": "structure.bp_i_b",
        "bp_II_in": "structure.bp_ii_in",
        "bp_II_A_in": "structure.bp_ii_a_in",
        "bp_II_B_in": "structure.bp_ii_b_in",
        "bp_II_out": "structure.bp_ii_out",
        "bp_II_B": "structure.bp_ii_b",
        "bp_III": "structure.bp_iii",
        "bp_IV": "structure.bp_iv",
        "bp_V": "structure.bp_v",
    },
    # Ligands.get_bioactivity_list_id()
    "bioactivities": {
        "pref_name": "kinase.pref_name",
        "accession": "kinase.uniprot",
        "organism": "species.chembl",
        "standard_type": "ligand.bioactivity_standard_type",
        "standard_relation": "ligand.bioactivity_standard_relation",
        "standard_value": "ligand.bioactivity_standard_value",
        "standard_units": "ligand.bioactivity_standard_units",
        "pchembl_value": "ligand.bioactivity_pchembl_value",
    },
    # Interactions.get_interactions_get_IFP()
    "interactions": {"structure_ID": "structure.id", "IFP": "interaction.fingerprint",},
    # Interactions.get_interactions_get_types()
    "interaction_types": {"position": "interaction.id", "name": "interaction.name",},
    # Interactions.get_interactions_match_residues()
    "pockets": {
        "index": "residue.klifs_id",
        "Xray_position": "residue.pdb_id",
        "KLIFS_position": "residue.klifs_region",
    },
}

COLUMN_NAMES = {
    "kinase_groups": ["kinase.group"],
    "kinase_families": ["kinase.family"],
    "kinases_all": ["kinase.id", "kinase.name", "kinase.name_full", "species.klifs",],
    "kinases": [
        "kinase.id",
        "kinase.name",
        "kinase.hgnc",
        "kinase.family",
        "kinase.group",
        "kinase.class",
        "species.klifs",
        "kinase.name_full",
        "kinase.uniprot",
        "kinase.iuphar",
        "kinase.pocket",
    ],
    "ligands": ["ligand.id", "ligand.pdb", "ligand.name", "ligand.smiles", "ligand.inchikey"],
    "structures": [
        "structure.id",
        "structure.pdb",
        "structure.alternate_model",
        "structure.chain",
        "species.klifs",
        "kinase.id",
        "kinase.name",
        # "kinase.name_all",  # Excluded, otherwise operations like drop_duplicates() do not work
        "kinase.family",
        "kinase.group",
        "kinase.pocket",
        "ligand.pdb",
        "ligand.pdb_allosteric",
        "ligand.name",
        "ligand.name_allosteric",
        "structure.dfg",
        "structure.ac_helix",
        "structure.resolution",
        "structure.qualityscore",
        "structure.missing_residues",
        "structure.missing_atoms",
        "structure.rmsd1",
        "structure.rmsd2",
        "structure.front",
        "structure.gate",
        "structure.back",
        "structure.fp_i",
        "structure.fp_ii",
        "structure.bp_i_a",
        "structure.bp_i_b",
        "structure.bp_ii_in",
        "structure.bp_ii_a_in",
        "structure.bp_ii_b_in",
        "structure.bp_ii_out",
        "structure.bp_ii_b",
        "structure.bp_iii",
        "structure.bp_iv",
        "structure.bp_v",
        "structure.grich_distance",
        "structure.grich_angle",
        "structure.grich_rotation",
        "structure.filepath",
    ],
    "bioactivities": [
        # "kinase.id"  # Add if added to KLIFS API? TODO
        "kinase.pref_name",
        "kinase.uniprot",
        # "ligand.id"  # Add if added to KLIFS API? TODO
        "ligand.bioactivity_standard_type",
        "ligand.bioactivity_standard_relation",
        "ligand.bioactivity_standard_value",
        "ligand.bioactivity_standard_units",
        "ligand.bioactivity_pchembl_value",
        "species.chembl",
    ],
    "interactions": ["structure.id", "interaction.fingerprint"],
    "interaction_types": ["interaction.id", "interaction.name"],
    "pockets": ["residue.klifs_id", "residue.pdb_id", "residue.klifs_region"],
    "coordinates": [
        "atom.id",
        "atom.name",
        "atom.x",
        "atom.y",
        "atom.z",
        "atom.type",
        "residue.subst_id",
        "residue.subst_name",
        "atom.charge",
        "atom.backbone",
        "residue.name",
        "residue.pdb_id",
        "residue.klifs_id",
        "residue.klifs_region",
    ],
}

POCKET_KLIFS_REGIONS = {
    1: "I.1",
    2: "I.2",
    3: "I.3",
    4: "g.l.4",
    5: "g.l.5",
    6: "g.l.6",
    7: "g.l.7",
    8: "g.l.8",
    9: "g.l.9",
    10: "II.10",
    11: "II.11",
    12: "II.12",
    13: "II.13",
    14: "III.14",
    15: "III.15",
    16: "III.16",
    17: "III.17",
    18: "III.18",
    19: "III.19",
    20: "αC.20",
    21: "αC.21",
    22: "αC.22",
    23: "αC.23",
    24: "αC.24",
    25: "αC.25",
    26: "αC.26",
    27: "αC.27",
    28: "αC.28",
    29: "αC.29",
    30: "αC.30",
    31: "b.l.31",
    32: "b.l.32",
    33: "b.l.33",
    34: "b.l.34",
    35: "b.l.35",
    36: "b.l.36",
    37: "b.l.37",
    38: "IV.38",
    39: "IV.39",
    40: "IV.40",
    41: "IV.41",
    42: "V.42",
    43: "V.43",
    44: "V.44",
    45: "GK.45",
    46: "hinge.46",
    47: "hinge.47",
    48: "hinge.48",
    49: "linker.49",
    50: "linker.50",
    51: "linker.51",
    52: "linker.52",
    53: "αD.53",
    54: "αD.54",
    55: "αD.55",
    56: "αD.56",
    57: "αD.57",
    58: "αD.58",
    59: "αD.59",
    60: "αE.60",
    61: "αE.61",
    62: "αE.62",
    63: "αE.63",
    64: "αE.64",
    65: "VI.65",
    66: "VI.66",
    67: "VI.67",
    68: "c.l.68",
    69: "c.l.69",
    70: "c.l.70",
    71: "c.l.71",
    72: "c.l.72",
    73: "c.l.73",
    74: "c.l.74",
    75: "c.l.75",
    76: "VII.76",
    77: "VII.77",
    78: "VII.78",
    79: "VIII.79",
    80: "xDFG.80",
    81: "xDFG.81",
    82: "xDFG.82",
    83: "xDFG.83",
    84: "a.l.84",
    85: "a.l.85",
}