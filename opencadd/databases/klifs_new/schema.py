"""
schema.py

Defines schema used accross the klifs module.
"""

KINASE_GROUPS = [
    "AGC",
    "CAMK",
    "CK1",
    "CMGC",
    "Other",
    "STE",
    "TK",
    "TKL",
]

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
        "index": "structure.pocket_klifs_numbering",
        "Xray_position": "structure.pocket_pdb_numbering",
        "KLIFS_position": "structure.pocket_klifs_regions",
    },
}

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom_id", int),
        1: ("atom_name", str),
        2: ("x", float),
        3: ("y", float),
        4: ("z", float),
        5: ("atom_type", str),
        6: ("subst_id", int),
        7: ("subst_name", str),
        8: ("charge", float),
        9: ("backbone", str),
    },
    "n_cols_9": {
        0: ("atom_id", int),
        1: ("atom_name", str),
        2: ("x", float),
        3: ("y", float),
        4: ("z", float),
        5: ("atom_type", str),
        6: ("subst_id", int),
        7: ("subst_name", str),
        8: ("charge", float),
    },
}

LOCAL_REMOTE_COLUMNS = {
    "kinase_groups": {"local": ["kinase.group"], "remote": ["kinase.group"]},
    "kinase_families": {"local": ["kinase.family"], "remote": ["kinase.family"]},
    "kinases_all": {
        "local": ["kinase.name", "species.klifs"],
        "remote": ["kinase.id", "kinase.name", "kinase.name_full", "species.klifs",],
    },
    "kinases": {
        "local": [
            "kinase.name",
            "kinase.group",
            "kinase.family",
            "species.klifs",
            "kinase.pocket",
        ],
        "remote": list(REMOTE_COLUMNS_MAPPING["kinases"].values()),
    },
    "ligands": {
        "local": ["ligand.pdb", "ligand.name"],
        "remote": list(REMOTE_COLUMNS_MAPPING["ligands"].values()),
    },
    "structures": {
        "local": [
            "kinase.name",
            # "kinase.name_all",
            "kinase.family",
            "kinase.group",
            "species.klifs",
            "structure.pdb",
            "structure.alternate_model",
            "structure.chain",
            "structure.rmsd1",
            "structure.rmsd2",
            "kinase.pocket",
            "structure.resolution",
            "structure.qualityscore",
            "structure.missing_residues",
            "structure.missing_atoms",
            "ligand.pdb",
            "ligand.pdb_allosteric",
            "ligand.name",
            "ligand.name_allosteric",
            "structure.dfg",
            "structure.ac_helix",
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
            "structure.filepath",
        ],
        "remote": list(REMOTE_COLUMNS_MAPPING["structures"].values()),
    },
    "bioactivities": {
        "local": [],  # Not available locally
        "remote": list(REMOTE_COLUMNS_MAPPING["bioactivities"].values()),
    },
    "interactions": {
        "local": [
            "structure.pdb",
            "structure.alternate_model",
            "structure.chain",
            "interaction.fingerprint",
        ],
        "remote": list(REMOTE_COLUMNS_MAPPING["interactions"].values()),
    },
    "interaction_types": {
        "local": [],  # Not available locally
        "remote": list(REMOTE_COLUMNS_MAPPING["interaction_types"].values()),
    },
    "pockets": {
        "local": [],  # Not available locally
        "remote": list(REMOTE_COLUMNS_MAPPING["pockets"].values()),
    },
    "coordinates": {"local": [], "remote": [],},
}
