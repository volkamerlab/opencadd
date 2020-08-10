"""
opencadd.databases.klifs.local.structures
Utility functions to work with KLIFS data (local)

Structure details.
"""

from pathlib import Path

# Redundancy here - think about alternatives.


def complex(klifs_download_path, kinase_name, pdb_id, alt=None, chain=None, species='HUMAN'):

    mol2_path = _mol2_path('complex', klifs_download_path, kinase_name, pdb_id, alt, chain, species)
    structure = _read_with_biopandas(mol2_path)

    return structure


def protein(klifs_download_path, kinase_name, pdb_id, alt=None, chain=None, species='HUMAN'):

    mol2_path = _mol2_path('protein', klifs_download_path, kinase_name, pdb_id, alt, chain, species)
    structure = _read_with_biopandas(mol2_path)

    return structure


def pocket(klifs_download_path, kinase_name, pdb_id, alt=None, chain=None, species='HUMAN'):

    mol2_path = _mol2_path('pocket', klifs_download_path, kinase_name, pdb_id, alt, chain, species)
    structure = _read_with_biopandas(mol2_path)

    return structure


def ligand(klifs_download_path, kinase_name, pdb_id, alt=None, chain=None, species='HUMAN'):

    mol2_path = _mol2_path('ligand', klifs_download_path, kinase_name, pdb_id, alt, chain, species)
    structure = _read_with_biopandas(mol2_path)

    return structure


def water(klifs_download_path, kinase_name, pdb_id, alt=None, chain=None, species='HUMAN'):

    mol2_path = _mol2_path('water', klifs_download_path, kinase_name, pdb_id, alt, chain, species)
    structure = _read_with_biopandas(mol2_path)

    return structure


def _mol2_path(structure_type, klifs_download_path, kinase_name, pdb_id, alt, chain, species):

    klifs_download_path = Path(klifs_download_path)
    species = species.upper()
    structure_name = f'{pdb_id}_alt{alt}_chain{chain}'  # Add exceptions (if alt/chain None)

    return klifs_download_path / species / kinase_name / structure_name / f'{structure_type}.mol2'


def _read_with_biopandas(mol2_path):
    pass


def _read_with_rdkit(mol2_path):
    pass
