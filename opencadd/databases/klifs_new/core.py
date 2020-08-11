"""
core.py

Defines core classes and functions.
"""


class KinasesFactory:
    """Class for kinases requests."""

    def __init__(self):
        pass

    @property
    def all_kinase_groups(self):
        """Return all kinase groups available remotely or locally."""
        pass

    def all_kinase_families(self, groups=None):
        """Return all kinase families available remotely or locally."""
        pass

    def all_kinase_names(self, groups=None, families=None, species=None):
        """Return all kinase names available remotely or locally."""
        pass

    def from_kinase_ids(self, kinase_ids):
        pass

    def from_kinase_names(self, kinase_names):
        pass


class LigandsFactory:
    """Class for ligands requests."""

    def __init__(self):
        pass

    @property
    def all(self):
        """Return all ligands available remotely or locally."""
        pass

    def from_kinase_ids(self, kinase_ids):
        pass

    def from_ligand_ids(self, ligand_ids):
        pass

    def from_kinase_names(self, kinase_names):
        pass

    def from_ligand_pdbs(self, ligand_pdbs):
        pass


class StructuresFactory:
    """Class for structures requests."""

    def __init__(self):
        pass

    @property
    def all(self):
        """Return all structures available remotely or locally."""
        pass

    def from_structure_ids(self, structure_ids):
        pass

    def from_ligand_ids(self, ligand_ids):
        pass

    def from_kinase_ids(self, kinase_ids):
        pass

    def from_structure_pdbs(self, structure_pdbs):
        pass

    def from_ligand_pdbs(self, ligand_pdbs):
        pass

    def from_kinase_names(self, kinase_names):
        pass


class BioactivitiesFactory:
    """Class for bioactivities requests."""

    def __init__(self):
        pass

    @property
    def all(self):
        """Return all bioactivities available remotely or locally."""
        pass

    def from_kinase_ids(self, kinase_ids):
        pass

    def from_ligand_ids(self, ligand_ids):
        pass


class InteractionsFactory:
    """Class for interactions requests."""

    def __init__(self):
        pass

    @property
    def interaction_types(self):
        pass

    @property
    def all(self):
        """Return all interactions available remotely or locally."""
        pass

    def from_structure_ids(self, structure_ids):
        pass

    def from_ligand_ids(self, ligand_ids):
        pass

    def from_kinase_ids(self, kinase_ids):
        pass


class CoordinatesFactory:
    """Class for coordinates requests."""

    def __init__(self):
        pass

    def from_structure_id(self, structure_id, entity, input_format, output_format):
        pass
