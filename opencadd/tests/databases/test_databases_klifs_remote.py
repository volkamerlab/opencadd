"""
Tests for opencadd.databases.klifs.remote
"""

from bravado_core.exception import SwaggerMappingError
import pandas as pd
import pytest

from opencadd.databases.klifs_new.api import setup_remote
from opencadd.databases.klifs_new.schema import KINASE_GROUPS, LOCAL_REMOTE_COLUMNS


class TestsAllQueries:
    """
    Test all class methods that return all entries for a query entity such as 
    kinases, ligands, structure, etc.
    """

    def test_all_kinase_groups(self):
        """
        Test request result for all kinase groups.
        """
        session = setup_remote()

        kinases = session.kinases.all_kinase_groups()
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_groups"]
        assert kinases.shape[0] == 8
        assert sorted(kinases["kinase.group"].to_list()) == KINASE_GROUPS

    @pytest.mark.parametrize("group", [None, *KINASE_GROUPS])
    def test_all_kinase_families(self, group):
        """
        Test request result for all kinase families.
        """
        session = setup_remote()

        kinases = session.kinases.all_kinase_families(group)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_families"]

    @pytest.mark.parametrize("group", ["XXX"])
    def test_all_kinase_families_raise(self, group):
        """
        Test request result for all kinase families: Error raised if input invalid?
        """
        session = setup_remote()

        with pytest.raises(SwaggerMappingError):
            session.kinases.all_kinase_families(group)

    @pytest.mark.parametrize(
        "group, family, species",
        [
            (None, None, None),
            ("TK", None, None),
            (None, "EGFR", None),
            (None, None, "Human"),
            (None, None, "HUMAN"),
            ("TK", "EGFR", "Human"),
        ],
    )
    def test_all_kinases(self, group, family, species):
        """
        Test request result for all kinases.
        """
        session = setup_remote()

        kinases = session.kinases.all_kinases(group, family, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases_all"]["remote"]

    @pytest.mark.parametrize(
        "group, family, species", [("XXX", None, None), (None, "XXX", None), ("XXX", None, "XXX")],
    )
    def test_all_kinases_raise(self, group, family, species):
        """
        Test request result for all kinases: Error raised if input invalid?
        """
        session = setup_remote()
        with pytest.raises(SwaggerMappingError):
            session.kinases.all_kinases(group, family, species)

    def test_all_ligands(self):
        """
        Test request result for all ligands.
        """
        session = setup_remote()

        ligands = session.ligands.all_ligands()
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

    def test_all_structures(self):
        """
        Test request result for all kinases.
        """
        session = setup_remote()

        structures = session.structures.all_structures()
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    def test_interaction_types(self):
        """
        Test request result for all interaction types.
        """
        session = setup_remote()

        interaction_types = session.interactions.interaction_types
        assert isinstance(interaction_types, pd.DataFrame)
        assert (
            interaction_types.columns.to_list()
            == LOCAL_REMOTE_COLUMNS["interaction_types"]["remote"]
        )

    def test_all_interactions(self):
        """
        Test request result for all kinases.
        """
        session = setup_remote()

        interactions = session.interactions.all_interactions()
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

    def test_all_bioactivities(self):

        """
        Test request result for all kinases.
        """
        session = setup_remote()

        # Usually this class method is used to get ALL bioactivities from ALL ligands
        # Since this query takes a few minutes, only the frist 3 ligands are used here for testing
        bioactivities = session.bioactivities.all_bioactivities(n=3)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"][
            "remote"
        ] + ["ligand.id (query)"]


class TestsFromKinaseIds:
    """
    Test all class methods with kinase IDs as input.
    """

    @pytest.mark.parametrize("kinase_ids", [1, [1, 2], [1, 2, 10000]])
    def test_from_kinase_ids(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input.
        """
        session = setup_remote()

        # Kinases
        kinases = session.kinases.from_kinase_ids(kinase_ids)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["remote"]

        # Ligands
        ligands = session.ligands.from_kinase_ids(kinase_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"] + [
            "kinase.id (query)"
        ]

        # Structures
        structures = session.structures.from_kinase_ids(kinase_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

        # Bioactivities
        bioactivities = session.bioactivities.from_kinase_ids(kinase_ids)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"][
            "remote"
        ] + ["ligand.id (query)"]

        # Interactions
        interactions = session.interactions.from_kinase_ids(kinase_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert (
            interactions.columns.to_list()
            == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]  # TODO Add "kinase.id (query)"?
        )

    @pytest.mark.parametrize("kinase_ids", [10000, "XXX"])
    def test_from_kinase_ids_raise(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input: Error raised if input invalid?
        """
        session = setup_remote()

        with pytest.raises(SwaggerMappingError):
            session.kinases.from_kinase_ids(kinase_ids)
            session.ligands.from_kinase_ids(kinase_ids)
            session.structures.from_kinase_ids(kinase_ids)
            session.bioactivities.from_kinase_ids(kinase_ids)
            session.interactions.from_kinase_ids(kinase_ids)


class TestsFromLigandIds:
    """
    Test all class methods with ligand IDs as input.
    """

    @pytest.mark.parametrize("ligand_ids", [10, [10, 20], [10, 20, 10000]])
    def test_from_ligand_ids(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input.
        """
        session = setup_remote()

        # Ligands
        ligands = session.ligands.from_ligand_ids(ligand_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

        # Structures
        structures = session.structures.from_ligand_ids(ligand_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert (
            structures.columns.to_list()
            == LOCAL_REMOTE_COLUMNS["structures"]["remote"]  # TODO Add "ligand.id (query)"?
        )

        # Bioactivities
        bioactivities = session.bioactivities.from_ligand_ids(ligand_ids)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"][
            "remote"
        ] + ["ligand.id (query)"]

        # Interactions
        interactions = session.interactions.from_ligand_ids(ligand_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert (
            interactions.columns.to_list()
            == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]  # TODO Add "ligand.id (query)"?
        )

    @pytest.mark.parametrize("ligand_ids", [10000, "XXX"])
    def test_from_ligand_ids_raise(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input: Error raised if input invalid?
        """
        session = setup_remote()

        with pytest.raises(SwaggerMappingError):
            session.ligands.from_ligand_ids(ligand_ids)
            session.structures.from_ligand_ids(ligand_ids)
            session.bioactivities.from_ligand_ids(ligand_ids)
            session.interactions.from_ligand_ids(ligand_ids)


class TestsFromStructureIds:
    """
    Test class methods with structure IDs as input.
    """

    @pytest.mark.parametrize("structure_ids", [100, [100, 100000]])
    def test_from_structure_ids(self, structure_ids):
        """
        Test class methods with structure IDs as input.
        """
        session = setup_remote()

        # Structures
        structures = session.structures.from_structure_ids(structure_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

        # Interactions
        interactions = session.interactions.from_structure_ids(structure_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

        # Pockets (takes only one structure ID as input!)
        if isinstance(structure_ids, int):
            structure_id = structure_ids
            pocket = session.pockets.from_structure_id(structure_id)
            assert isinstance(pocket, pd.DataFrame)
            assert pocket.columns.to_list() == LOCAL_REMOTE_COLUMNS["pockets"]["remote"]

    @pytest.mark.parametrize("structure_ids", [100000, "XXX"])
    def test_from_structure_ids_raise(self, structure_ids):
        """
        Test class methods with structure IDs as input: Error raised if input invalid?
        """
        session = setup_remote()

        with pytest.raises(SwaggerMappingError):
            session.structures.from_structure_ids(structure_ids)
            session.interactions.from_structure_ids(structure_ids)


class TestsFromKinaseNames:
    """
    Test class methods with kinase names as input.
    """

    @pytest.mark.parametrize(
        "kinase_names, species",
        [("EGFR", None), (["EGFR", "BMX"], None), (["EGFR", "XXX"], None)],
    )
    def test_from_kinase_names(self, kinase_names, species):
        """
        Test class methods with kinase names as input.
        """
        session = setup_remote()

        # Kinases
        kinases = session.kinases.from_kinase_names(kinase_names, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["remote"]

        # Ligands
        ligands = session.ligands.from_kinase_names(kinase_names)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"] + [
            "kinase.id (query)",
            "kinase.name (query)",
            "species.klifs (query)",
        ]

        # Structures
        structures = session.structures.from_kinase_names(kinase_names)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("kinase_names, species", [("XXX", None), (1, None), ("EGFR", "XXX")])
    def test_from_kinase_names_raise(self, kinase_names, species):
        """
        Test class methods with kinase names as input: Error raised if input invalid?
        """
        session = setup_remote()
        with pytest.raises(SwaggerMappingError):
            session.kinases.from_kinase_names(kinase_names, species)
            session.ligands.from_kinase_names(kinase_names)
            session.structures.from_kinase_names(kinase_names)


class TestsFromLigandPdbs:
    """
    Test class methods with ligand PDB IDs as input.
    """

    @pytest.mark.parametrize("ligand_pdbs", ["STU", ["STU", "STI"], ["STU", "XXX"]])
    def test_from_ligand_pdbs(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input.
        """
        session = setup_remote()

        # Ligands
        ligands = session.ligands.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

        # Structure
        structures = session.structures.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("ligand_pdbs", [1, "XXX"])
    def test_from_ligand_pdbs_raise(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input: Error raised if input invalid?
        """
        session = setup_remote()

        with pytest.raises(SwaggerMappingError):
            session.ligands.from_ligand_pdbs(ligand_pdbs)
            session.structures.from_ligand_pdbs(ligand_pdbs)


class TestsFromStructurePdbs:
    """
    Test class methods with structure PDB IDs as input.
    """

    @pytest.mark.parametrize("structure_pdbs", ["3sxr", ["3sxr", "1fpu", "xxxx"]])
    def test_from_structure_pdbs(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input.
        """
        session = setup_remote()

        # Structure
        structures = session.structures.from_structure_pdbs(structure_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("structure_pdbs", [1, "xxxx"])
    def test_from_structure_pdbs_raise(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input: Error raised if input invalid?
        """
        session = setup_remote()

        # Structure
        # Note: Integer input raises jsonschema.exceptions.ValidationError,
        # here tested using Exception
        with pytest.raises((SwaggerMappingError, Exception)):
            session.structures.from_structure_pdbs(structure_pdbs)


class TestsCoordinates:
    """
    Test remote Coordinates class.
    """

    def ttest_fetch(self):
        pass

    def ttest_save(self):
        pass
