"""
Tests for opencadd.databases.klifs.local
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.databases.klifs.api import setup_local
from opencadd.databases.klifs.schema import LOCAL_REMOTE_COLUMNS

PATH_TEST_DATA = Path(__file__).parent / "data" / "KLIFS_download"


class TestsAllQueries:
    """
    Test local Kinases class.
    """

    @pytest.mark.parametrize("result_groups", [["TK", "TKL"]])
    def test_all_kinase_groups(self, result_groups):
        """
        Test request result for all kinase groups.
        """
        session = setup_local(PATH_TEST_DATA)

        kinases = session.kinases.all_kinase_groups()
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_groups"]["local"]
        assert kinases["kinase.group"].to_list() == result_groups

    @pytest.mark.parametrize(
        "group, result_families", [(None, ["Tec", "RAF", "Abl"]), ("TK", ["Tec", "Abl"])],
    )
    def test_all_kinase_families(self, group, result_families):
        """
        Test request result for all kinase families.
        """
        session = setup_local(PATH_TEST_DATA)

        kinases = session.kinases.all_kinase_families(group)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_families"]["local"]
        assert kinases["kinase.family"].to_list() == result_families

    @pytest.mark.parametrize("group", ["XXX"])
    def test_all_kinase_families_raise(self, group):
        """
        Test request result for all kinase families: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)

        with pytest.raises(ValueError):
            session.kinases.all_kinase_families(group)

    @pytest.mark.parametrize(
        "group, family, species, result_kinases",
        [
            (None, None, None, [["BMX", "BRAF", "ABL1"], ["Human", "Human", "Mouse"]]),
            ("TK", None, None, [["BMX", "ABL1"], ["Human", "Mouse"]]),
            (None, "Tec", None, [["BMX"], ["Human"]]),
            (None, None, "Human", [["BMX", "BRAF"], ["Human", "Human"]]),
            (None, None, "HUMAN", [["BMX", "BRAF"], ["Human", "Human"]]),
            ("TK", "Tec", "Human", [["BMX"], ["Human"]]),
        ],
    )
    def test_all_kinases(self, group, family, species, result_kinases):
        """
        Test request result for all kinases.
        """
        session = setup_local(PATH_TEST_DATA)

        kinases = session.kinases.all_kinases(group, family, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases_all"]["local"]

    @pytest.mark.parametrize(
        "group, family, species", [("XXX", None, None), (None, "XXX", None), ("XXX", None, "XXX")],
    )
    def test_all_kinases_raise(self, group, family, species):
        """
        Test request result for all kinases: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)

        with pytest.raises(ValueError):
            session.kinases.all_kinases(group, family, species)

    def test_all_ligands(self):
        """
        Test request result for all ligands.
        """
        session = setup_local(PATH_TEST_DATA)

        ligands = session.ligands.all_ligands()
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["local"]
        assert ligands["ligand.pdb"].to_list() == ["1N1", "QH1", "PRC"]

    def test_all_structures(self):
        """
        Test request result for all structures.
        """
        session = setup_local(PATH_TEST_DATA)

        structures = session.structures.all_structures()
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]
        assert structures["structure.id"].to_list() == [3482, 12347, 5728, 5705]

    def test_interaction_types(self):
        """
        Test request result for all interaction types.
        """
        session = setup_local(PATH_TEST_DATA)

        with pytest.raises(NotImplementedError):
            session.interactions.interaction_types

    def test_all_interactions(self):
        """
        Test request result for all kinases.
        """
        session = setup_local(PATH_TEST_DATA)

        interactions = session.interactions.all_interactions()
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["local"]

    def test_all_bioactivities(self):

        """
        Test request result for all kinases.
        """
        session = setup_local(PATH_TEST_DATA)
        with pytest.raises(NotImplementedError):
            session.bioactivities.all_bioactivities()


class TestsFromKinaseIds:
    """
    Test all class methods with kinase IDs as input.
    """

    @pytest.mark.parametrize("kinase_ids", [472, [472, 509], [472, 509, 10000]])
    def test_from_kinase_ids(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input.
        """
        session = setup_local(PATH_TEST_DATA)

        # Kinases
        kinases = session.kinases.from_kinase_ids(kinase_ids)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["local"]

        # Ligands
        ligands = session.ligands.from_kinase_ids(kinase_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["local"] + [
            "kinase.id (query)"
        ]

        # Structures
        structures = session.structures.from_kinase_ids(kinase_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]

        # Bioactivities
        with pytest.raises(NotImplementedError):
            session.bioactivities.from_kinase_ids(kinase_ids)

        # Interactions
        interactions = session.interactions.from_kinase_ids(kinase_ids)
        assert isinstance(kinases, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["local"] + [
            "kinase.id (query)"
        ]

    @pytest.mark.parametrize("kinase_ids", [10000, "XXX"])
    def test_from_kinase_ids_raise(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)
        with pytest.raises(ValueError):
            session.kinases.from_kinase_ids(kinase_ids)
            session.ligands.from_kinase_ids(kinase_ids)
            session.structures.from_kinase_ids(kinase_ids)
            session.interactions.from_kinase_ids(kinase_ids)


class TestFromStructureIds:
    """
    Test class methods with structure IDs as input.
    """

    @pytest.mark.parametrize("structure_ids", [12347, [12347, 100000]])
    def test_from_structure_ids(self, structure_ids):
        """
        Test class methods with structure IDs as input.
        """
        session = setup_local(PATH_TEST_DATA)

        # Structures
        structures = session.structures.from_structure_ids(structure_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]

        # Interactions
        interactions = session.interactions.from_structure_ids(structure_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["local"]

        # Pockets (takes only one structure ID as input!)
        if isinstance(structure_ids, int):
            structure_id = structure_ids
            pocket = session.pockets.from_structure_id(structure_id)
            assert isinstance(pocket, pd.DataFrame)
            assert pocket.columns.to_list() == LOCAL_REMOTE_COLUMNS["pockets"]["local"]

    @pytest.mark.parametrize("structure_ids", [100000, "XXX"])
    def test_from_structure_ids_raise(self, structure_ids):
        """
        Test class methods with structure IDs as input: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)

        with pytest.raises(ValueError):
            session.structures.from_structure_ids(structure_ids)
            session.interactions.from_structure_ids(structure_ids)
            if isinstance(structure_ids, int):
                structure_id = structure_ids
                session.pockets.from_structure_id(structure_id)


class TestsFromKinaseNames:
    """
    Test class methods with kinase names as input.
    """

    @pytest.mark.parametrize("kinase_names, species", [("BMX", None), (["BMX", "BRAF"], None)])
    def test_from_kinase_names(self, kinase_names, species):
        """
        Test class methods with kinase names as input.
        """
        session = setup_local(PATH_TEST_DATA)

        # Kinases
        kinases = session.kinases.from_kinase_names(kinase_names, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["local"]

        # Ligands
        ligands = session.ligands.from_kinase_names(kinase_names)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["local"] + [
            "kinase.name (query)",
            "species.klifs (query)",
        ]

        # Structures
        structures = session.structures.from_kinase_names(kinase_names)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]

    @pytest.mark.parametrize("kinase_names, species", [("XXX", None), (1, None), ("EGFR", "XXX")])
    def test_from_kinase_names_raise(self, kinase_names, species):
        """
        Test class methods with kinase names as input: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)
        with pytest.raises(ValueError):
            session.kinases.from_kinase_names(kinase_names, species)
            session.ligands.from_kinase_names(kinase_names)
            session.structures.from_kinase_names(kinase_names)


class TestsFromLigandPdbs:
    """
    Test class methods with ligand PDB IDs as input.
    """

    @pytest.mark.parametrize("ligand_pdbs", ["PRC", ["PRC", "1N1"], ["PRC", "1N1", "XXX"]])
    def test_from_ligand_pdbs(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input.
        """
        session = setup_local(PATH_TEST_DATA)

        # Ligands
        ligands = session.ligands.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["local"]

        # Structure
        structures = session.structures.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]

    @pytest.mark.parametrize("ligand_pdbs", [1, "XXX"])
    def test_from_ligand_pdbs_raise(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)

        with pytest.raises(ValueError):
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
        session = setup_local(PATH_TEST_DATA)

        # Structure
        structures = session.structures.from_structure_pdbs(structure_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["local"]

    @pytest.mark.parametrize("structure_pdbs", [1, "xxxx"])
    def test_from_structure_pdbs_raise(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input: Error raised if input invalid?
        """
        session = setup_local(PATH_TEST_DATA)

        # Structure
        with pytest.raises(ValueError):
            session.structures.from_structure_pdbs(structure_pdbs)
