"""
Tests for opencadd.databases.klifs.remote
"""

from bravado_core.exception import SwaggerMappingError
import pandas as pd
import pytest
from rdkit import Chem

from opencadd.databases.klifs.api import setup_remote
from opencadd.databases.klifs.schema import KINASE_GROUPS, LOCAL_REMOTE_COLUMNS


class TestsAllQueries:
    """
    Test all class methods that return all entries for a query entity such as 
    kinases, ligands, structure, etc.
    """

    def test_all_kinase_groups(self):
        """
        Test request result for all kinase groups.
        """
        remote = setup_remote()

        kinases = remote.kinases.all_kinase_groups()
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_groups"]["remote"]
        assert sorted(kinases["kinase.group"].to_list()) == KINASE_GROUPS

    @pytest.mark.parametrize("group", [None, *KINASE_GROUPS])
    def test_all_kinase_families(self, group):
        """
        Test request result for all kinase families.
        """
        remote = setup_remote()

        kinases = remote.kinases.all_kinase_families(group)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinase_families"]["remote"]

    @pytest.mark.parametrize("group", ["XXX"])
    def test_all_kinase_families_raise(self, group):
        """
        Test request result for all kinase families: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(SwaggerMappingError):
            remote.kinases.all_kinase_families(group)

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
        remote = setup_remote()

        kinases = remote.kinases.all_kinases(group, family, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases_all"]["remote"]

    @pytest.mark.parametrize(
        "group, family, species", [("XXX", None, None), (None, "XXX", None), ("XXX", None, "XXX")],
    )
    def test_all_kinases_raise(self, group, family, species):
        """
        Test request result for all kinases: Error raised if input invalid?
        """
        remote = setup_remote()
        with pytest.raises(SwaggerMappingError):
            remote.kinases.all_kinases(group, family, species)

    def test_all_ligands(self):
        """
        Test request result for all ligands.
        """
        remote = setup_remote()

        ligands = remote.ligands.all_ligands()
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

    def test_all_structures(self):
        """
        Test request result for all structures.
        """
        remote = setup_remote()

        structures = remote.structures.all_structures()
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    def test_interaction_types(self):
        """
        Test request result for all interaction types.
        """
        remote = setup_remote()

        interaction_types = remote.interactions.interaction_types
        assert isinstance(interaction_types, pd.DataFrame)
        assert (
            interaction_types.columns.to_list()
            == LOCAL_REMOTE_COLUMNS["interaction_types"]["remote"]
        )

    def test_all_interactions(self):
        """
        Test request result for all kinases.
        """
        remote = setup_remote()

        interactions = remote.interactions.all_interactions()
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

    def test_all_bioactivities(self):

        """
        Test request result for all kinases.
        """
        remote = setup_remote()

        # Usually this class method is used to get ALL bioactivities from ALL ligands
        # Since this query takes a few minutes, only the frist 3 ligands are used here for testing
        bioactivities = remote.bioactivities.all_bioactivities(n=3)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"]["remote"]


class TestsFromKinaseIds:
    """
    Test all class methods with kinase IDs as input.
    """

    @pytest.mark.parametrize("kinase_ids", [1, [1, 2], [1, 2, 10000]])
    def test_from_kinase_ids(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input.
        """
        remote = setup_remote()

        # Kinases
        kinases = remote.kinases.from_kinase_ids(kinase_ids)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["remote"]

        # Ligands
        ligands = remote.ligands.from_kinase_ids(kinase_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"] + [
            "kinase.id (query)"
        ]

        # Structures
        structures = remote.structures.from_kinase_ids(kinase_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

        # Bioactivities
        bioactivities = remote.bioactivities.from_kinase_ids(kinase_ids)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"]["remote"]

        # Interactions
        interactions = remote.interactions.from_kinase_ids(kinase_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

    @pytest.mark.parametrize("kinase_ids", [10000, "XXX"])
    def test_from_kinase_ids_raise(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(SwaggerMappingError):
            remote.kinases.from_kinase_ids(kinase_ids)
            remote.ligands.from_kinase_ids(kinase_ids)
            remote.structures.from_kinase_ids(kinase_ids)
            remote.bioactivities.from_kinase_ids(kinase_ids)
            remote.interactions.from_kinase_ids(kinase_ids)


class TestsFromLigandIds:
    """
    Test all class methods with ligand IDs as input.
    """

    @pytest.mark.parametrize("ligand_ids", [10, [10, 20], [10, 20, 10000]])
    def test_from_ligand_ids(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input.
        """
        remote = setup_remote()

        # Ligands
        ligands = remote.ligands.from_ligand_ids(ligand_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

        # Structures
        structures = remote.structures.from_ligand_ids(ligand_ids)
        assert isinstance(ligands, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

        # Bioactivities
        bioactivities = remote.bioactivities.from_ligand_ids(ligand_ids)
        assert isinstance(bioactivities, pd.DataFrame)
        assert bioactivities.columns.to_list() == LOCAL_REMOTE_COLUMNS["bioactivities"]["remote"]

        # Interactions
        interactions = remote.interactions.from_ligand_ids(ligand_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

    @pytest.mark.parametrize("ligand_ids", [10000, "XXX"])
    def test_from_ligand_ids_raise(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(SwaggerMappingError):
            remote.bioactivities.from_ligand_ids(ligand_ids)
            remote.interactions.from_ligand_ids(ligand_ids)

        with pytest.raises(ValueError):
            remote.ligands.from_ligand_ids(ligand_ids)
            remote.structures.from_ligand_ids(ligand_ids)


class TestsFromStructureIds:
    """
    Test class methods with structure IDs as input.
    """

    @pytest.mark.parametrize("structure_ids", [12347, [12347, 100000]])
    def test_from_structure_ids(self, structure_ids):
        """
        Test class methods with structure IDs as input.
        """
        remote = setup_remote()

        # Structures
        structures = remote.structures.from_structure_ids(structure_ids)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

        # Interactions
        interactions = remote.interactions.from_structure_ids(structure_ids)
        assert isinstance(interactions, pd.DataFrame)
        assert interactions.columns.to_list() == LOCAL_REMOTE_COLUMNS["interactions"]["remote"]

        # Pockets (takes only one structure ID as input!)
        if isinstance(structure_ids, int):
            structure_id = structure_ids
            pocket = remote.pockets.from_structure_id(structure_id)
            assert isinstance(pocket, pd.DataFrame)
            assert pocket.columns.to_list() == LOCAL_REMOTE_COLUMNS["pockets"]["remote"]

    @pytest.mark.parametrize("structure_ids", [100000, "XXX"])
    def test_from_structure_ids_raise(self, structure_ids):
        """
        Test class methods with structure IDs as input: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(SwaggerMappingError):
            remote.structures.from_structure_ids(structure_ids)
            remote.interactions.from_structure_ids(structure_ids)
            if isinstance(structure_ids, int):
                structure_id = structure_ids
                remote.pockets.from_structure_id(structure_id)


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
        remote = setup_remote()

        # Kinases
        kinases = remote.kinases.from_kinase_names(kinase_names, species)
        assert isinstance(kinases, pd.DataFrame)
        assert kinases.columns.to_list() == LOCAL_REMOTE_COLUMNS["kinases"]["remote"]

        # Ligands
        ligands = remote.ligands.from_kinase_names(kinase_names)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"] + [
            "kinase.id (query)",
            "kinase.name (query)",
            "species.klifs (query)",
        ]

        # Structures
        structures = remote.structures.from_kinase_names(kinase_names)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("kinase_names, species", [("XXX", None), (1, None), ("XXX", "XXX")])
    def test_from_kinase_names_raise(self, kinase_names, species):
        """
        Test class methods with kinase names as input: Error raised if input invalid?
        """
        remote = setup_remote()
        with pytest.raises(SwaggerMappingError):
            remote.kinases.from_kinase_names(kinase_names, species)
            remote.ligands.from_kinase_names(kinase_names)

        with pytest.raises(ValueError):
            remote.structures.from_kinase_names(kinase_names)


class TestsFromLigandPdbs:
    """
    Test class methods with ligand PDB IDs as input.
    """

    @pytest.mark.parametrize("ligand_pdbs", ["PRC", ["PRC", "1N1"], ["PRC", "1N1", "XXX"]])
    def test_from_ligand_pdbs(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input.
        """
        remote = setup_remote()

        # Ligands
        ligands = remote.ligands.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(ligands, pd.DataFrame)
        assert ligands.columns.to_list() == LOCAL_REMOTE_COLUMNS["ligands"]["remote"]

        # Structure
        structures = remote.structures.from_ligand_pdbs(ligand_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("ligand_pdbs", [1, "XXX"])
    def test_from_ligand_pdbs_raise(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(ValueError):
            remote.ligands.from_ligand_pdbs(ligand_pdbs)
            remote.structures.from_ligand_pdbs(ligand_pdbs)


class TestsFromStructurePdbs:
    """
    Test class methods with structure PDB IDs as input.
    """

    @pytest.mark.parametrize("structure_pdbs", ["3sxr", ["3sxr", "1fpu", "xxxx"]])
    def test_from_structure_pdbs(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input.
        """
        remote = setup_remote()

        # Structure
        structures = remote.structures.from_structure_pdbs(structure_pdbs)
        assert isinstance(structures, pd.DataFrame)
        assert structures.columns.to_list() == LOCAL_REMOTE_COLUMNS["structures"]["remote"]

    @pytest.mark.parametrize("structure_pdbs", [1, "xxxx"])
    def test_from_structure_pdbs_raise(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input: Error raised if input invalid?
        """
        remote = setup_remote()

        # Structure
        # Note: Integer input raises jsonschema.exceptions.ValidationError,
        # here tested using Exception
        with pytest.raises((SwaggerMappingError, Exception)):
            remote.structures.from_structure_pdbs(structure_pdbs)


class TestsCoordinates:
    """
    Test remote Coordinates class.
    """

    @pytest.mark.parametrize(
        "structure_id, entity, input_format, output_format, n_atoms, centroid",
        [
            (12347, "complex", "mol2", "biopandas", 3604, [-3.996449, 17.509910, 31.077763]),
            (12347, "ligand", "mol2", "biopandas", 49, [2.291216, 20.590290, 39.074586]),
            (12347, "ligand", "mol2", "rdkit", None, None),
            (12347, "pocket", "mol2", "biopandas", 1156, [0.308657, 21.880768, 35.903844]),
            (12347, "protein", "mol2", "biopandas", 3552, [-4.061604, 17.472405, 30.976719]),
        ],
    )
    def test_from_structure_id(
        self, structure_id, entity, input_format, output_format, n_atoms, centroid
    ):
        """
        Test retrieval of coordinates data from structure ID.
        """
        remote = setup_remote()

        coordinates = remote.coordinates.from_structure_id(
            structure_id, entity, input_format, output_format
        )
        if output_format == "biopandas":
            assert isinstance(coordinates, pd.DataFrame)
            assert coordinates.columns.to_list() == LOCAL_REMOTE_COLUMNS["coordinates"]
            assert coordinates.shape[0] == n_atoms
            assert pytest.approx(coordinates["atom.x"].mean(), centroid[0], abs=1.0e-6)
            assert pytest.approx(coordinates["atom.y"].mean(), centroid[1], abs=1.0e-6)
            assert pytest.approx(coordinates["atom.z"].mean(), centroid[2], abs=1.0e-6)
        elif output_format == "rdkit":
            assert isinstance(coordinates, Chem.rdchem.Mol)

    @pytest.mark.parametrize(
        "structure_id, entity, input_format, output_format",
        [
            (12347, "complex", "mol2", "rdkit"),
            (12347, "ligand", "pdb", "biopandas"),
            (12347, "water", "mol2", "biopandas"),
        ],
    )
    def test_from_structure_id_raise(self, structure_id, entity, input_format, output_format):
        """
        Test retrieval of coordinates data from structure ID: Error raised if input invalid?
        """
        remote = setup_remote()

        with pytest.raises(ValueError):
            remote.coordinates.from_structure_id(structure_id, entity, input_format, output_format)

