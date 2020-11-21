"""
Tests for opencadd.databases.klifs.remote
"""

from pathlib import Path

from bravado_core.exception import SwaggerMappingError
import pandas as pd
import pytest
from rdkit import Chem

from opencadd.databases.klifs.api import setup_local, setup_remote
from opencadd.databases.klifs.schema import COLUMN_NAMES
from opencadd.utils import enter_temp_directory

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/klifs"

# Set local and remote session
REMOTE = setup_remote()
LOCAL = setup_local(PATH_TEST_DATA)


def check_dataframe(dataframe, column_names):
    """
    Base function that tests if input is a DataFrame, if column names and index is correct.
    """
    # Is input a DataFrame?
    assert isinstance(dataframe, pd.DataFrame)

    # Are DataFrame column names and their order correct?
    assert dataframe.columns.to_list() == column_names

    # Are DataFrame indices enumerated starting from 0 to length of DataFrame - 1?
    assert dataframe.index.to_list() == list(range(0, len(dataframe)))


class TestsAllQueries:
    """
    Test all class methods that return all entries for a query entity such as
    kinases, ligands, structure, etc.
    """

    def test_all_kinase_groups(self):
        """
        Test request result for all kinase groups.
        """

        result_remote = REMOTE.kinases.all_kinase_groups()
        result_local = LOCAL.kinases.all_kinase_groups()

        check_dataframe(result_remote, COLUMN_NAMES["kinase_groups"])
        check_dataframe(result_local, COLUMN_NAMES["kinase_groups"])

        assert sorted(result_remote["kinase.group"].to_list()) == [
            "AGC",
            "CAMK",
            "CK1",
            "CMGC",
            "Other",
            "STE",
            "TK",
            "TKL",
        ]
        assert sorted(result_local["kinase.group"].to_list()) == ["TK", "TKL"]

    @pytest.mark.parametrize(
        "group, local_families",
        [(None, ["Tec", "RAF", "Abl"]), ("TK", ["Tec", "Abl"])],
    )
    def test_all_kinase_families(self, group, local_families):
        """
        Test request result for all kinase families.
        """

        result_remote = REMOTE.kinases.all_kinase_families(group)
        result_local = LOCAL.kinases.all_kinase_families(group)

        check_dataframe(result_remote, COLUMN_NAMES["kinase_families"])
        check_dataframe(result_local, COLUMN_NAMES["kinase_families"])

        assert result_local["kinase.family"].to_list() == local_families
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    @pytest.mark.parametrize("group", ["XXX"])
    def test_all_kinase_families_raise(self, group):
        """
        Test request result for all kinase families: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.all_kinase_families(group)
        with pytest.raises(ValueError):
            LOCAL.kinases.all_kinase_families(group)

    @pytest.mark.parametrize(
        "group, family, species, local_kinases",
        [
            (None, None, None, [["BMX", "BRAF", "Abl1"], ["Human", "Human", "Mouse"]]),
            ("TK", None, None, [["BMX", "Abl1"], ["Human", "Mouse"]]),
            (None, "Tec", None, [["BMX"], ["Human"]]),
            (None, None, "Human", [["BMX", "BRAF"], ["Human", "Human"]]),
            (None, None, "HUMAN", [["BMX", "BRAF"], ["Human", "Human"]]),
            ("TK", "Tec", "Human", [["BMX"], ["Human"]]),
        ],
    )
    def test_all_kinases(self, group, family, species, local_kinases):
        """
        Test request result for all kinases.
        """

        result_remote = REMOTE.kinases.all_kinases(group, family, species)
        result_local = LOCAL.kinases.all_kinases(group, family, species)

        check_dataframe(result_remote, COLUMN_NAMES["kinases_all"])
        check_dataframe(result_local, COLUMN_NAMES["kinases_all"])

        assert result_local["kinase.gene_name"].to_list() == local_kinases[0]
        assert result_local["species.klifs"].to_list() == local_kinases[1]
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    @pytest.mark.parametrize(
        "group, family, species",
        [("XXX", None, None), (None, "XXX", None), ("XXX", None, "XXX")],
    )
    def test_all_kinases_raise(self, group, family, species):
        """
        Test request result for all kinases: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.all_kinases(group, family, species)
        with pytest.raises(ValueError):
            LOCAL.kinases.all_kinases(group, family, species)

    def test_all_ligands(self):
        """
        Test request result for all ligands.
        """

        result_remote = REMOTE.ligands.all_ligands()
        result_local = LOCAL.ligands.all_ligands()

        check_dataframe(result_remote, COLUMN_NAMES["ligands"])
        check_dataframe(result_local, COLUMN_NAMES["ligands"])

        assert result_local["ligand.expo_id"].to_list() == ["1N1", "QH1", "PRC"]
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    def test_all_structures(self):
        """
        Test request result for all structures.
        """

        result_remote = REMOTE.structures.all_structures()
        result_local = LOCAL.structures.all_structures()

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

        assert result_local["structure.klifs_id"].to_list() == [3482, 12347, 5728, 5705]
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    def test_interaction_types(self):
        """
        Test request result for all interaction types.
        """

        result_remote = REMOTE.interactions.interaction_types
        check_dataframe(result_remote, COLUMN_NAMES["interaction_types"])

        with pytest.raises(NotImplementedError):
            LOCAL.interactions.interaction_types()

    def test_all_interactions(self):
        """
        Test request result for all kinases.
        """

        result_remote = REMOTE.interactions.all_interactions()
        result_local = LOCAL.interactions.all_interactions()

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])
        check_dataframe(result_local, COLUMN_NAMES["interactions"])

    def test_all_bioactivities(self):

        """
        Test request result for all kinases.
        """

        # Usually this class method is used to get ALL bioactivities from ALL ligands
        # Since this query takes a few minutes, only the frist 3 ligands are used here for testing
        result_remote = REMOTE.bioactivities.all_bioactivities(_top_n=3)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.all_bioactivities()

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])


class TestsFromKinaseIds:
    """
    Test all class methods with kinase IDs as input.
    """

    @pytest.mark.parametrize("kinase_klifs_ids", [472, [472, 509], [472, 509, 10000]])
    def test_by_kinase_klifs_id(self, kinase_klifs_ids):
        """
        Test all class methods with kinase IDs as input.
        """

        # Kinases
        result_remote = REMOTE.kinases.by_kinase_klifs_id(kinase_klifs_ids)
        result_local = LOCAL.kinases.by_kinase_klifs_id(kinase_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["kinases"])
        check_dataframe(result_local, COLUMN_NAMES["kinases"])

        # Ligands
        result_remote = REMOTE.ligands.by_kinase_klifs_id(kinase_klifs_ids)
        result_local = LOCAL.ligands.by_kinase_klifs_id(kinase_klifs_ids)
        check_dataframe(result_remote, COLUMN_NAMES["ligands"] + ["kinase.klifs_id (query)"])
        check_dataframe(result_local, COLUMN_NAMES["ligands"] + ["kinase.klifs_id (query)"])

        # Structures
        result_remote = REMOTE.structures.by_kinase_klifs_id(kinase_klifs_ids)
        result_local = LOCAL.structures.by_kinase_klifs_id(kinase_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

        # Bioactivities
        result_remote = REMOTE.bioactivities.by_kinase_klifs_id(kinase_klifs_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.by_kinase_klifs_id(kinase_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])

        # Interactions
        result_remote = REMOTE.interactions.by_kinase_klifs_id(kinase_klifs_ids)
        result_local = LOCAL.interactions.by_kinase_klifs_id(kinase_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])
        check_dataframe(result_local, COLUMN_NAMES["interactions"] + ["kinase.klifs_id (query)"])

    @pytest.mark.parametrize("kinase_klifs_ids", [10000, "XXX"])
    def test_by_kinase_klifs_id_raise(self, kinase_klifs_ids):
        """
        Test all class methods with kinase IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.by_kinase_klifs_id(kinase_klifs_ids)
            REMOTE.ligands.by_kinase_klifs_id(kinase_klifs_ids)
            REMOTE.structures.by_kinase_klifs_id(kinase_klifs_ids)
            REMOTE.bioactivities.by_kinase_klifs_id(kinase_klifs_ids)
            REMOTE.interactions.by_kinase_klifs_id(kinase_klifs_ids)

        with pytest.raises(ValueError):
            LOCAL.kinases.by_kinase_klifs_id(kinase_klifs_ids)
            LOCAL.ligands.by_kinase_klifs_id(kinase_klifs_ids)
            LOCAL.structures.by_kinase_klifs_id(kinase_klifs_ids)
            LOCAL.interactions.by_kinase_klifs_id(kinase_klifs_ids)


class TestsFromLigandIds:
    """
    Test all class methods with ligand IDs as input.
    """

    @pytest.mark.parametrize("ligand_klifs_ids", [10, [10, 20], [10, 20, 10000]])
    def test_by_ligand_klifs_id(self, ligand_klifs_ids):
        """
        Test all class methods with ligand IDs as input.
        """

        # Ligands
        result_remote = REMOTE.ligands.by_ligand_klifs_id(ligand_klifs_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.ligands.by_ligand_klifs_id(ligand_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["ligands"])

        # Structures
        result_remote = REMOTE.structures.by_ligand_klifs_id(ligand_klifs_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.structures.by_ligand_klifs_id(ligand_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])

        # Bioactivities
        result_remote = REMOTE.bioactivities.by_ligand_klifs_id(ligand_klifs_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.by_ligand_klifs_id(ligand_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])

        # Interactions
        result_remote = REMOTE.interactions.by_ligand_klifs_id(ligand_klifs_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.interactions.by_ligand_klifs_id(ligand_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])

    @pytest.mark.parametrize("ligand_klifs_ids", [10000, "XXX"])
    def test_by_ligand_klifs_id_raise(self, ligand_klifs_ids):
        """
        Test all class methods with ligand IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.bioactivities.by_ligand_klifs_id(ligand_klifs_ids)
            REMOTE.interactions.by_ligand_klifs_id(ligand_klifs_ids)

        with pytest.raises(ValueError):
            REMOTE.ligands.by_ligand_klifs_id(ligand_klifs_ids)
            REMOTE.structures.by_ligand_klifs_id(ligand_klifs_ids)


class TestsFromStructureIds:
    """
    Test class methods with structure IDs as input.
    """

    @pytest.mark.parametrize("structure_klifs_ids", [12347, [12347, 100000]])
    def test_by_structure_klifs_id(self, structure_klifs_ids):
        """
        Test class methods with structure IDs as input.
        """

        # Structures
        result_remote = REMOTE.structures.by_structure_klifs_id(structure_klifs_ids)
        result_local = LOCAL.structures.by_structure_klifs_id(structure_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

        # Interactions
        result_remote = REMOTE.interactions.by_structure_klifs_id(structure_klifs_ids)
        result_local = LOCAL.interactions.by_structure_klifs_id(structure_klifs_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])
        check_dataframe(result_local, COLUMN_NAMES["interactions"])

        # Pockets (takes only one structure ID as input!)
        if isinstance(structure_klifs_ids, int):
            structure_klifs_id = structure_klifs_ids

            result_remote = REMOTE.pockets.by_structure_klifs_id(structure_klifs_ids)
            result_local = LOCAL.pockets.by_structure_klifs_id(structure_klifs_ids)

            check_dataframe(result_remote, COLUMN_NAMES["pockets"])
            check_dataframe(result_local, COLUMN_NAMES["pockets"])

            assert all(result_local == result_remote)

    @pytest.mark.parametrize("structure_klifs_ids", [100000, "XXX"])
    def test_by_structure_klifs_id_raise(self, structure_klifs_ids):
        """
        Test class methods with structure IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.structures.by_structure_klifs_id(structure_klifs_ids)
            REMOTE.interactions.by_structure_klifs_id(structure_klifs_ids)
            if isinstance(structure_klifs_ids, int):
                structure_klifs_id = structure_klifs_ids
                REMOTE.pockets.by_structure_klifs_id(structure_klifs_id)

        with pytest.raises(ValueError):
            LOCAL.structures.by_structure_klifs_id(structure_klifs_ids)
            LOCAL.interactions.by_structure_klifs_id(structure_klifs_ids)
            if isinstance(structure_klifs_ids, int):
                structure_klifs_id = structure_klifs_ids


class TestsFromKinaseNames:
    """
    Test class methods with kinase names as input.
    """

    @pytest.mark.parametrize(
        "kinase_names, species",
        [("BMX", None), (["BMX", "BRAF"], None)],
    )
    def test_by_kinase_name(self, kinase_names, species):
        """
        Test class methods with kinase names as input.
        """

        # Kinases
        result_remote = REMOTE.kinases.by_kinase_name(kinase_names, species)
        result_local = LOCAL.kinases.by_kinase_name(kinase_names, species)

        check_dataframe(result_remote, COLUMN_NAMES["kinases"])
        check_dataframe(result_local, COLUMN_NAMES["kinases"])

        # Ligands
        result_remote = REMOTE.ligands.by_kinase_name(kinase_names)
        result_local = LOCAL.ligands.by_kinase_name(kinase_names)

        check_dataframe(
            result_remote,
            COLUMN_NAMES["ligands"]
            + [
                "kinase.klifs_id (query)",
                "kinase.klifs_name (query)",
                "kinase.gene_name (query)",
                "species.klifs (query)",
            ],
        )
        check_dataframe(
            result_local,
            COLUMN_NAMES["ligands"]
            + [
                "kinase.klifs_name (query)",
                "kinase.gene_name (query)",
                "species.klifs (query)",
            ],
        )

        # Structures
        result_remote = REMOTE.structures.by_kinase_name(kinase_names)
        result_local = LOCAL.structures.by_kinase_name(kinase_names)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("kinase_names, species", [("XXX", None), ("XXX", "XXX")])
    def test_by_kinase_name_raise(self, kinase_names, species):
        """
        Test class methods with kinase names as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.by_kinase_name(kinase_names, species)
            REMOTE.ligands.by_kinase_name(kinase_names)

        with pytest.raises(ValueError):
            REMOTE.structures.by_kinase_name(kinase_names)
            LOCAL.kinases.by_kinase_name(kinase_names, species)
            LOCAL.ligands.by_kinase_name(kinase_names)
            LOCAL.structures.by_kinase_name(kinase_names)


class TestsFromLigandPdbs:
    """
    Test class methods with Ligand Expo IDs as input.
    """

    @pytest.mark.parametrize("ligand_expo_ids", ["PRC", ["PRC", "1N1"], ["PRC", "1N1", "XXX"]])
    def test_by_ligand_expo_id(self, ligand_expo_ids):
        """
        Test class methods with Ligand Expo IDs as input.
        """

        # Ligands
        result_remote = REMOTE.ligands.by_ligand_expo_id(ligand_expo_ids)
        result_local = LOCAL.ligands.by_ligand_expo_id(ligand_expo_ids)

        check_dataframe(result_remote, COLUMN_NAMES["ligands"])
        check_dataframe(result_local, COLUMN_NAMES["ligands"])

        # Structure
        result_remote = REMOTE.structures.by_ligand_expo_id(ligand_expo_ids)
        result_local = LOCAL.structures.by_ligand_expo_id(ligand_expo_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("ligand_expo_ids", [1, "XXX"])
    def test_by_ligand_expo_id_raise(self, ligand_expo_ids):
        """
        Test class methods with Ligand Expo IDs as input: Error raised if input invalid?
        """

        with pytest.raises(ValueError):
            REMOTE.ligands.by_ligand_expo_id(ligand_expo_ids)
            REMOTE.structures.by_ligand_expo_id(ligand_expo_ids)
            LOCAL.ligands.by_ligand_expo_id(ligand_expo_ids)
            LOCAL.structures.by_ligand_expo_id(ligand_expo_ids)


class TestsFromStructurePdbs:
    """
    Test class methods with structure PDB IDs as input.
    """

    @pytest.mark.parametrize("structure_pdb_ids", ["3sxr", ["3sxr", "1fpu", "xxxx"]])
    def test_by_structure_pdb_id(self, structure_pdb_ids):
        """
        Test class methods with structure PDB IDs as input.
        """

        # Structure
        result_remote = REMOTE.structures.by_structure_pdb_id(structure_pdb_ids)
        result_local = LOCAL.structures.by_structure_pdb_id(structure_pdb_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("structure_pdb_ids", [1, "xxxx"])
    def test_by_structure_pdb_id_raise(self, structure_pdb_ids):
        """
        Test class methods with structure PDB IDs as input: Error raised if input invalid?
        """

        # Structure
        # Note: Integer input raises jsonschema.exceptions.ValidationError,
        # here tested using Exception
        with pytest.raises((SwaggerMappingError, Exception)):
            REMOTE.structures.by_structure_pdb_id(structure_pdb_ids)
        with pytest.raises(ValueError):
            LOCAL.structures.by_structure_pdb_id(structure_pdb_ids)


class TestsCoordinates:
    """
    Test remote Coordinates class.
    """

    @pytest.mark.parametrize(
        "structure_klifs_id, entity, extension, n_atoms, centroid",
        [
            (12347, "complex", "mol2", 3604, [-3.996449, 17.509910, 31.077763]),
            (12347, "complex", "pdb", 1819, [-3.903167, 17.447048, 31.263985]),
            (12347, "protein", "mol2", 3552, [-4.061604, 17.472406, 30.976721]),
            (12347, "pocket", "mol2", 1156, [0.308657, 21.880768, 35.903843]),
            (12347, "ligand", "mol2", 49, [2.291216, 20.590290, 39.074585]),
        ],
    )
    def test_remote(self, structure_klifs_id, entity, extension, n_atoms, centroid):
        """
        Test remote retrieval of coordinates data from structure ID.
        """

        # Load coordinates as DataFrame
        dataframe = REMOTE.coordinates.to_dataframe(structure_klifs_id, entity, extension)
        self._test_to_dataframe(dataframe, n_atoms, centroid)

        # Load coordinates as RDKit molecule
        if entity == "ligand" and extension == "mol2":
            rdkit_molecule = REMOTE.coordinates.to_rdkit(structure_klifs_id, entity, extension)
            self._test_to_rdkit(rdkit_molecule)

        # Save coordinates to file (to temporary directory)
        with enter_temp_directory():
            if extension == "mol2":
                filepath = REMOTE.coordinates.to_mol2(structure_klifs_id, ".", entity)
                assert filepath.exists()
            if extension == "pdb":
                filepath = REMOTE.coordinates.to_pdb(structure_klifs_id, ".", entity)
                assert filepath.exists()

    @pytest.mark.parametrize(
        "structure_klifs_id, entity, extension",
        [(12347, "XXX", "mol2"), (12347, "complex", "XXX")],
    )
    def test_remote_raise(self, structure_klifs_id, entity, extension):
        """
        Test remote retrieval of coordinates data from structure ID:
        Error raised if input invalid?
        """

        with pytest.raises(ValueError):
            REMOTE.coordinates.to_dataframe(structure_klifs_id, entity, extension)
            REMOTE.coordinates.to_rdkit(structure_klifs_id, entity, extension)

    @pytest.mark.parametrize(
        "structure_klifs_id, entity, extension, n_atoms, centroid",
        [
            (12347, "complex", "mol2", 3604, [-3.996449, 17.509910, 31.077763]),
            (12347, "complex", "pdb", 1819, [-3.903167, 17.447048, 31.263985]),
            (12347, "protein", "mol2", 3552, [-4.061604, 17.472406, 30.976721]),
            (12347, "pocket", "mol2", 1156, [0.308657, 21.880768, 35.903843]),
            (12347, "pocket", "pdb", 1156, [0.308650, 21.880758, 35.903843]),
            (12347, "ligand", "mol2", 49, [2.291216, 20.590290, 39.074585]),
            (12347, "ligand", "pdb", 31, [2.427290, 20.193062, 39.483578]),
            (12347, "water", "mol2", 3, [-29.550966, 11.602567, 20.098101]),
        ],
    )
    def test_local(self, structure_klifs_id, entity, extension, n_atoms, centroid):
        """
        Test local retrieval of coordinates data from structure ID.
        Coordinates can also be loaded locally from file, which is not tested individually, since
        loading from structure ID means (1) convert structure ID > filepath and (2) load
        coordinates from file.
        """

        # Load coordinates as DataFrame
        dataframe = LOCAL.coordinates.to_dataframe(structure_klifs_id, entity, extension)
        self._test_to_dataframe(dataframe, n_atoms, centroid)

        # Load coordinates as RDKit molecule
        if entity == "ligand":
            rdkit_molecule = LOCAL.coordinates.to_rdkit(structure_klifs_id, entity, extension)
            self._test_to_rdkit(rdkit_molecule)

    @pytest.mark.parametrize(
        "structure_klifs_id, entity, extension",
        [(12347, "XXX", "mol2"), (12347, "complex", "XXX")],
    )
    def test_local_raise(self, structure_klifs_id, entity, extension):
        """
        Test local retrieval of coordinates data from structure ID or file:
        Error raised if input invalid?
        """

        with pytest.raises(FileNotFoundError):
            LOCAL.coordinates.to_dataframe(structure_klifs_id, entity, extension)
            LOCAL.coordinates.to_rdkit(structure_klifs_id, entity, extension)

    @staticmethod
    def _test_to_dataframe(dataframe, n_atoms, centroid):
        """
        Test coordinates DataFrame that was loaded locally or remotely.
        """

        assert isinstance(dataframe, pd.DataFrame)
        assert dataframe.columns.to_list() == COLUMN_NAMES["coordinates"]
        assert dataframe.shape[0] == n_atoms
        assert centroid[0] == pytest.approx(dataframe["atom.x"].mean(), abs=1.0e-6)
        assert centroid[1] == pytest.approx(dataframe["atom.y"].mean(), abs=1.0e-6)
        assert centroid[2] == pytest.approx(dataframe["atom.z"].mean(), abs=1.0e-6)

    @staticmethod
    def _test_to_rdkit(rdkit_molecule):
        """
        Test RDKit molecule that was loaded locally or remotely.
        """

        assert isinstance(rdkit_molecule, Chem.rdchem.Mol)
