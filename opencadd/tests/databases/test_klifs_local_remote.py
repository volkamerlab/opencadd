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

# Set local and remote session
REMOTE = setup_remote()
LOCAL = setup_local(Path(__name__).parent / "opencadd" / "tests" / "data" / "klifs")


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
        "group, local_families", [(None, ["Tec", "RAF", "Abl"]), ("TK", ["Tec", "Abl"])],
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

        assert result_local["kinase.name"].to_list() == local_kinases[0]
        assert result_local["species.klifs"].to_list() == local_kinases[1]
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    @pytest.mark.parametrize(
        "group, family, species", [("XXX", None, None), (None, "XXX", None), ("XXX", None, "XXX")],
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

        assert result_local["ligand.pdb"].to_list() == ["1N1", "QH1", "PRC"]
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

        assert result_local["structure.id"].to_list() == [3482, 12347, 5728, 5705]
        # Do not test remote,
        # since too many and may vary if structures are added to KLIFS.

    def test_interaction_types(self):
        """
        Test request result for all interaction types.
        """

        result_remote = REMOTE.interactions.interaction_types
        check_dataframe(result_remote, COLUMN_NAMES["interaction_types"])

        with pytest.raises(NotImplementedError):
            LOCAL.interactions.interaction_types

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
        result_remote = REMOTE.bioactivities.all_bioactivities(n=3)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.all_bioactivities()

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])


class TestsFromKinaseIds:
    """
    Test all class methods with kinase IDs as input.
    """

    @pytest.mark.parametrize("kinase_ids", [472, [472, 509], [472, 509, 10000]])
    def test_from_kinase_ids(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input.
        """

        # Kinases
        result_remote = REMOTE.kinases.from_kinase_ids(kinase_ids)
        result_local = LOCAL.kinases.from_kinase_ids(kinase_ids)

        check_dataframe(result_remote, COLUMN_NAMES["kinases"])
        check_dataframe(result_local, COLUMN_NAMES["kinases"])

        # Ligands
        result_remote = REMOTE.ligands.from_kinase_ids(kinase_ids)
        result_local = LOCAL.ligands.from_kinase_ids(kinase_ids)
        check_dataframe(result_remote, COLUMN_NAMES["ligands"] + ["kinase.id (query)"])
        check_dataframe(result_local, COLUMN_NAMES["ligands"] + ["kinase.id (query)"])

        # Structures
        result_remote = REMOTE.structures.from_kinase_ids(kinase_ids)
        result_local = LOCAL.structures.from_kinase_ids(kinase_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

        # Bioactivities
        result_remote = REMOTE.bioactivities.from_kinase_ids(kinase_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.from_kinase_ids(kinase_ids)

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])

        # Interactions
        result_remote = REMOTE.interactions.from_kinase_ids(kinase_ids)
        result_local = LOCAL.interactions.from_kinase_ids(kinase_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])
        check_dataframe(result_local, COLUMN_NAMES["interactions"] + ["kinase.id (query)"])

    @pytest.mark.parametrize("kinase_ids", [10000, "XXX"])
    def test_from_kinase_ids_raise(self, kinase_ids):
        """
        Test all class methods with kinase IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.from_kinase_ids(kinase_ids)
            REMOTE.ligands.from_kinase_ids(kinase_ids)
            REMOTE.structures.from_kinase_ids(kinase_ids)
            REMOTE.bioactivities.from_kinase_ids(kinase_ids)
            REMOTE.interactions.from_kinase_ids(kinase_ids)

        with pytest.raises(ValueError):
            LOCAL.kinases.from_kinase_ids(kinase_ids)
            LOCAL.ligands.from_kinase_ids(kinase_ids)
            LOCAL.structures.from_kinase_ids(kinase_ids)
            LOCAL.interactions.from_kinase_ids(kinase_ids)


class TestsFromLigandIds:
    """
    Test all class methods with ligand IDs as input.
    """

    @pytest.mark.parametrize("ligand_ids", [10, [10, 20], [10, 20, 10000]])
    def test_from_ligand_ids(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input.
        """

        # Ligands
        result_remote = REMOTE.ligands.from_ligand_ids(ligand_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.ligands.from_ligand_ids(ligand_ids)

        check_dataframe(result_remote, COLUMN_NAMES["ligands"])

        # Structures
        result_remote = REMOTE.structures.from_ligand_ids(ligand_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.structures.from_ligand_ids(ligand_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])

        # Bioactivities
        result_remote = REMOTE.bioactivities.from_ligand_ids(ligand_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.bioactivities.from_ligand_ids(ligand_ids)

        check_dataframe(result_remote, COLUMN_NAMES["bioactivities"])

        # Interactions
        result_remote = REMOTE.interactions.from_ligand_ids(ligand_ids)
        with pytest.raises(NotImplementedError):
            LOCAL.interactions.from_ligand_ids(ligand_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])

    @pytest.mark.parametrize("ligand_ids", [10000, "XXX"])
    def test_from_ligand_ids_raise(self, ligand_ids):
        """
        Test all class methods with ligand IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.bioactivities.from_ligand_ids(ligand_ids)
            REMOTE.interactions.from_ligand_ids(ligand_ids)

        with pytest.raises(ValueError):
            REMOTE.ligands.from_ligand_ids(ligand_ids)
            REMOTE.structures.from_ligand_ids(ligand_ids)


class TestsFromStructureIds:
    """
    Test class methods with structure IDs as input.
    """

    @pytest.mark.parametrize("structure_ids", [12347, [12347, 100000]])
    def test_from_structure_ids(self, structure_ids):
        """
        Test class methods with structure IDs as input.
        """

        # Structures
        result_remote = REMOTE.structures.from_structure_ids(structure_ids)
        result_local = LOCAL.structures.from_structure_ids(structure_ids)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

        # Interactions
        result_remote = REMOTE.interactions.from_structure_ids(structure_ids)
        result_local = LOCAL.interactions.from_structure_ids(structure_ids)

        check_dataframe(result_remote, COLUMN_NAMES["interactions"])
        check_dataframe(result_local, COLUMN_NAMES["interactions"])

        # Pockets (takes only one structure ID as input!)
        if isinstance(structure_ids, int):
            structure_id = structure_ids

            result_remote = REMOTE.pockets.from_structure_id(structure_ids)
            result_local = LOCAL.pockets.from_structure_id(structure_ids)

            check_dataframe(result_remote, COLUMN_NAMES["pockets"])
            check_dataframe(result_local, COLUMN_NAMES["pockets"])

            assert all(result_local == result_remote)

    @pytest.mark.parametrize("structure_ids", [100000, "XXX"])
    def test_from_structure_ids_raise(self, structure_ids):
        """
        Test class methods with structure IDs as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.structures.from_structure_ids(structure_ids)
            REMOTE.interactions.from_structure_ids(structure_ids)
            if isinstance(structure_ids, int):
                structure_id = structure_ids
                REMOTE.pockets.from_structure_id(structure_id)

        with pytest.raises(ValueError):
            LOCAL.structures.from_structure_ids(structure_ids)
            LOCAL.interactions.from_structure_ids(structure_ids)
            if isinstance(structure_ids, int):
                structure_id = structure_ids


class TestsFromKinaseNames:
    """
    Test class methods with kinase names as input.
    """

    @pytest.mark.parametrize(
        "kinase_names, species", [("BMX", None), (["BMX", "BRAF"], None)],
    )
    def test_from_kinase_names(self, kinase_names, species):
        """
        Test class methods with kinase names as input.
        """

        # Kinases
        result_remote = REMOTE.kinases.from_kinase_names(kinase_names, species)
        result_local = LOCAL.kinases.from_kinase_names(kinase_names, species)

        check_dataframe(result_remote, COLUMN_NAMES["kinases"])
        check_dataframe(result_local, COLUMN_NAMES["kinases"])

        # Ligands
        result_remote = REMOTE.ligands.from_kinase_names(kinase_names)
        result_local = LOCAL.ligands.from_kinase_names(kinase_names)

        check_dataframe(
            result_remote,
            COLUMN_NAMES["ligands"]
            + ["kinase.id (query)", "kinase.name (query)", "species.klifs (query)",],
        )
        check_dataframe(
            result_local,
            COLUMN_NAMES["ligands"] + ["kinase.name (query)", "species.klifs (query)",],
        )

        # Structures
        result_remote = REMOTE.structures.from_kinase_names(kinase_names)
        result_local = LOCAL.structures.from_kinase_names(kinase_names)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("kinase_names, species", [("XXX", None), (1, None), ("XXX", "XXX")])
    def test_from_kinase_names_raise(self, kinase_names, species):
        """
        Test class methods with kinase names as input: Error raised if input invalid?
        """

        with pytest.raises(SwaggerMappingError):
            REMOTE.kinases.from_kinase_names(kinase_names, species)
            REMOTE.ligands.from_kinase_names(kinase_names)

        with pytest.raises(ValueError):
            REMOTE.structures.from_kinase_names(kinase_names)
            LOCAL.kinases.from_kinase_names(kinase_names, species)
            LOCAL.ligands.from_kinase_names(kinase_names)
            LOCAL.structures.from_kinase_names(kinase_names)


class TestsFromLigandPdbs:
    """
    Test class methods with ligand PDB IDs as input.
    """

    @pytest.mark.parametrize("ligand_pdbs", ["PRC", ["PRC", "1N1"], ["PRC", "1N1", "XXX"]])
    def test_from_ligand_pdbs(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input.
        """

        # Ligands
        result_remote = REMOTE.ligands.from_ligand_pdbs(ligand_pdbs)
        result_local = LOCAL.ligands.from_ligand_pdbs(ligand_pdbs)

        check_dataframe(result_remote, COLUMN_NAMES["ligands"])
        check_dataframe(result_local, COLUMN_NAMES["ligands"])

        # Structure
        result_remote = REMOTE.structures.from_ligand_pdbs(ligand_pdbs)
        result_local = LOCAL.structures.from_ligand_pdbs(ligand_pdbs)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("ligand_pdbs", [1, "XXX"])
    def test_from_ligand_pdbs_raise(self, ligand_pdbs):
        """
        Test class methods with ligand PDB IDs as input: Error raised if input invalid?
        """

        with pytest.raises(ValueError):
            REMOTE.ligands.from_ligand_pdbs(ligand_pdbs)
            REMOTE.structures.from_ligand_pdbs(ligand_pdbs)
            LOCAL.ligands.from_ligand_pdbs(ligand_pdbs)
            LOCAL.structures.from_ligand_pdbs(ligand_pdbs)


class TestsFromStructurePdbs:
    """
    Test class methods with structure PDB IDs as input.
    """

    @pytest.mark.parametrize("structure_pdbs", ["3sxr", ["3sxr", "1fpu", "xxxx"]])
    def test_from_structure_pdbs(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input.
        """

        # Structure
        result_remote = REMOTE.structures.from_structure_pdbs(structure_pdbs)
        result_local = LOCAL.structures.from_structure_pdbs(structure_pdbs)

        check_dataframe(result_remote, COLUMN_NAMES["structures"])
        check_dataframe(result_local, COLUMN_NAMES["structures"])

    @pytest.mark.parametrize("structure_pdbs", [1, "xxxx"])
    def test_from_structure_pdbs_raise(self, structure_pdbs):
        """
        Test class methods with structure PDB IDs as input: Error raised if input invalid?
        """

        # Structure
        # Note: Integer input raises jsonschema.exceptions.ValidationError,
        # here tested using Exception
        with pytest.raises((SwaggerMappingError, Exception)):
            REMOTE.structures.from_structure_pdbs(structure_pdbs)
        with pytest.raises(ValueError):
            LOCAL.structures.from_structure_pdbs(structure_pdbs)


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

        result_remote = REMOTE.coordinates.from_structure_id(
            structure_id, entity, input_format, output_format
        )
        result_local = LOCAL.coordinates.from_structure_id(
            structure_id, entity, input_format, output_format
        )

        self._test_from_structure_id(result_remote, output_format, n_atoms, centroid)
        self._test_from_structure_id(result_local, output_format, n_atoms, centroid)

    @staticmethod
    def _test_from_structure_id(coordinates, output_format, n_atoms, centroid):

        if output_format == "biopandas":
            assert isinstance(coordinates, pd.DataFrame)
            assert coordinates.columns.to_list() == COLUMN_NAMES["coordinates"]
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

        with pytest.raises(ValueError):
            REMOTE.coordinates.from_structure_id(structure_id, entity, input_format, output_format)
            LOCAL.coordinates.from_structure_id(structure_id, entity, input_format, output_format)
