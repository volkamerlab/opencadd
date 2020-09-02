"""
Tests for opencadd.structure.subpockets.core
"""

from pathlib import Path

from bravado.client import SwaggerClient
import pandas as pd
import pytest

from opencadd.io.dataframes import Mol2ToDataFrame
from opencadd.structure.subpockets.core import Base, AnchorResidue, Pocket


def load_dataframe_protein(klifs_structure_id=None):
    """
    Load protein DataFrame by KLIFS structure ID. If ID is None, load simple DataFrame.
    """

    if klifs_structure_id:

        KLIFS_API_DEFINITIONS = "http://klifs.vu-compmedchem.nl/swagger/swagger.json"
        KLIFS_CLIENT = SwaggerClient.from_url(
            KLIFS_API_DEFINITIONS, config={"validate_responses": False}
        )

        result = (
            KLIFS_CLIENT.Structures.get_structure_get_protein(structure_ID=klifs_structure_id)
            .response()
            .result
        )

        with open("protein.mol2", "w") as f:
            f.write(result)

        parser = Mol2ToDataFrame()
        dataframe = parser.from_file("protein.mol2")

        # Delete file
        Path("protein.mol2").unlink()

    else:

        dataframe = pd.DataFrame(
            {
                "residue.pdb_id": ["1", "2", "3", "7", "8", "9", "11"],
                "atom.name": ["CA", "CA", "CA", "CA", "CA", "CA", "CA"],
                "atom.x": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
                "atom.y": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
                "atom.z": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
            }
        )

    return dataframe


class TestsBase:
    """
    Test Base class methods.
    """

    @pytest.mark.parametrize("color_name, color_rgb", [("red", (1, 0, 0))])
    def test_format_color(self, color_name, color_rgb):
        base = Base()
        name, rgb = base._format_color(color_name)
        assert name == color_name
        assert rgb == color_rgb

    @pytest.mark.parametrize("color_name", ["xxx"])
    def test_format_color_raise(self, color_name):

        with pytest.raises(ValueError):
            base = Base()
            base._format_color(color_name)

    @pytest.mark.parametrize(
        "residue_pdb_ids, residue_labels, residue_pdb_ids_formatted, residue_labels_formatted",
        [
            ([1, 2, 3], [1, 2, 3], ["1", "2", "3"], ["1", "2", "3"]),
            ([1, 2], None, ["1", "2"], [None, None]),
        ],
    )
    def test_format_residue_pdb_ids_and_labels(
        self, residue_pdb_ids, residue_labels, residue_pdb_ids_formatted, residue_labels_formatted
    ):

        base = Base()
        residue_pdb_ids, residue_labels = base._format_residue_pdb_ids_and_labels(
            residue_pdb_ids, residue_labels
        )
        assert residue_pdb_ids == residue_pdb_ids_formatted
        assert residue_labels == residue_labels_formatted


class TestsAnchorResidue:
    """
    Test AnchorResidue class methods.
    """

    @pytest.mark.parametrize(
        "residue_pdb_id, residue_pdb_id_alternative, residue_center",
        [
            (1, None, [1, 1, 1]),
            (10, ["9", "11"], [35, 35, 35]),
            (4, ["3"], [3, 3, 3]),
            (6, ["7"], [4, 4, 4]),
            (5, None, None),
        ],
    )
    def test_from_dataframe(self, residue_pdb_id, residue_pdb_id_alternative, residue_center):
        """
        Test anchor residue cases.
        """

        dataframe = load_dataframe_protein()

        residue = AnchorResidue()
        residue.from_dataframe(dataframe, residue_pdb_id)

        assert residue.pdb_id_alternative == residue_pdb_id_alternative
        if residue_center:
            assert pytest.approx(residue.center, residue_center, abs=1.0e-6)


class TestsPocket:
    """
    Test Pocket class methods.
    """

    @pytest.mark.parametrize(
        "name, residue_pdb_ids, residue_labels, centroid",
        [("example kinase", [127, 128, 129], [1, 2, 3], [-1.178433, 23.859733, 45.091933])],
    )
    def test_init(self, name, residue_pdb_ids, residue_labels, centroid):

        dataframe = load_dataframe_protein(3834)
        pocket = Pocket(dataframe, name, residue_pdb_ids, residue_labels)

        assert pocket.name == name
        assert pocket.residues.equals(
            pd.DataFrame(
                {
                    "residue.pdb_id": [str(i) for i in residue_pdb_ids],
                    "residue.label": [str(i) for i in residue_labels],
                }
            )
        )
        assert pytest.approx(pocket.centroid, centroid, abs=1.0e-6)

    @pytest.mark.parametrize(
        "name, color, residue_pdb_ids, residue_labels",
        [("hinge", "magenta", [127, 128, 129], [46, 47, 48])],
    )
    def test_add_regions(self, name, color, residue_pdb_ids, residue_labels):

        dataframe = load_dataframe_protein(3834)
        pocket = Pocket(dataframe, "", residue_pdb_ids, residue_labels)

        pocket.add_region(name, color, residue_pdb_ids, residue_labels)

        n_region_residues = len(residue_pdb_ids)
        region = pd.DataFrame(
            {
                "region.name": [name] * n_region_residues,
                "region.color_name": [color] * n_region_residues,
                "residue.pdb_id": [str(i) for i in residue_pdb_ids],
                "residue.label": [str(i) for i in residue_labels],
            }
        )
        assert pocket.regions.drop("region.color_rgb", axis=1).equals(region)

    @pytest.mark.parametrize(
        "name, color, residue_pdb_ids, residue_labels, center",
        [
            (
                "hinge_region",
                "magenta",
                [73, 128, 193],
                [16, 47, 80],
                [1.957633, 21.923767, 41.69003333],
            )
        ],
    )
    def test_add_subpocket(self, name, color, residue_pdb_ids, residue_labels, center):

        dataframe = load_dataframe_protein(3834)
        pocket = Pocket(dataframe, "", residue_pdb_ids, residue_labels)

        pocket.add_subpocket(name, color, residue_pdb_ids, residue_labels)
        subpocket = pocket._subpockets[0]

        assert subpocket.name == name
        assert subpocket.color_name == color
        assert pytest.approx(subpocket.center, center, abs=1.0e-6)
        assert pocket.subpockets.columns.to_list() == [
            "subpocket.name",
            "subpocket.color_name",
            "subpocket.color_rgb",
            "subpocket.center",
        ]