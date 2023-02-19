"""
Test if opencadd is up-to-date with KLIFS database (website and download).
If errors are raised, it is time to update opencadd.
"""

from opencadd.db.klifs.remote import KLIFS_CLIENT
from opencadd.db.klifs.schema import FIELDS


class TestSyncKlifsSwaggerWithOpencadd:
    """
    Test if opencadd is up-to-date with the KLIFS OpenAPI (remote!).
    """

    def _test_klifs_model(self, data_opencadd, data_klifs):
        """
        Check if opencadd is up-to-date with KLIFS models.
        """

        # Get kinases details keys in opencadd
        keys_opencadd = set(sorted(data_opencadd.keys()))
        # Get kinases details keys in KLIFS
        result = data_klifs.response().result[0]
        keys_klifs = set(sorted(list(result)))

        assert keys_opencadd == keys_klifs

    def test_all_kinases(self):
        """
        Check if opencadd is up-to-date with XXX model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("kinases_all"), KLIFS_CLIENT.Information.get_kinase_names()
        )

    def test_kinases(self):
        """
        Check if opencadd is up-to-date with KinaseInformation model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("kinases"),
            KLIFS_CLIENT.Information.get_kinase_information(),
        )

    def test_ligands(self):
        """
        Check if opencadd is up-to-date with ligandDetails model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("ligands"),
            KLIFS_CLIENT.Ligands.get_ligands_list(kinase_ID=[1]),
        )

    def test_structures(self):
        """
        Check if opencadd is up-to-date with ligandDetails model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("structures"),
            KLIFS_CLIENT.Structures.get_structure_list(structure_ID=[1]),
        )

    def test_bioactivities(self):
        """
        Check if opencadd is up-to-date with ligandDetails model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("bioactivities"),
            KLIFS_CLIENT.Ligands.get_bioactivity_list_id(ligand_ID=2),
        )

    def test_interaction_types(self):
        """
        Check if opencadd is up-to-date with InteractionList model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("interaction_types"),
            KLIFS_CLIENT.Interactions.get_interactions_get_types(),
        )

    def test_pockets(self):
        """
        Check if opencadd is up-to-date with MatchList model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("pockets"),
            KLIFS_CLIENT.Interactions.get_interactions_match_residues(structure_ID=100),
        )

    def test_interactions(self):
        """
        Check if opencadd is up-to-date with InteractionList model.
        """

        self._test_klifs_model(
            FIELDS.remote_to_oc_names("interactions"),
            KLIFS_CLIENT.Interactions.get_interactions_get_IFP(structure_ID=[1]),
        )
