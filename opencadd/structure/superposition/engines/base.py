"""
Base class for all the superposition engines.
"""
import logging

_logger = logging.getLogger(__name__)


class BaseAligner:
    """ """

    def calculate(self, structures, selections, *args, **kwargs):
        """
        Compute superposition for given structures

        Parameters
        ----------
        structures : list of opencadd.core.Structure
        selections : list of MDAnalysis.core.groups.AtomGroup
            The selection is done by user input and the api provides the selection in form of AtomGroups

        Returns
        -------
        dict
            Results with following keys

            - ``superposed``: ``Structure`` objects, aligned in place
            - ``scores``: Dictionary of relevant metrics. Usually includes ``rmsd``.
            - ``metadata``: Contextual information provided by the method

        """
        self._safety_checks()

        message = "This method can only be used for two structures at the same time, for now"
        assert len(structures) == 2, message

        return self._calculate(structures, selections, *args, **kwargs)

    def _calculate(self, structures, selections, *args, **kwargs):
        raise NotImplementedError("Reimplement in your subclass")

    def _safety_checks(self):
        raise NotImplementedError("Reimplement in your subclass")
