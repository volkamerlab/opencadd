"""
Base class for all the superposition engines.
"""
import logging

_logger = logging.getLogger(__name__)


class BaseAligner:
    """ """

    def calculate(self, structures, *args, **kwargs):
        """
        Compute superposition for given structures

        Parameters
        ----------
        structures : list of opencadd.core.Structure

        Returns
        -------
        dict
            Results with following keys

            - ``superposed``: ``Structure`` objects, aligned in place
            - ``scores``: Dictionary of relevant metrics. Usually includes ``rmsd``.
            - ``metadata``: Contextual information provided by the method

        """
        self._safety_checks()
        assert (
            len(structures) == 2,
            "This method can only be used for two structures at the same time, for now",
        )
        return self._calculate(structures, *args, **kwargs)

    def _calculate(self, structures, *args, **kwargs):
        raise NotImplementedError("Reimplement in your subclass")

    def _safety_checks(self):
        raise NotImplementedError("Reimplement in your subclass")
