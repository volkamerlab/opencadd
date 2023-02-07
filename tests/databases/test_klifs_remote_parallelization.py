"""
Tests for opencadd.databases.klifs.remote parallelization!
"""

from itertools import repeat
from multiprocessing import Pool

import pytest

from opencadd.data.klifs import setup_remote

# Set local session
REMOTE = setup_remote()


@pytest.mark.parametrize(
    "klifs_session, structure_klifs_ids, n_cores",
    [(REMOTE, [111, 113], 2)],
)
def test_remote_parallelization(klifs_session, structure_klifs_ids, n_cores):

    pool = Pool(processes=n_cores)

    pool.starmap(_parallelize_klifs_session, zip(structure_klifs_ids, repeat(klifs_session)))
    pool.close()
    pool.join()


def _parallelize_klifs_session(structure_klifs_id, klifs_session):
    """
    Dummy function that simulates a function that takes the KLIFS session as input;
    shall be used in our test where we want to pass a single klifs_session to multiple functions
    that are called in a parallelized process.
    """
    return klifs_session.structures.by_structure_klifs_id(structure_klifs_id)
