import numpy as np

from opencadd.spacetime import spatial


def test_grid_distance():
    grid = np.zeros((3, 3, 3))
    grid[1, 1, 1] = 1
    dists = spatial.grid_distance(
        grid,
        directions=spatial.GRID_DIRS,
        start_indices=np.argwhere(grid),
        target_value=0,
        max_dist_in_dir=np.ones(spatial.GRID_DIRS.shape[0])
    )
    assert np.array_equal(dists, np.ones(shape=(1, 26)))