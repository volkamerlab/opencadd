from __future__ import annotations
from typing import Sequence, Union, Optional, Any, Tuple, Self, NamedTuple, Literal, List, NoReturn

import numpy as np
import scipy as sp
import jax.numpy as jnp
import numpy.typing as npt

import opencadd as oc
from opencadd.typing import ArrayLike


class DynamicPointCloud:
    """

    """

    __slots__ = (
        "_data",
        "_data_view_2d",
        "_kdtree_combined",
        "_kdtrees_per_instance",
    )

    def __init__(self, data: npt.ArrayLike):
        # Declare attributes
        self._data: jnp.ndarray
        self._data_view_2d: jnp.ndarray
        self._kdtree_combined: sp.spatial.KDTree
        self._kdtrees_per_instance: List[sp.spatial.KDTree]
        # Check for errors in `data`:
        data_tensor = jnp.asarray(data)
        if data_tensor.ndim != 3:
            raise ValueError(
                "Parameter `data` expects a 3-dimensional array. "
                f"Input argument had: {data_tensor.ndim}"
            )

        # Assign attributes
        self._data = data_tensor
        self._data_view_2d = self._data.reshape(-1, self._data.shape[-1])
        self._kdtree_combined = None
        self._kdtrees_per_instance = None
        return

    @property
    def count_instances(self) -> int:
        return self._data.shape[0]

    @property
    def count_points_per_instance(self) -> int:
        return self._data.shape[1]

    @property
    def count_points_total(self) -> int:
        return np.prod(self._data.shape[:2])

    @property
    def dimension_points(self) -> int:
        return self._data.shape[2]

    @property
    def _kdtree_total(self) -> sp.spatial.KDTree:
        if self._kdtree_combined is None:
            self._kdtree_combined = sp.spatial.KDTree(self._data_view_2d)
        return self._kdtree_combined

    @property
    def _kdtrees_list(self) -> List[sp.spatial.KDTree]:
        if self._kdtrees_per_instance is None:
            self._kdtrees_per_instance = [sp.spatial.KDTree(data=points) for points in self._data]
        return self._kdtrees_per_instance

    def axis_aligned_minimum_bounding_box(
            self,
            per_instance: bool = True,
    ) -> oc.spacetime.volume.RectangularCuboid:
        """
        Axis-aligned minimum bounding box of the point cloud, either per instance or as a whole.

        Parameters
        ----------
        per_instance : bool, optional, default: True
            Whether to calculate one bounding box per instance, or one for all instances
            superposed.

        Returns
        -------
        opencadd.spacetime.volume.RectangularCuboid
            A `RectangularCuboid` volume object representing the axis-aligned minimum bounding
            box(es).

        References
        ----------
        * https://en.wikipedia.org/wiki/Minimum_bounding_box
        """
        if per_instance:
            mins = jnp.min(self._data, axis=1)
            maxes = jnp.max(self._data, axis=1)
        else:
            mins = jnp.expand_dims(jnp.min(self._data, axis=(0, 1)), axis=0),
            maxes = jnp.expand_dims(jnp.max(self._data, axis=(0, 1)), axis=0)
        return oc.spacetime.volume.RectangularCuboid(lower_bounds=mins, upper_bounds=maxes)

    def toxelate(
            self,
            resolution: float,
            radius_points: Union[float, npt.ArrayLike],
            padding: float = 0,
    ) -> oc.spacetime.volume.Toxel:
        # Get the bounding box of all instances superposed.
        total_bounding_box = self.axis_aligned_minimum_bounding_box(per_instance=False)
        # Create a grid the size of the total bounding box, with given resolution
        grid = oc.spacetime.grid.from_bounds_and_spacing(
            lower_bounds=total_bounding_box.lower_bounds - padding,
            upper_bounds=total_bounding_box.upper_bounds + padding,
            spacings=resolution
        )
        # If `radius_points` is a scalar (i.e. int or float), it means all points have the same
        # radius, and thus we only need to query for the first nearest neighbor of each point:
        if (
                np.issubdtype(type(radius_points), np.floating)
                or np.issubdtype(type(radius_points), np.integer)
        ):
            # Radius must be positive:
            if radius_points <= 0:
                raise ValueError("Parameter `radius_points` expects positive values.")

            dists, indices = self.nearest_neighbors(
                points=grid.coordinates_2d,
                num_neaerst_neighbors=1,
                per_instance=True,
                distance_upper_bound=radius_points
            )
            # Each toxel on the grid is occupied when the nearest point in self is within
            # `radius_points`, i.e. when it is not `np.inf`, since `self.nearest_neighbors` returns
            # `np.inf` for points where the nearest distance is larger than the upper bound.
            toxel_tensor = dists != np.inf  # True when toxel is occupied
            # Create `ToxelField`:
            toxel_field = oc.spacetime.field.ToxelField(
                field_tensor=toxel_tensor,
                grid=grid,
            )
            # Create Toxel volume from field and return
            return oc.spacetime.volume.Toxel(field=toxel_field)
        # If `radius_points` is an array of values, then we cannot rely only on the distances to
        # first nearest neighbors, since it is possible that the first k nearest neighbors have
        # small radii and do not overlap with the toxel, while the (k+1)-th neighbor has a large
        # enough radius to overlap.
        else:
            radii_array = np.asarray(radius_points)
            max_radius = radii_array.max()
            min_radius = radii_array.min()
            toxel_tensor = np.zeros(
                shape=(self.count_instances, grid.shape, 1),
                dtype=np.bool_
            )
            ind_self, ind_gird, dists = self.distance_matrix_sparse(
                points=grid, max_distance=max_radius
            )
            filter_definitely_occupied = dists <= min_radius
            grid_inds_occupied = np.unravel_index(
                ind_gird[1][filter_definitely_occupied],
                shape=grid.shape
            )
            toxel_tensor[(ind_self[0][filter_definitely_occupied], *grid_inds_occupied)] = True
            filter_maybe_occupied = jnp.logical_not(filter_definitely_occupied)
            dists_from_surface = dists[filter_maybe_occupied] - radii_array[ind_self[1][
                filter_maybe_occupied]]
            occupied = dists_from_surface <= 0
            grid_inds_occupied2 = np.unravel_index(
                ind_gird[1][filter_maybe_occupied][occupied],
                shape=grid.shape
            )
            toxel_tensor[
                (ind_self[0][filter_maybe_occupied][occupied], *grid_inds_occupied2)
            ] = True
            toxel_field = oc.spacetime.field.ToxelField(
                field_tensor=toxel_tensor,
                grid=grid,
            )
            # Create Toxel volume from field and return
            return oc.spacetime.volume.Toxel(field=toxel_field)

    def distance_matrix_sparse(
            self,
            points: DynamicPointCloud,
            max_distance: float,
            p_norm: float = 2,
            output_type: Literal[
                'nd_unraveled', 'dok_matrix', 'coo_matrix', 'dict', 'ndarray'
            ] = 'nd_unraveled'
    ):
        dist_matrix = self._kdtree_total.sparse_distance_matrix(
            other=points._kdtree_total,
            max_distance=max_distance,
            p=p_norm,
            output_type='ndarray' if output_type == 'nd_unraveled' else output_type
        )
        if output_type != 'nd_unraveled':
            return dist_matrix
        else:
            indices_self = np.unravel_index(
                dist_matrix["i"],
                shape=(self.count_instances, self.count_points_per_instance)
            )
            indices_other = np.unravel_index(
                dist_matrix["j"],
                shape=(points.count_instances, points.count_points_per_instance)
            )
            return indices_self, indices_other, dist_matrix["v"]






    def find_point_pairs_within_radius(
            self,

    ):
        pass

    def count_neighbors_within_radius(self, points: Self):
        pass


    def minimize_bounding_box(self) -> Self:
        """
        Rotate the observations to minimize the bounding box volume.

        Returns
        -------

        References
        ----------
        Google: how to rotate a shape so that enclosing box volume is minimized
        * https://en.wikipedia.org/wiki/Minimum_bounding_rectangle
        * https://gis.stackexchange.com/questions/22895/finding-minimum-area-rectangle-for-given-points
        * https://math.stackexchange.com/questions/2342844/how-to-find-the-rotation-which-minimizes-the-volume-of-the-bounding-box
        * https://perso.uclouvain.be/chia-tche.chang/resources/CGM11_paper.pdf
        """
        raise NotImplementedError

    def nearest_neighbors(
            self,
            points: npt.ArrayLike,
            num_neaerst_neighbors: Union[int, Sequence[int]] = 1,
            per_instance: bool = True,
            error_tolerance: float = 0,
            p_norm: float = 2,
            distance_upper_bound: float = np.inf,
            distance_dtype: npt.DTypeLike = np.single,
    ):
        """
        For each point in `coordinates`, find the distances to, and indices of, a given number of
        nearest points in self.

        Parameters
        ----------
        points : numpy.ndarray, shape: (d_1, d_2, ..., d_{m-1}, self.dimension_points)
            Coordinates of n points (n = d_1 * d_2 * ... * d_{m-1}),
            for which the nearest points in self must be found.
        num_neaerst_neighbors : int | Sequence[int], optional, default: 1
            Either the number of nearest neighbors (as an integer), or a sequence of
            the k-th nearest neighbors to find.
        error_tolerance : float, optional, default: 0
            Tolerance for error in finding the nearest atoms. The k-th nearest atom will be
            within (1 + eps) times the distance to the real k-th nearest atom.
        p_norm : float, range: [1, inf), optional, default: 2
            The Minkowski p-norm to use, e.g.:
                * 1: Manhattan distance, i.e. sum-of-absolute-values distance.
                * 2: Euclidean distance.
                * inf: Maximum-coordinate-difference distance.
        distance_upper_bound : float, range: [0, inf), optional, default: inf
            Prune the search tree to return only neighbors within this range.

        Returns
        -------
        distances, indices : Tuple[ndarray, ndarray], shape: (d_1, d_2, ..., d_{m-1}, k)
            Distances to, and indices of the k nearest points in self, to each point in
            `coordinates`. Both returned arrays match the dimensions of `coordinates` along all
            axes but last; the last axis has k elements, corresponding to the k nearest neighbors.
        """
        if np.issubdtype(type(num_neaerst_neighbors), np.integer):
            if num_neaerst_neighbors < 1:
                raise ValueError(
                    "Parameter `num_neaerst_neighbors` expects positive (nonzero) integers. "
                    f"Input argument was: {num_neaerst_neighbors}"
                )
            k_neighbors = tuple(range(1, num_neaerst_neighbors + 1))
        else:
            num_neaerst_neighbors_array = np.asarray(num_neaerst_neighbors)
            if not np.issubdtype(num_neaerst_neighbors_array.dtype, np.integer):
                raise ValueError("Parameter `num_neaerst_neighbors` expects integers.")
            if num_neaerst_neighbors_array.ndim != 1:
                raise ValueError("Parameter `num_neaerst_neighbors` expects 1-dimensional arrays.")
            k_neighbors = tuple(num_neaerst_neighbors_array)

        points_array = jnp.asarray(points)

        if per_instance:
            shape_distances = (self.count_instances, *points_array.shape[:-1], len(k_neighbors))
            shape_indices = (self.count_instances, *points_array.shape[:-1], len(k_neighbors), 2)
            distances = np.empty(shape=shape_distances, dtype=distance_dtype)
            indices = np.empty(
                shape=shape_indices,
                dtype=oc.typing.smallest_integer_dtype_for_range(
                    min_val=0,
                    max_val=self.count_points_per_instance
                )
            )
            for idx_instance, kdtree in enumerate(self._kdtrees_per_instance):
                indices[idx_instance, ..., 0] = idx_instance
                distances[idx_instance], indices[idx_instance, ..., 1] = kdtree.query(
                    x=points_array,
                    k=k_neighbors,
                    eps=error_tolerance,
                    p=p_norm,
                    distance_upper_bound=distance_upper_bound,
                    workers=-1
                )
            return distances, indices
        else:
            raise NotImplementedError


    def count_neighbors(self):
        pass


    def cluster__common_nearest_neighbor(
            self,
            radius_neighborhood: float,
            min_samples: int,
            metric: Literal = 'euclidean',
            metric_params=None,
            leaf_size=30,
            p_norm: int = 2,
    ):
        """

        Parameters
        ----------
        radius_neighborhood
        min_samples
        metric
        metric_params
        leaf_size
        p_norm

        Returns
        -------

        References
        ----------
        Documentation on scikit-learn-extra, with examples:
        * https://scikit-learn-extra.readthedocs.io/en/stable/modules/cluster.html#common-nearest-neighbors-clustering
        * https://scikit-learn-extra.readthedocs.io/en/stable/auto_examples/plot_commonnn.html#sphx-glr-auto-examples-plot-commonnn-py
        * https://scikit-learn-extra.readthedocs.io/en/stable/auto_examples/cluster/plot_commonnn_data_sets.html
        Source code of the scikit-learn-extra implementation:
        * https://github.com/scikit-learn-contrib/scikit-learn-extra/tree/main/sklearn_extra/cluster
        GitHub and documentation of the independent package:
        * https://github.com/bkellerlab/CommonNNClustering
        * https://bkellerlab.github.io/CommonNNClustering
        Publications:
        * https://doi.org/10.1063/1.3301140
        * https://doi.org/10.3390/a11020019
        * https://doi.org/10.1063/1.4965440
        """
