from __future__ import annotations
from typing import Sequence, Union, Optional, Any, Tuple, NamedTuple, Literal, List, NoReturn
from dataclasses import dataclass

import numpy as np
import scipy as sp
import jax.numpy as jnp
import numpy.typing as npt

import opencadd as oc
from opencadd._typing import ArrayLike


class MinimumBoundingBox(NamedTuple):
    """


    Coordinates (x, y, z) of the origin point of a right rectangular prism defining
    the boundaries of the set of coordinates.
    The origin point is the vertex with smallest (x, y, z) values.
    """
    lower_bounds: jnp.ndarray
    upper_bounds: jnp.ndarray


    @property
    def lengths(self) -> jnp.ndarray:
        return jnp.abs(self.upper_bounds - self.lower_bounds)

    @property
    def volume(self) -> jnp.ndarray:
        return jnp.prod(self.lengths, axis=-1)


class TensorDataSet:
    """
    A set of features observed at different instances for each individual in a sample/population.
    The observed features can be a tensor of any rank, e.g. a scalar or a vector.
    For example, instances may be consecutive timepoints, the individuals in the sample may be a
    set of particles, and the observed feature may be the position of particles in 3-dimensional
    space, in which case the dataset represents a trajectory.
    """

    __slots__ = (
        "_data",
        "_first_axes",
        "_title_instance",
        "_title_sample",
        "_title_observation",
        "_labels_instance",
        "_labels_sample",
        "_labels_observation",
        "_kdtrees",
    )

    def __init__(
            self,
            data: npt.ArrayLike,
            first_axes: Tuple[int, int] = (1, 2),
            data_type: npt.DTypeLike = np.single,
            title_instance: Optional[Any] = None,
            labels_instance: Optional[npt.ArrayLike] = None,
            title_sample: Optional[Any] = None,
            labels_sample: Optional[npt.ArrayLike] = None,
            title_observation: Optional[Any] = None,
            labels_observation: Optional[npt.ArrayLike] = None
    ):
        """
        Parameters
        ----------
        data : ArrayLike, ndim >= 3, dtype: numeric
            An n-dimensional array (n >= 3) representing the dataset, where the first i dimensions
            represents the different instances, the next j dimensions represents the individuals
            in the sample/population, and the remaining dimensions represent the observation.
            The observation must at least be one-dimensional, i.e., even if the observation is a
            single scalar value, it must still be represented as a one-dimensional vector of
            length 1. For example, if the dataset represents the observation of scalar value for a
            sample size k, at t different timepoints, then the shape of the array must be (t, k, 1).
        data_type : DTypeLike (numeric), optional, default: numpy.single (float32)
            Type of the data.
        title_instance : Any, optional, default: None
            Type of instances in the dataset, i.e. the title of the first dimension (e.g. 'time').
            This is for the user's own future reference, and does not have any effect other than in
            the representation of the object (i.e. its `repr` and `str` dunder methods).
        labels_instance
        title_sample
        labels_sample
        title_observation
        labels_observation
        """
        self._data: jnp.ndarray
        self._title_instance: Any
        self._title_sample: Any
        self._title_observation: Any
        self._labels_instance: np.ndarray
        self._labels_sample: np.ndarray
        self._labels_observation: np.ndarray
        self._kdtrees: np.ndarray

        # Check for errors in `data`:
        data_tensor = jnp.asarray(data, dtype=data_type)
        if data_tensor.ndim < 3:
            raise ValueError(
                "Parameter `data` expects at least a 3-dimensional array. "
                f"Input argument had: {data_tensor.ndim}"
            )
        # Check for errors in first axes:
        first_axes_types = tuple(type(first_axis) for first_axis in first_axes)
        if first_axes_types != (int, int):
            raise ValueError(
                "Parameter `first_axes` expects a tuple of two integers. "
                f"Intput argument was: {first_axes} : {first_axes_types}."
            )
        if not 1 <= first_axes[0] < first_axes[1] < data_tensor.ndim:
            raise ValueError(
                "Parameter `firs_axes` expects a tuple of two strictly increasing integers "
                "in the half-open interval $[1, n_d)$, where $n_d$ is the number of dimensions "
                f"in argument `data`. Input argument was {first_axes}, "
                f"and `data` had {data_tensor.ndim} dimensions."
            )

        self._data = data_tensor
        self._first_axes = np.asarray(first_axes).astype(np.ubyte)

        num_instances = np.prod(self._data.shape[:self._first_axes[0]])
        num_individuals = np.prod(self._data.shape[self._first_axes[0]:self._first_axes[1]])
        num_features = np.prod(self._data.shape[self._first_axes[1]:])

        view_total = self._data.reshape(-1, num_features)
        view_per_instance = self._data.reshape(num_instances, -1, num_features)
        view_per_individual = self._data.reshape(num_individuals, )


        first_axes_array = np.asarray(first_axes).astype(np.uby)


        self._title_instance = title_instance
        self._title_sample = title_sample
        self._title_observation = title_observation
        if labels_instance is None:
            self._labels_instance = np.arange(self._data.shape[0])
        else:
            labels_instance_array = np.asarray(labels_instance)
            if labels_instance_array.shape[0] != self._data.shape[0]:
                raise ValueError("Shape of `labels_instance` along first dimension must match "
                                 "the number of instances in the dataset.")
            else:
                self._labels_instance = labels_instance_array



        self._kdtrees = [sp.spatial.KDTree(data=data_tensor)]
        return

    def __repr__(self):
        return (
            "TensorDataSet("
            f"instances({self.title_instance}): {self.count_instance}, "
            f"samples({self.title_sample}): {self.count_sample}, "
            f"observations({self.title_observation}): {self.shape_observation}"
            ")"
        )

    def __str__(self):
        pass

    def __getitem__(self, item):
        return self._data.__getitem__(item)

    @property
    def title_instance(self) -> Any:
        """Title of the instance type of data."""
        return self._title_instance

    @property
    def title_sample(self) -> Any:
        """Title of the sample type of data."""
        return self._title_sample

    @property
    def title_observation(self) -> Any:
        """Title of the observation type of data."""
        return self._title_observation

    @property
    def labels_instance(self) -> np.ndarray:
        return self._labels_instance

    @property
    def labels_sample(self) -> np.ndarray:
        return self._labels_sample

    @property
    def labels_observation(self) -> np.ndarray:
        return self._labels_observation

    @property
    def count_instance(self) -> int:
        """Number of instances in the data."""
        return self._data.shape[0]

    @property
    def count_sample(self) -> int:
        """Sample size of the data."""
        return self._data.shape[1]

    @property
    def shape_observation(self) -> Tuple[int]:
        """Shape of the observation tensor."""
        return self._data.shape[2:]

    @property
    def count_points(self) -> int:
        """
        Number of coordinates (points).
        """
        return self._data.shape[0]

    def axis_aligned_minimum_bounding_box(
            self,
            space: Literal["instance", "individual", "feature"] = "instance",
    ) -> jnp.ndarray:
        """
        Lengths (along x, y, and z axes, respectively) of a right rectangular prism defining
        the boundaries of the set of coordinates.
        """
        axes = {
            "instance": tuple(range(*self._first_axes)),
            "individual": tuple(range(self._first_axes[0])),
            "feature": tuple(range(self._first_axes[1]))
        }
        mins = jnp.min(self._data, axis=axes[space])
        maxes = jnp.max(self._data, axis=axes[space])
        return self._kdtree.maxes - self._kdtree.mins


