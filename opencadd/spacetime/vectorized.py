from typing import Tuple, Optional
import operator
import numpy as np
import jax
import jax.numpy as jnp


def filter_array(
        array: np.ndarray,
        minmax_vals: Tuple[float, float],
        reduction_op: str,
        reduction_axis: int = -1,
        invert: bool = False
) -> np.ndarray:
    """
    Filter points by their distances to their nearest atoms.

    Parameters
    ----------
    array
    dist_range : tuple[float, float]
        Lower and upper bounds of distance to the nearest atom.
    invert : bool, optional, default: False
        If False (default), the points whose nearest atoms lie within the distance range are
        returned. If True, the points whose nearest atoms lie outside the range are returned.

    Returns
    -------
    boolean_mask : ndarray
    """
    red_funcs_pre = {"max": np.max, "min": np.min, "avg": np.mean, "sum": np.sum}
    red_funcs_post = {"all": np.all, "any": np.any, "one": np.logical_xor.reduce}

    if reduction_op in red_funcs_pre:
        array_values = red_funcs_pre[reduction_op](array, axis=reduction_axis)
        array_mask = filter_array_by_range(
            array=array_values,
            minmax_vals=minmax_vals,
            invert=invert
        )
        return array_mask, (array_values.min(), array_values.max())
    elif reduction_op in red_funcs_post:
        mask_array = filter_array_by_range(
            array=array,
            minmax_vals=minmax_vals,
            invert=invert
        )
        reduced_mask_array = red_funcs_post[reduction_op](mask_array, axis=reduction_axis)
        return reduced_mask_array, (array.min(), array().max())
    else:
        raise ValueError("Reduction operation not recognized.")


def filter_array_by_range(
        array: np.ndarray,
        minmax_vals: Tuple[float, float],
        invert: bool = False
):
    op_lower, op_upper = (operator.le, operator.ge) if invert else (operator.ge, operator.le)
    op_combine = np.logical_or if invert else np.logical_and
    mask = op_combine(
        op_lower(array, minmax_vals[0]),
        op_upper(array, minmax_vals[1]),
    )
    return mask


def _rmsd_unweighted(p1: jax.Array, p2: jax.Array) -> jax.Array:
    """
    Root-mean-square deviation (RMSD) between two sets of points in n-dimensional euclidean space.

    Parameters
    ----------
    p1 : jax.Array, shape: (n_points, n_dim)
        Coordinates of the first set of points.
    p2 : jax.Array, shape: (n_points, n_dim)
        Coordinates of the second set of points.

    Returns
    -------
    jax.Array
        RMSD between the two sets as a 0-dimensional array (i.e. a scalar)
    """
    return jnp.sqrt(jnp.sum((p1 - p2) ** 2) / p1.shape[0])


def _rmsd_weighted(p1: jax.Array, p2: jax.Array, weights: jax.Array) -> jax.Array:
    n_dim = p1.shape[1]
    sd = (p1 - p2) ** 2
    w_sd = sd * weights
    w_msd = w_sd.sum() / (weights.sum() / n_dim)
    w_rmsd = jnp.sqrt(w_msd)
    return w_rmsd
