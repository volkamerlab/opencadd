from typing import Literal

import jax
import jax.numpy as jnp

import opencadd as oc
from opencadd import _typing


@jax.jit
def _rmsd_1_1_s(
    p0: jax.Array,
    p1: jax.Array
) -> jax.Array:
    """
    RMSD between two instances.

    Parameters
    ----------
    p0, p1 : jax.Array
        An instance as a 2-Dimensional array of shape (n_samples, n_features).

    Returns
    -------
    rmsd : jax.Array
        RMSD value between the two sets, as a 0-dim. array (i.e. a float).
    """
    return jnp.sqrt(jnp.sum((p0 - p1) ** 2) / p0.shape[0])


_rmsd_1_1_v = jax.vmap(_rmsd_1_1_s, in_axes=(0, 0))
"""
Pairwise RMSD between two arrays of instances.

Parameters
----------
p0, p1 : jax.Array
    3-Dimensional array of shape (n_instances, n_samples, n_features).

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of size n_instances, where rmsd[i] is the RMSD value between p0[i] and p1[i]. 
"""


_rmsd_1_2_s = jax.vmap(_rmsd_1_1_s, in_axes=(None, 0))
"""
Pairwise RMSD between a reference instance, and an array of instances.

Parameters
----------
p0 : jax.Array, shape: (n_samples, n_features)
    The reference instance, as a 2D array of shape (n_samples, n_features).
p1 : jax.Array, shape: (n_instances, n_samples, n_features)
    Array of instances, as a 3D array of shape (n_instances, n_samples, n_features).

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of shape (n_instances,), where rmsd[i] is the RMSD value between p0 and p1[i]. 
"""


_rmsd_1_2_v = jax.vmap(_rmsd_1_2_s, in_axes=(0, 0))
"""
Pairwise RMSD between a reference instance, and an array of instances, over a batch dimension.

Parameters
----------
p0 : jax.Array, shape: (n_batch_instances, n_samples, n_features)
    The reference instances, as a 3D array of shape (n_batch_instances, n_samples, n_features).
p1 : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features)
    Array of instances, as a 4D array of shape (n_batch_instances, n_instances, n_samples, n_features).

Returns
-------
rmsd : jax.Array, shape: (n_batch_instances, n_instances)
    2D array of shape (n_batch_instances, n_instances), where rmsd[b, i] 
    is the RMSD value between p0[b] and p1[b, i]. 
"""


_rmsd_2_2_s = jax.vmap(_rmsd_1_2_s, in_axes=(0, None))
"""
Pairwise RMSD between all permutations of two array of instances.

Parameters
----------
p0 : jax.Array, shape: (n_instances_0, n_samples, n_features)
    The first instances, as a 3D array of shape (n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_instances_1, n_samples, n_features)
    The second instances, as a 3D array of shape (n_instances_1, n_samples, n_features).

Returns
-------
rmsd : jax.Array, shape: (n_instances_0, n_instances_1)
    2D array of shape (n_instances_0, n_instances_1), where rmsd[i, j] 
    is the RMSD value between p0[i] and p1[j]. 
"""


_rmsd_2_2_v = jax.vmap(_rmsd_2_2_s, in_axes=(0, 0))
"""
Pairwise RMSD between all permutations of two array of instances, over a batch dimension.

Parameters
----------
p0 : jax.Array, shape: (n_batches, n_instances_0, n_samples, n_features)
    The first instances, as a 4D array of shape (n_batches, n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_batches, n_instances_1, n_samples, n_features)
    The second instances, as a 4D array of shape (n_batches, n_instances_1, n_samples, n_features).

Returns
-------
rmsd : jax.Array, shape: (n_batches, n_instances_0, n_instances_1)
    3D array of shape (n_batches, n_instances_0, n_instances_1), where rmsd[b, i, j] 
    is the RMSD value between p0[b, i] and p1[b, j]. 
"""


@jax.jit
def _wrmsd_1_1_s(
    p0: jax.Array,
    p1: jax.Array,
    w: jax.Array,
) -> jax.Array:
    """
    Weighted RMSD between two instances.

    Parameters
    ----------
    p0, p1 : jax.Array, shape: (n_samples, n_features)
        Instances.
    w : jax.Array, shape: (n_samples, n_features) or (n_samples, 1) or (1, n_features).
        Weights.

    Returns
    -------
    rmsd : jax.Array
        Weighted RMSD value between the two sets, as a 0-dim. array (i.e. a float).
    """
    denom = jnp.sum(w)
    if w.shape == p0.shape:
        # Full weights were provided, thus sum(w) is the full sum;
        #  it only needs to be divided by n_features.
        denom /= p0.shape[1]
    elif w.shape[0] != p0.shape[0]:
        # w has shape (1, n_feature), thus sum(w) was only for one sample;
        #  it must be multiplied by the n_samples, and then divided by n_features.
        denom *= p0.shape[0] / p0.shape[1]
    # The last possible case is that w has shape (n_samples, 1); in this case
    #  sum(w) must be multiplied by n_features to get the full sum,
    #  but then it must also be divided by n_features as well;
    #  thus, no action is required.
    return jnp.sqrt(jnp.sum(w * ((p0 - p1) ** 2)) / denom)


_wrmsd_1_1_v_s = jax.vmap(_wrmsd_1_1_s, in_axes=(0, 0, None))
"""
Pairwise weighted RMSD between two arrays of instances, with constant weights.

Parameters
----------
p0, p1 : jax.Array
    3-Dimensional array of shape (n_instances, n_samples, n_features).
w : jax.Array, shape: (n_samples, n_features) or (n_samples, 1) or (1, n_features).
    Weights used for all pairs of instances.

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of size n_instances, where rmsd[i] is the weighted RMSD value between p0[i] and p1[i]. 
"""

_wrmsd_1_1_v_v = jax.vmap(_wrmsd_1_1_s, in_axes=(0, 0, 0))
"""
Pairwise weighted RMSD between two arrays of instances, with different weights per instance.

Parameters
----------
p0, p1 : jax.Array
    3-Dimensional array of shape (n_instances, n_samples, n_features).
w : jax.Array, shape: (n_instances, n_samples, n_features) or (n_instances, n_samples, 1) 
or (n_instances, 1, n_features).
    Weights for each pair of instances, i.e. w[i] is the weights for the pair p0[i] and p1[i].

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of size n_instances, where rmsd[i] is the RMSD value between p0[i] and p1[i]. 
"""

_wrmsd_1_2_s_s = jax.vmap(_wrmsd_1_1_s, in_axes=(None, 0, None))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, with constant weights.

Parameters
----------
p0 : jax.Array, shape: (n_samples, n_features)
    The reference instance, as a 2D array of shape (n_samples, n_features).
p1 : jax.Array, shape: (n_instances, n_samples, n_features)
    Array of instances, as a 3D array of shape (n_instances, n_samples, n_features).
w : jax.Array, shape: (n_samples, n_features) or (n_samples, 1) or (1, n_features).
    Weights used for all pairs of instances.

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of shape (n_instances,), where rmsd[i] is the weighted RMSD value between p0 and p1[i]. 
"""

_wrmsd_1_2_s_v = jax.vmap(_wrmsd_1_1_s, in_axes=(None, 0, 0))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, 
with different weights for each instance.

Parameters
----------
p0 : jax.Array, shape: (n_samples, n_features)
    The reference instance, as a 2D array of shape (n_samples, n_features).
p1 : jax.Array, shape: (n_instances, n_samples, n_features)
    Array of instances, as a 3D array of shape (n_instances, n_samples, n_features).
w : jax.Array, shape: (n_instances, n_samples, n_features) or (n_instances, n_samples, 1) 
or (n_instances, 1, n_features).
    Weights for each pair of instances, i.e. w[i] is the weights for the pair p0 and p1[i].

Returns
-------
rmsd : jax.Array, shape: (n_instances,)
    1D array of shape (n_instances,), where rmsd[i] is the weighted RMSD value between p0 and p1[i]. 
"""

_wrmsd_1_2_v_ss = jax.vmap(_wrmsd_1_2_s_s, in_axes=(0, 0, None))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, over a batch dimension, 
with constant weights.

Parameters
----------
p0 : jax.Array, shape: (n_batch_instances, n_samples, n_features)
    The reference instances, as a 3D array of shape (n_batch_instances, n_samples, n_features).
p1 : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features)
    Array of instances, as a 4D array of shape (n_batch_instances, n_instances, n_samples, n_features).
w : jax.Array, shape: (n_samples, n_features) or (n_samples, 1) or (1, n_features).
    Weights used for all pairs of instances.

Returns
-------
rmsd : jax.Array, shape: (n_batch_instances, n_instances)
    2D array of shape (n_batch_instances, n_instances), where rmsd[b, i] 
    is the weighted RMSD value between p0[b] and p1[b, i]. 
"""


_wrmsd_1_2_v_vs = jax.vmap(_wrmsd_1_2_s_s, in_axes=(0, 0, 0))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, over a batch dimension, 
with different weights for each batch.

Parameters
----------
p0 : jax.Array, shape: (n_batch_instances, n_samples, n_features)
    The reference instances, as a 3D array of shape (n_batch_instances, n_samples, n_features).
p1 : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features)
    Array of instances, as a 4D array of shape (n_batch_instances, n_instances, n_samples, n_features).
w : jax.Array, shape: (n_batch_instances, n_samples, n_features) 
or (n_batch_instances, n_samples, 1) or (n_batch_instances, 1, n_features).
    Weights for each pair of instances, i.e. w[i] is the weights for all pairs between p0[i] 
    and instances in p1[i].

Returns
-------
rmsd : jax.Array, shape: (n_batch_instances, n_instances)
    2D array of shape (n_batch_instances, n_instances), where rmsd[b, i] 
    is the weighted RMSD value between p0[b] and p1[b, i]. 
"""

_wrmsd_1_2_v_sv = jax.vmap(_wrmsd_1_2_s_v, in_axes=(0, 0, None))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, over a batch dimension, 
with different weights for each instance in a batch, but same for all batches..

Parameters
----------
p0 : jax.Array, shape: (n_batch_instances, n_samples, n_features)
    The reference instances, as a 3D array of shape (n_batch_instances, n_samples, n_features).
p1 : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features)
    Array of instances, as a 4D array of shape (n_batch_instances, n_instances, n_samples, n_features).
w : jax.Array, shape: (n_instances, n_samples, n_features) 
or (n_instances, n_samples, 1) or (n_instances, 1, n_features).
    Weights for each pair of instances, i.e. w[i] is the weights for pairs between p0[b] 
    and p1[b, i] in all batches b.

Returns
-------
rmsd : jax.Array, shape: (n_batch_instances, n_instances)
    2D array of shape (n_batch_instances, n_instances), where rmsd[b, i] 
    is the weighted RMSD value between p0[b] and p1[b, i]. 
"""

_wrmsd_1_2_v_vv = jax.vmap(_wrmsd_1_2_s_v, in_axes=(0, 0, 0))
"""
Pairwise weighted RMSD between a reference instance, and an array of instances, over a batch dimension, 
with different weights for each batch, and for each instance in a batch.

Parameters
----------
p0 : jax.Array, shape: (n_batch_instances, n_samples, n_features)
    The reference instances, as a 3D array of shape (n_batch_instances, n_samples, n_features).
p1 : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features)
    Array of instances, as a 4D array of shape (n_batch_instances, n_instances, n_samples, n_features).
w : jax.Array, shape: (n_batch_instances, n_instances, n_samples, n_features) 
or (n_batch_instances, n_instances, n_samples, 1) or (n_batch_instances, n_instances, 1, n_features).
    Weights for each pair of instances, i.e. w[b, i] is the weights for the pair between p0[b] 
    and p1[b, i].

Returns
-------
rmsd : jax.Array, shape: (n_batch_instances, n_instances)
    2D array of shape (n_batch_instances, n_instances), where rmsd[b, i] 
    is the weighted RMSD value between p0[b] and p1[b, i]. 
"""


_wrmsd_2_2_s_ss = jax.vmap(_wrmsd_1_2_s_s, in_axes=(0, None, None))
"""
Pairwise weighted RMSD between all permutations of two array of instances, with constant weights.

Parameters
----------
p0 : jax.Array, shape: (n_instances_0, n_samples, n_features)
    The first instances, as a 3D array of shape (n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_instances_1, n_samples, n_features)
    The second instances, as a 3D array of shape (n_instances_1, n_samples, n_features).
w : jax.Array, shape: (n_samples, n_features) or (n_samples, 1) or (1, n_features).
    Weights used for all pairs of instances.
    
Returns
-------
rmsd : jax.Array, shape: (n_instances_0, n_instances_1)
    2D array of shape (n_instances_0, n_instances_1), where rmsd[i, j] 
    is the weighted RMSD value between p0[i] and p1[j]. 
"""


_wrmsd_2_2_s_vs = jax.vmap(_wrmsd_1_2_s_s, in_axes=(0, None, 0))
"""
Pairwise weighted RMSD between all permutations of two array of instances, 
with different weights for each instance in p0.

Parameters
----------
p0 : jax.Array, shape: (n_instances_0, n_samples, n_features)
    The first instances, as a 3D array of shape (n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_instances_1, n_samples, n_features)
    The second instances, as a 3D array of shape (n_instances_1, n_samples, n_features).
w : jax.Array, shape: (n_instances_0, n_samples, n_features) 
or (n_instances_0, n_samples, 1) or (n_instances_0, 1, n_features).
    Weights for each pair of instances with p1, i.e. w[i] is the weights for all pairs between p0[i] 
    and instances in p1.

Returns
-------
rmsd : jax.Array, shape: (n_instances_0, n_instances_1)
    2D array of shape (n_instances_0, n_instances_1), where rmsd[i, j] 
    is the weighted RMSD value between p0[i] and p1[j]. 
"""


_wrmsd_2_2_s_sv = jax.vmap(_wrmsd_1_2_s_v, in_axes=(0, None, None))
"""
Pairwise weighted RMSD between all permutations of two array of instances, 
with different weights for each instance in p1.

Parameters
----------
p0 : jax.Array, shape: (n_instances_0, n_samples, n_features)
    The first instances, as a 3D array of shape (n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_instances_1, n_samples, n_features)
    The second instances, as a 3D array of shape (n_instances_1, n_samples, n_features).
w : jax.Array, shape: (n_instances_1, n_samples, n_features) 
or (n_instances_1, n_samples, 1) or (n_instances_1, 1, n_features).
    Weights for each pair of instances with p0, i.e. w[i] is the weights for all pairs 
    between instances of p0 with p1[i].

Returns
-------
rmsd : jax.Array, shape: (n_instances_0, n_instances_1)
    2D array of shape (n_instances_0, n_instances_1), where rmsd[i, j] 
    is the weighted RMSD value between p0[i] and p1[j]. 
"""

_wrmsd_2_2_s_vv = jax.vmap(_wrmsd_1_2_s_v, in_axes=(0, None, 0))
"""
Pairwise weighted RMSD between all permutations of two array of instances, 
with different weights for each pair.

Parameters
----------
p0 : jax.Array, shape: (n_instances_0, n_samples, n_features)
    The first instances, as a 3D array of shape (n_instances_0, n_samples, n_features).
p1 : jax.Array, shape: (n_instances_1, n_samples, n_features)
    The second instances, as a 3D array of shape (n_instances_1, n_samples, n_features).
w : jax.Array, shape: (n_instances_0, n_instances_1, n_samples, n_features) 
or (n_instances_0, n_instances_1, n_samples, 1) or (n_instances_0, n_instances_1, 1, n_features).
    Weights for each pair of instances, i.e. w[i, j] is the weights for pairs p0[i] and p1[j].  

Returns
-------
rmsd : jax.Array, shape: (n_instances_0, n_instances_1)
    2D array of shape (n_instances_0, n_instances_1), where rmsd[i, j] 
    is the weighted RMSD value between p0[i] and p1[j]. 
"""


def rmsd(
        p0: oc._typing.ArrayLike,
        p1: oc._typing.ArrayLike,
        weights: oc._typing.ArrayLike = None,
        stacked: bool = False,
        weights_axis: Literal[0, 1] = 0
) -> jax.Array:
    """
    Root-mean-square deviation (RMSD) between datasets.

    Both weighted and unweighted RMSD values can be calculated for different combinations of datasets,
    for example, between two datasets, between one dataset and a series of datasets, or between two
    series of datasets.

    A dataset is a 2D array of shape (n, d), where n is the number of samples (i.e. number of datapoints)
    and d is the number of features (i.e. dimension of each datapoint, which may also be 1 for scalar values).
    For example, the dataset may represent atomic coordinates of a molecule, in which case n is equal to the
    number of atoms in the molecule, and d is equal to 3, corresponding to the x-, y- and z-coordinates of
    each atom.




    The weights must be at least two-dimensional, and the last two dimensions of weights must be broadcastable
    to the data arrays, meaning they must either be (n, d), (n, 1), or (1, d).

    Parameters
    ----------
    p0
    p1
    weights
    stacked
    weights_axis

    Returns
    -------

    """

    p0 = jnp.asarray(p0)
    p1 = jnp.asarray(p1)
    p_dims = (p0.ndim, p1.ndim)
    for p_dim in p_dims:
        if not 2 <= p_dim <= 4:
            raise ValueError()
    if weights is not None:
        w = jnp.asarray(weights)
        w_shape = w.shape[-2:]
        p_shape = p0.shape[-2:]
        if w_shape != p_shape and (
            w_shape == (1, 1)
            or w_shape[0] not in (1, p_shape[0])
            or w_shape[1] not in (1, p_shape[1])
        ):
            raise ValueError()
    else:
        w = None

    if p_dims == (2, 2):
        if p0.shape != p1.shape:
            raise ValueError()
        if w is None:
            return _rmsd_1_1_s(p0, p1)
        if w.ndim != 2:
            raise ValueError()
        return _wrmsd_1_1_s(p0, p1, w)

    elif p_dims in ((2, 3), (3, 2)):
        if p_dims == (3, 2):
            p0, p1 = p1, p0
        if p0.shape != p1.shape[1:]:
            raise ValueError()
        if w is None:
            return _rmsd_1_2_s(p0, p1)
        if w.ndim == 2:
            return _wrmsd_1_2_s_s(p0, p1, w)
        elif w.ndim == 3:
            if w.shape[0] != p1.shape[0]:
                raise ValueError()
            return _wrmsd_1_2_s_v(p0, p1, w)
        raise ValueError("Dimension Mismatch")

    elif p_dims == (3, 3):
        if p0.shape[1:] != p1.shape[1:]:
            raise ValueError
        if stacked:
            if p0.shape[0] != p1.shape[0]:
                raise ValueError
            if w is None:
                return _rmsd_1_1_v(p0, p1)
            if w.ndim == 2:
                return _wrmsd_1_1_v_s(p0, p1, w)
            elif w.ndim == 3:
                if w.shape[0] != p1.shape[0]:
                    raise ValueError()
                return _wrmsd_1_1_v_v(p0, p1, w)
            raise ValueError("Dimension Mismatch")
        if w is None:
            return _rmsd_2_2_s(p0, p1)
        if w.ndim == 2:
            return _wrmsd_2_2_s_ss(p0, p1, w)
        if w.ndim == 3:
            w_shape = w.shape[0]
            fits_p0 = w_shape == p0.shape[0]
            fits_p1 = w_shape == p1.shape[0]
            if fits_p0 and not fits_p1:
                return _wrmsd_2_2_s_vs(p0, p1, w)
            elif fits_p1 and not fits_p0:
                return _wrmsd_2_2_s_sv(p0, p1, w)
            elif fits_p0 and fits_p1:
                if weights_axis == 0:
                    return _wrmsd_2_2_s_vs(p0, p1, w)
                if weights_axis == 1:
                    return _wrmsd_2_2_s_sv(p0, p1, w)
                raise ValueError
            raise ValueError
        if w.ndim == 4:
            if w.shape[:2] != (p0.shape[0], p1.shape[0]):
                raise ValueError
            return _wrmsd_2_2_s_vv(p0, p1, w)
        raise ValueError

    elif p_dims in ((3, 4), (4, 3)):
        if p_dims == (4, 3):
            p0, p1 = p1, p0
        if p0.shape != (p1.shape[0], *p1.shape[2:]):
            raise ValueError()
        if w is None:
            vals = _rmsd_1_2_v(p0, p1)
        elif w.ndim == 2:
            vals = _wrmsd_1_2_v_ss(p0, p1, w)
        elif w.ndim == 3:
            w_shape = w.shape[0]
            fits_p0 = w_shape == p0.shape[0]
            fits_p1 = w_shape == p1.shape[1]
            if fits_p0 and not fits_p1:
                vals = _wrmsd_1_2_v_vs(p0, p1, w)
            elif fits_p1 and not fits_p0:
                vals = _wrmsd_1_2_v_sv(p0, p1, w)
            elif fits_p0 and fits_p1:
                if weights_axis == 0:
                    vals = _wrmsd_1_2_v_vs(p0, p1, w)
                elif weights_axis == 1:
                    vals = _wrmsd_1_2_v_sv(p0, p1, w)
                else:
                    raise ValueError
            else:
                raise ValueError
        elif w.ndim == 4:
            if w.shape[:2] != (p0.shape[0], p1.shape[0]):
                raise ValueError
            vals = _wrmsd_1_2_v_vv(p0, p1, w)
        else:
            raise ValueError("Dimension Mismatch")
        return vals if p_dims == (3, 4) else vals.T

    elif p_dims == (4, 4):
        if (p0.shape[0], *p0.shape[2:]) != (p1.shape[0], *p1.shape[2:]):
            raise ValueError
        if w is None:
            return _rmsd_2_2_v(p0, p1)
        else:
            raise NotImplementedError

    raise ValueError





