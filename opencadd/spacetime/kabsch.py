"""
Kabsch algorithm for alignment of two point clouds via rigid transformations (i.e. rotation and translation).
"""
from typing import Tuple
import jax
import jax.numpy as jnp


__author__ = "Armin Ariamajd"


@jax.jit
def kabsch_unweighted(p0: jax.Array, p1: jax.Array) -> Tuple[jax.Array, jax.Array, jax.Array]:
    """
    Align two point clouds by a rigid transformation (i.e. only rotation and translation),
    using singular value decomposition (SVD) according to the Kabsch algorithm.

    Parameters
    ----------
    p0 : jax.Array, shape: (n_samples, n_features)
        The reference point cloud.
    p1 : jax.Array, shape: (n_samples, n_features)
        The point cloud to align.

    Returns
    -------
    r, t, rmsd : 3-tuple of jax.Array
        r : jax.Array, shape: (n_features, n_features)
            Rotation matrix to align p1 onto p0.
        t : jax.Array, shape: (n_features, )
            Translation vector to align p1 onto p0.
        rmsd : jax.Array, shape: ()
            RMSD value between p0 and p1 after alignment.

        The coordinates of the aligned point cloud can be obtained as:
        `p1_aligned = p1 @ r + t`
    """
    c_p0 = jnp.mean(p0, axis=0)  # geometric center of p0
    c_p1 = jnp.mean(p1, axis=0)  # geometric center of p1
    p0c = p0 - c_p0  # p0 centered
    p1c = p1 - c_p1  # p1 centered
    # Calculate SVD of the cross-covariance matrix
    h = jnp.einsum('ji,jk->ik', p0c, p1c)  # equivalent to `p0c.T @ p1c`
    u, s, vt = jnp.linalg.svd(h)
    # Calculate rotation matrix
    rot = jnp.einsum('ij,jk->ki', u, vt)  # equivalent to `(u @ vt).T`
    det_sign = jnp.sign(jnp.linalg.det(rot))  # Sign of the determinant of the rotation matrix
    # Note: since u @ vt is a rotation matrix, it's determinant is either 1, or -1 in case of
    #  an improper rotation. In the later case, the improper rotation is corrected to a proper rotation
    #  by manipulating the s and u values:
    s = s.at[-1].multiply(det_sign)
    u = u.at[:, -1].multiply(det_sign)
    rot = jnp.einsum('ij,jk->ki', u, vt)  # final rotation matrix
    # Notice that the above code could have be written as a conditional as well:
    #     rot = jnp.einsum('ij,jk->ki', u, vt)
    #     if jnp.linalg.det(rot) < 0:
    #         s[-1] *= -1
    #         u[:, -1] *= -1
    #         rot = np.einsum('ij,jk->ki', u, vt)
    #  However, then the function could not have been jitted.
    trans = c_p0 - jnp.dot(c_p1, rot)  # translation vector
    e0 = jnp.sum(p0c ** 2 + p1c ** 2)  # initial residual
    rmsd = jnp.sqrt(jnp.maximum((e0 - (2 * jnp.sum(s))) / p0.shape[0], 0))  # rmsd after alignment
    return rot, trans, rmsd


@jax.jit
def kabsch_unweighted_transform(
        p0: jax.Array, p1: jax.Array
) -> Tuple[jax.Array, jax.Array, jax.Array, jax.Array]:
    rot, trans, rmsd = kabsch_unweighted(p0, p1)
    return rot, trans, rmsd, p1 @ rot + trans


@jax.jit
def kabsch_weighted(p0, p1, w):
    denom = jnp.sum(w)
    if w.shape == p0.shape:
        # Full weights were provided, thus sum(w) is the full sum;
        #  it only needs to be divided by n_features.
        denom /= p0.shape[1]
        denom2 = jnp.sum(w, axis=0)
    elif w.shape[0] != p0.shape[0]:
        # w has shape (1, n_feature), thus sum(w) was only for one sample;
        #  it must be multiplied by the n_samples, and then divided by n_features.
        denom *= p0.shape[0] / p0.shape[1]
        denom2 = w[0] * p0.shape[0]
    else:
        denom2 = denom

    c_p0 = jnp.sum(p0 * w, axis=0) / denom2
    c_p1 = jnp.sum(p1 * w, axis=0) / denom2
    p0c = p0 - c_p0
    p1c = p1 - c_p1
    h = jnp.einsum('ji,jk->ik', p0c * w, p1c)
    u, s, vt = jnp.linalg.svd(h)
    det_sign = jnp.sign(jnp.linalg.det(jnp.einsum('ij,jk->ki', u, vt)))
    s = s.at[-1].multiply(det_sign)
    u = u.at[:, -1].multiply(det_sign)
    rot = jnp.einsum('ij,jk->ki', u, vt)
    trans = c_p0 - jnp.dot(c_p1, rot)
    e0 = jnp.sum(w * (p0c ** 2 + p1c ** 2))
    rmsd = jnp.sqrt(jnp.maximum((e0 - (2 * jnp.sum(s))) / denom, 0))
    return rot, trans, rmsd
