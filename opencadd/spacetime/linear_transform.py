
import jax
import jax.numpy as jnp


@jax.jit
def affine_map(p: jax.Array, linear_map: jax.Array, translation: jax.Array) -> jax.Array:
    """
    Apply an affine map (linear map and translation) to an n-dimensional array.

    Parameters
    ----------
    p : jax.Array
        Array to be mapped.
    linear_map : jax.Array
        Array representation of the linear map.
    trans : jax.Array
        Translation vector.

    Returns
    -------
    p_transformed : jax.Array
    """
    return p @ linear_map + translation
