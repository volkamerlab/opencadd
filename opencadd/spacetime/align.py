import numpy as np


def scipy(p0, p1, weights=None):
    if weights is None:
        weights = np.ones(len(p1))
    else:
        weights = np.asarray(weights)

    h = np.einsum('ji,jk->ik', weights[:, None] * p0, p1)

    u, s, vt = np.linalg.svd(h)

    if np.linalg.det(u @ vt) < 0:
        s[-1] = -s[-1]
        u[:, -1] = -u[:, -1]

    rot = np.dot(u, vt)

    rssd = np.sqrt(max(
        np.sum(weights * np.sum(p1 ** 2 + p0 ** 2, axis=1)) - 2 * np.sum(s),
        0))

    return rot, rssd


def optAlign(p0, p1):
    """
    optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
    Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
    Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA

    Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
    PyMol based upon your selections.

    By default, this program will optimally align the ALPHA CARBONS of the selections provided.
    To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.

    @param p0: First PyMol selection with N-atoms
    @param p1: Second PyMol selection with N-atoms
    """

    # must alway center the two proteins to avoid
    # affine transformations.  Center the two proteins to their selections.


    # Initial residual, see Kabsch.
    E0 = np.sum(np.sum(p0_centered**2, axis=0), axis=0) + np.sum(np.sum(p1_centered**2, axis=0), axis=0)

    u, s, vt = np.linalg.svd(np.dot(p1_centered.T, p0_centered))

    reflect = float(str(float(np.linalg.det(u) * np.linalg.det(vt))))

    if reflect == -1.0:
        s[-1] = -s[-1]
        u[:, -1] = -u[:, -1]

    RMSD = E0 - (2.0 * sum(s))
    RMSD = np.sqrt(abs(RMSD / L))

    # U is simply V*Wt
    rot = np.dot(u, vt)

    # rotate and translate the molecule
    p1 = np.dot((mol2 - center_p1), rot)
    # center the molecule
    p0 = mol1 - center_p0
