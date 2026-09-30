"""
_bands.py
    Shared helpers for the band-comparison scripts (compare_bands_qe.py, compare_bands_w90_qe.py,
    validate_wannier_bands.py): graphene high-symmetry points, path-corner detection/labelling, and
    the two ways of putting bands on a common zero (Dirac point).

    The raw QE bands.dat reader and the Cartesian->reduced conversion live in qe_io
    (get_qe_bands, kcart_to_kred).
"""
import numpy as np

# hexagonal high-symmetry points (reduced coords) + symmetry images, to label path corners.
# graphene convention: K = (2/3, 1/3, 0), M = (1/2, 0, 0), Gamma = (0, 0, 0).
SYM_POINTS = {
    r"$\Gamma$": [(0, 0, 0)],
    "K": [(2/3, 1/3, 0), (1/3, 2/3, 0), (-1/3, 1/3, 0), (1/3, -1/3, 0), (-1/3, -2/3, 0), (2/3, -1/3, 0)],
    "M": [(1/2, 0, 0), (0, 1/2, 0), (1/2, 1/2, 0), (-1/2, 1/2, 0), (1/2, -1/2, 0), (0, -1/2, 0)],
}


def path_corners(kcart):
    """Indices of the path endpoints and direction-change vertices (high-symmetry corners)."""
    seg = np.diff(kcart, axis=0)
    n = np.linalg.norm(seg, axis=1, keepdims=True)
    t = seg / np.where(n > 0, n, 1.0)
    corners = [0]
    for i in range(1, len(t)):
        if np.dot(t[i], t[i - 1]) < 0.999:        # tangent direction changes -> a corner
            corners.append(i)
    corners.append(len(kcart) - 1)
    return corners


def label_corner(kred, tol=2e-2):
    """Label a reduced-coord k-point with its high-symmetry name (modulo a reciprocal vector), or ''."""
    for name, imgs in SYM_POINTS.items():
        for s in imgs:
            d = kred - np.array(s)
            if np.linalg.norm(d - np.round(d)) < tol:
                return name
    return ""


def nearest_index(kred, targets):
    """Index of the path point closest (modulo a reciprocal vector) to any of `targets`."""
    best_i, best_d = 0, np.inf
    for i, kr in enumerate(kred):
        for s in targets:
            d = kr - np.array(s)
            dd = np.linalg.norm(d - np.round(d))
            if dd < best_d:
                best_d, best_i = dd, i
    return best_i


def dirac_energy(eig, nocc=4):
    """
    Dirac-point energy = midpoint of the smallest top-valence/bottom-conduction gap along the path.

    Robust to the path parametrization and to a rigid energy-reference offset; needs only the number
    of occupied bands (graphene: 4, the pi/pi* cone touches between 0-indexed bands 3 and 4). Returns
    NaN if there are too few bands to bracket the cone.
    """
    eig = np.asarray(eig)
    if eig.shape[1] <= nocc:
        return np.nan
    es = np.sort(eig, axis=1)
    gap = es[:, nocc] - es[:, nocc - 1]
    ik = int(np.argmin(gap))
    return 0.5 * (es[ik, nocc] + es[ik, nocc - 1])
