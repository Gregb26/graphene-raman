"""
alignment.py
    Potential alignment between a defective and a pristine supercell (R4, 2026-09-25).

    far_atom_alignment -- Lu et al. (2019)-type alignment: average of the local KS potential (pp.x plot_num=1)
                          in a sphere around the atom farthest from the vacancy, in the defective and in the
                          pristine supercell; the difference is the rigid shift of the defective energy reference.
    rigid_shift_fit    -- rigid shift between two spectra, fitted on the states below an energy (outside the window).
    Units: potentials and energies as given (the caller converts), positions reduced, radii in the unit of A.
"""
import numpy as np


def min_image_dist(x_red, center_red, A_cols):
    """Minimum-image Cartesian distances of reduced positions x_red (n, 3) to center_red."""
    d = (np.asarray(x_red) - np.asarray(center_red) + 0.5) % 1.0 - 0.5
    return np.linalg.norm(d @ A_cols.T, axis=1)


def vacancy_site(x_red_p, x_red_d, A_cols):
    """Reduced position of the pristine atom without counterpart in the defective supercell (the vacancy)."""
    dmin = np.array([min_image_dist(x_red_d, p, A_cols).min() for p in x_red_p])
    i = int(np.argmax(dmin))
    return np.mod(x_red_p[i], 1.0), i, float(dmin[i])


def far_atom(x_red_d, s_vac, A_cols):
    """Index and distance of the atom of the defective cell farthest (minimum image) from the vacancy."""
    d = min_image_dist(x_red_d, s_vac, A_cols)
    i = int(np.argmax(d))
    return i, float(d[i])


def sphere_average(V, A_cols, center_red, radius):
    """
    Mean of V[ix, iy, iz] over the grid points within a minimum-image sphere of `radius` around center_red.
    A_cols[:, i] = a_i in the same unit as radius. Loops over z planes to keep memory small.
    Returns (mean, npoints).
    """
    n1, n2, n3 = V.shape
    i1 = (np.arange(n1) / n1 - center_red[0] + 0.5) % 1.0 - 0.5
    i2 = (np.arange(n2) / n2 - center_red[1] + 0.5) % 1.0 - 0.5
    D1, D2 = np.meshgrid(i1, i2, indexing="ij")
    xy = np.stack([D1 * A_cols[0, 0] + D2 * A_cols[0, 1], D1 * A_cols[1, 0] + D2 * A_cols[1, 1],
                   D1 * A_cols[2, 0] + D2 * A_cols[2, 1]], axis=-1)          # (n1, n2, 3) in-plane part
    tot = 0.0; cnt = 0
    r2 = radius * radius
    for iz in range(n3):
        dz = (iz / n3 - center_red[2] + 0.5) % 1.0 - 0.5
        r = xy + dz * A_cols[:, 2][None, None, :]
        m = (r * r).sum(-1) < r2
        if m.any():
            tot += float(V[:, :, iz][m].sum()); cnt += int(m.sum())
    return (tot / cnt if cnt else np.nan), cnt


def far_atom_alignment(V_d, V_p, x_red_d, x_red_p, A_cols, radii):
    """
    Lu et al. (2019) alignment. V_d, V_p: local potentials [ix, iy, iz] on the same grid (same unit);
    x_red_*: reduced atomic positions of each cell; radii: iterable of sphere radii (unit of A_cols).
    The reference atom is the atom of the defective cell farthest from the vacancy; in the pristine cell
    the atom nearest to that position is used (identical position for an unrelaxed geometry).
    Returns dict(s_vac, i_far_d, i_far_p, dist_far, shifts {radius: <V_d> - <V_p>}, means {radius: (<V_d>, <V_p>, n_d, n_p)}).
    """
    s_vac, _, _ = vacancy_site(x_red_p, x_red_d, A_cols)
    i_d, dist = far_atom(x_red_d, s_vac, A_cols)
    c_d = np.mod(x_red_d[i_d], 1.0)
    i_p = int(np.argmin(min_image_dist(x_red_p, c_d, A_cols)))
    c_p = np.mod(x_red_p[i_p], 1.0)
    out = dict(s_vac=s_vac, i_far_d=i_d, i_far_p=i_p, dist_far=dist, shifts={}, means={},
               far_displacement=float(min_image_dist(x_red_p[i_p:i_p + 1], c_d, A_cols)[0]))
    for r in radii:
        md, nd = sphere_average(V_d, A_cols, c_d, r)
        mp, npts = sphere_average(V_p, A_cols, c_p, r)
        out["shifts"][r] = md - mp
        out["means"][r] = (md, mp, nd, npts)
    return out


def rigid_shift_fit(eps_ref, eps_model, e_lo):
    """
    Rigid shift `s` such that eps_model + s ~ eps_ref on the states below e_lo (in the reference).
    The n lowest reference states below e_lo are matched, in order, with the n lowest model states.
    Returns (shift = median(ref - model), max |residual|, n).
    """
    ref = np.sort(np.asarray(eps_ref)); ref = ref[ref < e_lo]
    mod = np.sort(np.asarray(eps_model))[:len(ref)]
    if len(ref) == 0 or len(mod) < len(ref):
        return np.nan, np.nan, 0
    d = ref - mod
    s = float(np.median(d))
    return s, float(np.abs(d - s).max()), int(len(ref))
