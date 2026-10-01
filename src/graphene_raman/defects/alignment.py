"""
alignment.py
    Potential alignment between a defective and a pristine supercell (R4, 2026-09-25).

    far_atom_alignment -- single-site alignment (Lu et al. 2019): average of the local KS potential (pp.x plot_num=1)
                          in a sphere around the atom farthest from the vacancy, in the defective and in the
                          pristine supercell; the difference is the rigid shift of the defective energy reference.
    rigid_shift_fit    -- rigid shift between two spectra, fitted on the states below an energy (outside the window).
    atom_sphere_shifts -- (R9, 2026-09-28; source of the Kumagai-Oba C_N) the same sphere-averaged shift around EVERY atom of the defective cell, with its true
                          minimum-image distance to the vacancy: profile Delta_a(d) of V_d - V_p away from the defect.
    Units: potentials and energies as given (the caller converts), positions reduced, radii in the unit of A.

    Note (R9): min_image_dist (and therefore far_atom) reduces each reduced coordinate to [-0.5, 0.5) separately. In a cell with
    a1, a2 at 60 degrees this is not always the shortest image (9x9: atom 161 is at 18.51 A per axis, 11.12 A true minimum image).
    Sphere averages (radius << cell) are unaffected; distances to the vacancy are. atom_sphere_shifts reports both.
"""
import itertools

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


# ---------------------------------------------------------------------------------------------------------------- R9 (P1)
def true_min_image_dist(x_red, center_red, A_cols):
    """
    True minimum-image Cartesian distances of reduced positions x_red (n, 3) to center_red: the per-axis reduced offset d in
    [-0.5, 0.5)^3 is completed by a search over the 27 neighbouring images d + s, s in {-1, 0, 1}^3, and the shortest is kept
    (exact for any cell whose per-axis reduction lands within one image of the shortest vector, e.g. a1, a2 at 60 degrees).
    """
    d = (np.asarray(x_red, float) - np.asarray(center_red, float) + 0.5) % 1.0 - 0.5            # (n, 3) per-axis reduction
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)), float)                      # (27, 3)
    r = (d[:, None, :] + shifts[None, :, :]) @ np.asarray(A_cols, float).T                        # (n, 27, 3) Cartesian
    return np.linalg.norm(r, axis=2).min(axis=1)


def _sphere_points(center_red, A_cols, shape, radius):
    """
    Flat indices of the grid points inside the sphere of `radius` around center_red, with EXACTLY the criterion of sphere_average
    (per-axis minimum-image offsets, same arithmetic order, strict inequality), evaluated only on a local index box around the
    centre instead of the whole grid. The box half-width along axis i is ceil(radius * ||row i of A^-1|| * n_i) + 1 points: a point
    inside the sphere has |d_i| <= radius * ||row i of A^-1|| (d = A^-1 r), so no point of the sphere is missed. An axis whose box
    would cover the whole period is taken whole (no duplicates).
    """
    A = np.asarray(A_cols, float); n = np.asarray(shape, int); c = np.asarray(center_red, float)
    reach = radius * np.linalg.norm(np.linalg.inv(A), axis=1)                                     # max |d_i| inside the sphere
    axes = []
    for i in range(3):
        h = int(np.ceil(reach[i] * n[i])) + 1
        ic = int(np.rint(c[i] * n[i]))
        idx = np.arange(n[i]) if 2 * h + 1 >= n[i] else np.mod(np.arange(ic - h, ic + h + 1), n[i])
        axes.append(idx)
    I1, I2, I3 = np.meshgrid(*axes, indexing="ij")
    d1 = (I1 / n[0] - c[0] + 0.5) % 1.0 - 0.5
    d2 = (I2 / n[1] - c[1] + 0.5) % 1.0 - 0.5
    d3 = (I3 / n[2] - c[2] + 0.5) % 1.0 - 0.5
    # same order of operations as sphere_average: in-plane part first, then + dz * a3, then sum of squares
    rx = (d1 * A[0, 0] + d2 * A[0, 1]) + d3 * A[0, 2]
    ry = (d1 * A[1, 0] + d2 * A[1, 1]) + d3 * A[1, 2]
    rz = (d1 * A[2, 0] + d2 * A[2, 1]) + d3 * A[2, 2]
    m = (rx * rx + ry * ry + rz * rz) < radius * radius
    return np.ravel_multi_index((I1[m], I2[m], I3[m]), tuple(n))


def atom_sphere_shifts(V_d, V_p, x_red_d, x_red_p, A_cols, radius):
    """
    Sphere-averaged potential shift around every atom of the defective cell (R9, A.1 (i); generalises far_atom_alignment).
    Its mean over the atoms at true distance >= 0.75 r_max is the Kumagai-Oba constant C_N (config block "alignment").

        Delta_a = <V_d>_{S(tau_a, rho)} - <V_p>_{S(tau'_p(a), rho)},
        S(tau, rho) = { grid points r : |r - tau|_(per-axis minimum image) < rho }   (3D sphere, all z planes),
        <V>_S = (1/|S|) sum_{r in S} V(r)   (arithmetic mean of the grid values),

    tau_a = position of atom a in the defective cell, p(a) = pristine atom nearest to it (true minimum image; same position for an
    unrelaxed geometry). The point set S is exactly the one of sphere_average, so the value at far_atom reproduces
    far_atom_alignment (to the summation order, ~1e-16 relative).

    Inputs:
        V_d, V_p: (n1, n2, n3) local potentials [ix, iy, iz] on the same grid (e.g. pp.x plot_num=1, eV).
        x_red_d, x_red_p: (n_d, 3), (n_p, 3) reduced atomic positions of the defective and pristine cells.
        A_cols: (3, 3) cell vectors as columns (A_cols[:, i] = a_i), in the unit of radius (Angstrom).
        radius: sphere radius (Angstrom).
    Returns dict:
        s_vac (3,)          reduced position of the vacancy (vacancy_site);
        dist (n_d,)         true minimum-image distance of each defective atom to the vacancy (unit of A_cols);
        dist_axis (n_d,)    per-axis minimum-image distance (min_image_dist, the convention of far_atom);
        shift (n_d,)        Delta_a (unit of V);  mean_d, mean_p (n_d,) the two sphere means;
        npts_d, npts_p (n_d,) number of grid points in each sphere;  i_p (n_d,) index of the pristine partner;
        i_far_axis          index of far_atom (the atom used by far_atom_alignment).
    """
    x_d = np.asarray(x_red_d, float); x_p = np.asarray(x_red_p, float); A = np.asarray(A_cols, float)
    shape = V_d.shape
    if V_p.shape != shape:
        raise ValueError(f"atom_sphere_shifts: grids differ {V_d.shape} vs {V_p.shape}")
    s_vac, _, _ = vacancy_site(x_p, x_d, A)
    n_d = len(x_d)
    out = dict(s_vac=s_vac, dist=true_min_image_dist(x_d, s_vac, A), dist_axis=min_image_dist(x_d, s_vac, A),
               shift=np.zeros(n_d), mean_d=np.zeros(n_d), mean_p=np.zeros(n_d),
               npts_d=np.zeros(n_d, int), npts_p=np.zeros(n_d, int), i_p=np.zeros(n_d, int), i_far_axis=far_atom(x_d, s_vac, A)[0])
    Vd = V_d.reshape(-1); Vp = V_p.reshape(-1)
    for a in range(n_d):
        c_d = np.mod(x_d[a], 1.0)
        ip = int(np.argmin(true_min_image_dist(x_p, c_d, A)))
        c_p = np.mod(x_p[ip], 1.0)
        pd = _sphere_points(c_d, A, shape, radius)
        pp = pd if np.array_equal(c_d, c_p) else _sphere_points(c_p, A, shape, radius)
        out["mean_d"][a] = Vd[pd].mean() if len(pd) else np.nan
        out["mean_p"][a] = Vp[pp].mean() if len(pp) else np.nan
        out["npts_d"][a], out["npts_p"][a], out["i_p"][a] = len(pd), len(pp), ip
    out["shift"] = out["mean_d"] - out["mean_p"]
    return out
