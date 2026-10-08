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


# ---------------------------------------------------------------------------------------------------------------- R11 (2026-10-08)
PLATEAU_FRAC = 0.75     # Kumagai-Oba constant: atoms at true distance >= PLATEAU_FRAC * r_max from the vacancy (R9 decision, R10 p1_size)


def plateau_stats(values, mask):
    """Mean, rms, max |deviation| and count of values[mask] (the plateau statistics of r9/r10_driver)."""
    v = np.asarray(values, float)[mask]
    mu = float(v.mean())
    return dict(mean=mu, rms=float(np.sqrt(((v - mu) ** 2).mean())), max_abs_dev=float(np.abs(v - mu).max()), n=int(mask.sum()))


def kumagai_oba_constant(V_d, V_p, x_red_d, x_red_p, A_cols, radius=1.0, frac=PLATEAU_FRAC):
    """
    C_N of config "alignment" (R10 p1_size, approximation (i)): mean of the sphere shifts of atom_sphere_shifts (radius, Angstrom)
    over the atoms of the defective cell at true minimum-image distance >= frac * r_max from the vacancy; rms over those atoms =
    the uncertainty of tab:alignement (column A). Units of V (production: eV).
    Returns dict(C_N, rms, max_abs_dev, n, r_max, d_min = frac * r_max, mask (n_d,), shifts = the atom_sphere_shifts dict).
    """
    r = atom_sphere_shifts(V_d, V_p, x_red_d, x_red_p, A_cols, radius)
    r_max = float(r["dist"].max())
    mask = r["dist"] >= frac * r_max
    st = plateau_stats(r["shift"], mask)
    return dict(C_N=st["mean"], rms=st["rms"], max_abs_dev=st["max_abs_dev"], n=st["n"], r_max=r_max, d_min=frac * r_max,
                mask=mask, shifts=r)


def plane_radial_profile(dV, x_red_d, s_vac, A_cols, r_core=0.5, dr=0.05, r_max=None):
    """
    Core-masked azimuthal average of dV in the graphene plane (fig:Ved (c), formerly inline in scripts/t/analyze_Ved.py).
    On the grid plane iz = round(s_vac[2] * n3), the points at in-plane distance >= r_core (Angstrom) from every atom of the
    defective cell are binned by their true minimum-image distance r to the vacancy (bins of width dr up to r_max, default
    |a1| / 2, the half-width of the supercell) and averaged bin by bin. The vacancy site itself is not masked.
    dV (n1, n2, n3) [ix, iy, iz] ; x_red_d (n_d, 3) ; s_vac (3,) reduced (vacancy_site) ; A_cols (3, 3) Angstrom.
    Returns dict(rc (nbin,) bin centres, rad (nbin,) mean of dV per bin (nan when empty), cnt (nbin,), r_max, iz, masked_fraction).
    """
    A = np.asarray(A_cols, float)
    s = np.asarray(s_vac, float)
    n1, n2, n3 = dV.shape
    iz = int(np.round(s[2] * n3)) % n3
    plane = dV[:, :, iz]
    I1, I2 = np.meshgrid(np.arange(n1) / n1, np.arange(n2) / n2, indexing="ij")
    R = true_min_image_dist(np.stack([I1.ravel(), I2.ravel(), np.full(I1.size, s[2])], axis=1), s, A).reshape(I1.shape)
    mask = np.ones(plane.shape, bool)
    for s_at in np.asarray(x_red_d, float):
        d1 = (I1 - s_at[0] + 0.5) % 1.0 - 0.5
        d2 = (I2 - s_at[1] + 0.5) % 1.0 - 0.5
        ra = np.sqrt((d1 * A[0, 0] + d2 * A[0, 1]) ** 2 + (d1 * A[1, 0] + d2 * A[1, 1]) ** 2)
        mask &= ra >= r_core
    if r_max is None:
        r_max = 0.5 * float(np.linalg.norm(A[:, 0]))
    edges = np.arange(0, r_max + dr, dr)
    m = (R < r_max) & mask
    cnt, _ = np.histogram(R[m], edges)
    sm, _ = np.histogram(R[m], edges, weights=plane[m])
    rad = np.where(cnt > 0, sm / np.maximum(cnt, 1), np.nan)
    return dict(rc=0.5 * (edges[1:] + edges[:-1]), rad=rad, cnt=cnt, r_max=float(r_max), iz=iz, masked_fraction=float(1 - mask.mean()))


def z_profile(dV, A_cols, x_red_d, z_vac=5.0):
    """
    Plane-averaged profile <dV>_xy(z) and the vacuum offset (R10 c6 ; tab:alignement column E). z is measured from the sheet
    (median reduced z of the atoms) and a3 must be along z (c = A_cols[2, 2]); the vacuum offset is the mean and std of the
    profile over the planes with |z| > z_vac (Angstrom), i.e. between z_vac and c / 2.
    dV (n1, n2, n3) [ix, iy, iz] ; A_cols (3, 3) Angstrom ; x_red_d (n_d, 3).
    Returns dict(z (n3,) Angstrom, profile (n3,), vacuum_mean, vacuum_std, vacuum_min, vacuum_max, n_vacuum, mean3d,
                 sheet_plane (value on the plane nearest to the sheet), c).
    """
    A = np.asarray(A_cols, float)
    if abs(A[0, 2]) > 1e-8 or abs(A[1, 2]) > 1e-8:
        raise ValueError(f"z_profile: a3 must be along z, got {A[:, 2]}")
    n3 = dV.shape[2]
    c = float(A[2, 2])
    prof = dV.mean(axis=(0, 1))
    z_sheet = float(np.mod(np.median(np.asarray(x_red_d, float)[:, 2]), 1.0))
    z = ((np.arange(n3) / n3 - z_sheet + 0.5) % 1.0 - 0.5) * c
    vac = np.abs(z) > z_vac
    if not vac.any():
        raise ValueError(f"z_profile: no plane with |z| > {z_vac} A (c / 2 = {c / 2:.2f} A)")
    return dict(z=z, profile=prof, vacuum_mean=float(prof[vac].mean()), vacuum_std=float(prof[vac].std()),
                vacuum_min=float(prof[vac].min()), vacuum_max=float(prof[vac].max()), n_vacuum=int(vac.sum()),
                mean3d=float(prof.mean()), sheet_plane=float(prof[int(np.argmin(np.abs(z)))]), c=c)


def potential_background(V_d, V_p, x_red_d, x_red_p, A_cols, radius=1.0, frac=PLATEAU_FRAC, r_core=0.5, dr=0.05, z_vac=5.0,
                         r_eval=(2.0, 3.0)):
    """
    The three measures of the background of the defect potential dV = V_d - V_p on one supercell pair (R11):
      (1) C_N, the Kumagai-Oba constant (kumagai_oba_constant: spheres of `radius`, atoms at >= frac * r_max) ;
      (2) the in-plane background: core-masked azimuthal profile of dV in the graphene plane (plane_radial_profile, r_core, dr),
          reduced to `plane` = mean of the bins at r >= frac * (|a1| / 2) (the plateau of fig:Ved (c) near the boundary, rms over
          those bins) and to its values interpolated at r_eval (the columns of sampling_table.csv) ;
      (3) the vacuum offset: plane-averaged profile along z (z_profile), mean and std over |z| > z_vac.
    V_d, V_p (n1, n2, n3) [ix, iy, iz] on the same grid, same unit (production: eV) ; x_red_*: reduced positions of each cell ;
    A_cols (3, 3) cell vectors as columns, Angstrom.
    Returns dict: C_N, C_N_rms, C_N_n, r_max ; plane, plane_rms, plane_n, plane_at {r: value} ; vacuum, vacuum_std, mean3d ;
    s_vac ; and the three full results kumagai_oba, radial, z.
    """
    if V_d.shape != V_p.shape:
        raise ValueError(f"potential_background: grids differ {V_d.shape} vs {V_p.shape}")
    ko = kumagai_oba_constant(V_d, V_p, x_red_d, x_red_p, A_cols, radius, frac)
    dV = V_d - V_p
    rp = plane_radial_profile(dV, x_red_d, ko["shifts"]["s_vac"], A_cols, r_core, dr)
    far = (rp["rc"] >= frac * rp["r_max"]) & np.isfinite(rp["rad"])
    pl = plateau_stats(rp["rad"], far)
    ok = np.isfinite(rp["rad"])
    at = {float(r0): float(np.interp(r0, rp["rc"][ok], rp["rad"][ok])) for r0 in r_eval}
    zp = z_profile(dV, A_cols, x_red_d, z_vac)
    return dict(C_N=ko["C_N"], C_N_rms=ko["rms"], C_N_n=ko["n"], r_max=ko["r_max"],
                plane=pl["mean"], plane_rms=pl["rms"], plane_n=pl["n"], plane_at=at,
                vacuum=zp["vacuum_mean"], vacuum_std=zp["vacuum_std"], mean3d=zp["mean3d"],
                s_vac=ko["shifts"]["s_vac"], kumagai_oba=ko, radial=rp, z=zp)
