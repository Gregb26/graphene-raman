"""
qe_gamma_io.py
    Readers for Gamma-point pw.x runs written with the "gamma trick" (wfc attribute gamma_only = .TRUE.):
    the wfc*.hdf5 file stores only half of the G sphere (G = 0 and one member of each +-G pair) and the
    wavefunction is real in real space, psi(-G) = psi*(G). The production reader (qe_io._read_all_wfc)
    assumes a full sphere and one wfc file per k; it is left untouched. This module adds, for R4:

        read_wfc_gamma        -- one wfc file (wfc1 / wfcup1 / wfcdw1), a contiguous band slice
        get_eigenvalues_spin  -- KS eigenvalues per spin channel from data-file-schema.xml (LSDA block)
        mirror_parity_z       -- <psi| sigma_h |psi> from the plane-wave coefficients, sigma_h: z -> 2 z0 - z
        density_2d            -- |psi(r)|^2 summed over z on the in-plane FFT grid (gamma trick expanded)
        inplane_disc_mask     -- in-plane minimum-image disc around a reduced position (for local weights)

    Conventions: G . r = 2 pi (g . x_red), so the mirror phase for a plane at reduced height z0 is
    exp(-4 pi i g3 z0). Units: eigenvalues Hartree (as the XML), positions reduced.
"""
import os
import xml.etree.ElementTree as ET

import numpy as np
import h5py


def _floats(text):
    return np.array(text.split(), dtype=np.float64)


def wfc_path(save_dir, spin=None):
    """Path of the Gamma wfc file: spin None -> wfc1.hdf5, 'up'/1 -> wfcup1.hdf5, 'dw'/2 -> wfcdw1.hdf5."""
    name = {None: "wfc1.hdf5", "up": "wfcup1.hdf5", 1: "wfcup1.hdf5", "dw": "wfcdw1.hdf5", 2: "wfcdw1.hdf5"}[spin]
    p = os.path.join(save_dir, name)
    if not os.path.isfile(p):
        raise FileNotFoundError(p)
    return p


def read_wfc_gamma(save_dir, spin=None, bands=None):
    """
    Plane-wave coefficients of one Gamma-point wfc file.

    bands: None (all) or (b0, b1) -> contiguous slice [b0, b1) of band indices (0-based).
    Returns (C (nb, ngw) complex, mill (ngw, 3) int, gamma_only bool, attrs dict).
    The coefficients are those stored by QE (half sphere when gamma_only), scaled by scale_factor.
    """
    return read_wfc_file(wfc_path(save_dir, spin), bands=bands)


def read_wfc_file(path, bands=None):
    """Same as read_wfc_gamma for an explicit wfc<ik>.hdf5 path (also works for ordinary, full-sphere files)."""
    with h5py.File(path, "r") as f:
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
        gamma_only = str(attrs.get("gamma_only", ".FALSE.")).strip().upper() in (".TRUE.", "T", "TRUE")
        if int(attrs.get("npol", 1)) != 1:
            raise NotImplementedError("npol != 1 not supported")
        mill = np.asarray(f["MillerIndices"][:], dtype=np.int64)
        if bands is None:
            evc = np.asarray(f["evc"][:], dtype=np.float64)
        else:
            b0, b1 = int(bands[0]), int(bands[1])
            evc = np.asarray(f["evc"][b0:b1], dtype=np.float64)
        scale = float(attrs.get("scale_factor", 1.0))
    C = (evc[:, 0::2] + 1j * evc[:, 1::2]) * scale
    return C, mill, gamma_only, attrs


def get_eigenvalues_spin(save_dir):
    """
    KS eigenvalues (Hartree) as (nspin, nbnd) for a one-k-point (Gamma) run; nspin = 2 for LSDA, where the
    XML <eigenvalues> block of one <ks_energies> holds 2*nbnd values (spin up first, then spin down).
    Returns (eigs (nspin, nbnd), fermi (Hartree), lsda bool).
    """
    root = ET.parse(os.path.join(save_dir, "data-file-schema.xml")).getroot()
    bs = root.find(".//output/band_structure")
    lsda = (bs.find("lsda").text.strip().lower() == "true")
    fermi = float(bs.find("fermi_energy").text)
    ks = bs.findall("ks_energies")
    if len(ks) != 1:
        raise ValueError(f"{save_dir}: expected one k-point (Gamma), found {len(ks)}")
    e = _floats(ks[0].find("eigenvalues").text)
    if lsda:
        nb_up = int(bs.find("nbnd_up").text); nb_dw = int(bs.find("nbnd_dw").text)
        if nb_up != nb_dw or e.size != nb_up + nb_dw:
            raise ValueError(f"{save_dir}: LSDA eigenvalue block {e.size} vs nbnd_up/dw {nb_up}/{nb_dw}")
        return e.reshape(2, nb_up), fermi, True
    return e.reshape(1, -1), fermi, False


def _mill_keys(mill):
    """Injective integer key per Miller triplet (for lookups)."""
    off = int(np.abs(mill).max()) + 1
    B = 2 * off + 1
    return (mill[:, 0] + off) * B * B + (mill[:, 1] + off) * B + (mill[:, 2] + off), off, B


def _lookup(keys_sorted, order, query):
    """Index into the original array of each query key, -1 if absent."""
    pos = np.searchsorted(keys_sorted, query)
    pos = np.clip(pos, 0, len(keys_sorted) - 1)
    found = keys_sorted[pos] == query
    return np.where(found, order[pos], -1)


def mirror_parity_z(C, mill, z0_red=0.0, gamma_only=True):
    """
    <psi_b| sigma_h |psi_b> for each band b (rows of C), with sigma_h: z -> 2 z0 - z, z0 in reduced units
    of a3. In plane waves: (sigma psi)(G) = psi(g1, g2, -g3) exp(-4 pi i g3 z0), so
        <psi|sigma|psi> = sum_G C*(G) C(sigma G) exp(-4 pi i g3 z0).
    Gamma trick: the file holds half the sphere H; the -G terms are the complex conjugates of the +G terms,
    so the full-sphere sum is |C(0)|^2 + 2 Re sum_{G in H, G != 0} [...]. When sigma G is not stored,
    -sigma G is, and C(sigma G) = C*(-sigma G).
    Returns a complex array (nb,); for eigenstates of definite parity it is +-1 (real).
    """
    mill = np.asarray(mill, dtype=np.int64)
    keys, off, B = _mill_keys(mill)
    order = np.argsort(keys); ks = keys[order]
    sig = mill.copy(); sig[:, 2] = -sig[:, 2]
    kq = (sig[:, 0] + off) * B * B + (sig[:, 1] + off) * B + (sig[:, 2] + off)
    j = _lookup(ks, order, kq)
    phase = np.exp(-4j * np.pi * mill[:, 2] * float(z0_red))
    if gamma_only:
        kqn = (-sig[:, 0] + off) * B * B + (-sig[:, 1] + off) * B + (-sig[:, 2] + off)
        jn = _lookup(ks, order, kqn)
        missing = (j < 0) & (jn < 0)
        if missing.any():
            raise ValueError(f"{int(missing.sum())} mirrored G neither stored as +G nor -G (not a closed half sphere)")
        use_conj = j < 0
        jj = np.where(use_conj, jn, j)
        Cs = C[:, jj]
        Cs = np.where(use_conj[None, :], np.conj(Cs), Cs)
        w = np.where((mill == 0).all(axis=1), 1.0, 2.0)
        return np.sum(w[None, :] * np.real(np.conj(C) * Cs * phase[None, :]), axis=1).astype(complex)
    if (j < 0).any():
        raise ValueError(f"{int((j < 0).sum())} mirrored G not in the sphere (truncated sphere?)")
    return np.sum(np.conj(C) * C[:, j] * phase[None, :], axis=1)


def density_2d(C, mill, ngfft, gamma_only=True, workers=8):
    """
    |psi_b(r)|^2 summed over the third grid axis (z), on the (n1, n2) in-plane FFT grid, for each band.
    The full G sphere is rebuilt from the stored half (A[-G] = A*[G]) when gamma_only. Normalization:
    each map is divided by its own sum (fractions of the state).
    Returns rho2 (nb, n1, n2) float.
    """
    from scipy import fft as sfft
    n1, n2, n3 = (int(x) for x in ngfft)
    mill = np.asarray(mill, dtype=np.int64)
    idx = np.ravel_multi_index((mill[:, 0] % n1, mill[:, 1] % n2, mill[:, 2] % n3), (n1, n2, n3))
    nonzero = ~(mill == 0).all(axis=1)
    idxn = np.ravel_multi_index(((-mill[:, 0]) % n1, (-mill[:, 1]) % n2, (-mill[:, 2]) % n3), (n1, n2, n3))
    out = np.zeros((C.shape[0], n1, n2))
    A = np.zeros(n1 * n2 * n3, dtype=complex)
    for b in range(C.shape[0]):
        A[:] = 0.0
        A[idx] = C[b]
        if gamma_only:
            A[idxn[nonzero]] = np.conj(C[b][nonzero])
        psi = sfft.ifftn(A.reshape(n1, n2, n3), workers=workers)
        r2 = (psi.real ** 2 + psi.imag ** 2).sum(axis=2)
        out[b] = r2 / r2.sum()
    return out


def inplane_disc_mask(ngfft, A_cols, center_red, radius):
    """
    Boolean (n1, n2) mask of grid points whose in-plane minimum-image distance to center_red (reduced)
    is below radius (same length unit as the lattice vectors A_cols[:, i] = a_i).
    """
    n1, n2 = int(ngfft[0]), int(ngfft[1])
    i1 = (np.arange(n1) / n1 - center_red[0] + 0.5) % 1.0 - 0.5
    i2 = (np.arange(n2) / n2 - center_red[1] + 0.5) % 1.0 - 0.5
    D1, D2 = np.meshgrid(i1, i2, indexing="ij")
    x = D1 * A_cols[0, 0] + D2 * A_cols[0, 1]
    y = D1 * A_cols[1, 0] + D2 * A_cols[1, 1]
    return (x * x + y * y) < radius * radius
