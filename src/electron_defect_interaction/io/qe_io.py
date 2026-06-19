"""
qe_io.py
    Helper functions that extract the quantities needed by the matrix-element code from Quantum
    ESPRESSO (pw.x) outputs.

    A QE calculation writes everything we need inside the `prefix.save/` directory:
        - data-file-schema.xml : geometry, reciprocal lattice, ecut, k-points, atoms, eigenvalues, FFT grid
        - wfc<ik>.hdf5         : plane-wave coefficients C_nk(G) and Miller indices G, one file per k-point

    Conventions (verified against a graphene unit-cell run, QE 7.5):
        - The XML is written in *Hartree atomic units* (Units="Hartree atomic units"):
          energies in Ha, lengths in Bohr, ecutwfc in Ha.
        - Plane-wave coefficients are normalized so that sum_G |C_nk(G)|^2 = 1.
        - Miller indices are signed integer reduced coordinates of G, usable directly with
          utils.fft_utils.map_G_to_fft_grid (signed-integer FFT convention).

    All functions take the path to the `prefix.save/` directory.
"""

import os
import re
import glob
import xml.etree.ElementTree as ET

import numpy as np
import h5py

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

def _xml_path(save_dir):
    """Return the path to data-file-schema.xml inside a QE prefix.save directory."""
    p = os.path.join(save_dir, "data-file-schema.xml")
    if not os.path.isfile(p):
        raise FileNotFoundError(f"No data-file-schema.xml found in {save_dir!r}")
    return p

def _root(save_dir):
    """Parse and return the XML root element of data-file-schema.xml."""
    return ET.parse(_xml_path(save_dir)).getroot()

def _floats(text):
    """Parse a whitespace-separated block of Fortran-style floats into a 1d array."""
    return np.array(text.split(), dtype=np.float64)

def _wfc_files(save_dir):
    """
    Return the list of wfc<ik>.hdf5 files in `save_dir`, ordered by the QE k-point index `ik`
    (matching the order of the <ks_energies> blocks in the XML).
    """
    files = glob.glob(os.path.join(save_dir, "wfc*.hdf5"))
    if not files:
        raise FileNotFoundError(f"No wfc*.hdf5 files found in {save_dir!r}")

    def ik_of(path):
        with h5py.File(path, "r") as f:
            return int(f.attrs["ik"])

    return sorted(files, key=ik_of)

# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------

def get_A_volume(save_dir):
    """
    Extract the primitive direct lattice vectors and the cell volume.

    Returns
    -------
        A: (3,3) array of floats
            Primitive direct lattice vectors as columns A[:,i] = a_i, in Bohr.
        Omega: float
            Cell volume in Bohr^3.
    """
    root = _root(save_dir)
    cell = root.find(".//output/atomic_structure/cell")
    a1 = _floats(cell.find("a1").text)
    a2 = _floats(cell.find("a2").text)
    a3 = _floats(cell.find("a3").text)

    A = np.column_stack((a1, a2, a3))  # columns are a_i
    Omega = np.abs(np.linalg.det(A))

    return A, Omega

def get_B_volume(save_dir):
    """
    Extract the primitive reciprocal lattice vectors and reciprocal cell volume.

    Returns
    -------
        B: (3,3) array of floats
            Reciprocal primitive lattice vectors as columns B[:,i] = b_i, in 1/Bohr,
            satisfying a_i . b_j = 2*pi*delta_ij.
        Omega_G: float
            Reciprocal cell volume.
    """
    A, Omega = get_A_volume(save_dir)
    a1, a2, a3 = A[:, 0], A[:, 1], A[:, 2]

    b1 = 2 * np.pi * np.cross(a2, a3) / Omega
    b2 = 2 * np.pi * np.cross(a3, a1) / Omega
    b3 = 2 * np.pi * np.cross(a1, a2) / Omega

    B = np.column_stack((b1, b2, b3))
    Omega_G = np.abs(np.linalg.det(B))

    return B, Omega_G

def get_ecut(save_dir):
    """
    Plane-wave energy cutoff used in the calculation, in Hartree.

    Note: the QE XML stores ecutwfc in Hartree (Units="Hartree atomic units"), so the value is
    directly usable -- no Ry->Ha conversion is needed.
    """
    root = _root(save_dir)
    ecut = float(root.find(".//output/basis_set/ecutwfc").text)
    return np.asarray(ecut, dtype=float)

def get_ngfft(save_dir, smooth=False):
    """
    FFT grid used by QE, read directly from the XML.

    QE writes the exact FFT mesh, so downstream wavefunction reconstruction and
    local-potential handling should use this grid rather than fft_utils.fft_grid_from_G_red.

    Inputs
    ------
        smooth: bool
            If True return the smooth wavefunction grid <fft_smooth>, else the dense grid <fft_grid>.

    Returns
    -------
        ngfft: tuple (nr1, nr2, nr3) of ints
    """
    root = _root(save_dir)
    tag = "fft_smooth" if smooth else "fft_grid"
    g = root.find(f".//output/basis_set/{tag}")
    return (int(g.get("nr1")), int(g.get("nr2")), int(g.get("nr3")))

def get_x_red(save_dir):
    """
    Reduced (fractional) coordinates of the atoms.

    QE stores Cartesian positions (Bohr) in the XML; we convert to reduced coordinates so that
    red_to_cart(x_red, A) recovers the Cartesian positions.

    Returns
    -------
        x_red: (natom, 3) array of floats
    """
    root = _root(save_dir)
    A, _ = get_A_volume(save_dir)
    Ainv = np.linalg.inv(A)  # columns of A are a_i -> x_red = Ainv @ tau_cart

    pos = []
    for atom in root.findall(".//output/atomic_structure/atomic_positions/atom"):
        tau_cart = _floats(atom.text)
        pos.append(Ainv @ tau_cart)

    return np.asarray(pos, dtype=np.float64)

def get_typat(save_dir):
    """
    Integer atom-type index for each atom, 1-based.

    Returns
    -------
        typat: (natom,) array of ints
    """
    root = _root(save_dir)

    # Map species name -> 1-based type index, in declaration order
    species = [s.get("name") for s in root.findall(".//output/atomic_species/species")]
    type_of = {name: i + 1 for i, name in enumerate(species)}

    typat = [type_of[atom.get("name")]
             for atom in root.findall(".//output/atomic_structure/atomic_positions/atom")]

    return np.asarray(typat, dtype=int)

# -----------------------------------------------------------------------------
# k-points and eigenvalues
# -----------------------------------------------------------------------------

def get_k_red(save_dir):
    """
    k-points in reduced (fractional) coordinates.

    QE stores k-points in Cartesian coordinates in units of 2*pi/alat. We convert to reduced
    coordinates k_red such that k_cart = B @ k_red (B columns are b_i), i.e. k_red = inv(B) @ k_cart.

    Returns
    -------
        k_red: (nkpt, 3) array of floats
    """
    root = _root(save_dir)
    B, _ = get_B_volume(save_dir)
    Binv = np.linalg.inv(B)

    # alat: QE Cartesian k units are 2*pi/alat, with alat = |a1|.
    alat = float(root.find(".//output/atomic_structure").get("alat"))
    scale = 2 * np.pi / alat

    k_red = []
    for ks in root.findall(".//output/band_structure/ks_energies"):
        k_cart = _floats(ks.find("k_point").text) * scale  # 1/Bohr
        k_red.append(Binv @ k_cart)

    return np.asarray(k_red, dtype=np.float64)

def get_qe_bands(bands_dat_path):
    """
    Read a bands.x `bands.dat` file (raw format).

    Format: a header line '&plot nbnd=.., nks=.. /', then for each k-point a line with the Cartesian
    k (units 2*pi/alat) followed by nbnd eigenvalues in eV.

    Returns
    -------
        kcart: (nks, 3) array of floats -- Cartesian k in units of 2*pi/alat
        eig:   (nks, nbnd) array of floats -- eigenvalues in eV
    """
    with open(bands_dat_path) as f:
        header = f.readline()
        nbnd = int(re.search(r"nbnd=\s*(\d+)", header).group(1))
        nks = int(re.search(r"nks=\s*(\d+)", header).group(1))
        data = np.array(f.read().split(), dtype=float).reshape(nks, 3 + nbnd)
    return data[:, :3], data[:, 3:]

def kcart_to_kred(kcart, save_dir):
    """
    Convert QE Cartesian k (units 2*pi/alat, e.g. from get_qe_bands) to reduced coordinates, mirroring
    get_k_red (k_cart = B @ k_red, with alat = |a1|).

    Inputs
    ------
        kcart: (..., 3) array of floats -- Cartesian k in units of 2*pi/alat
        save_dir: str -- a prefix.save dir (for the reciprocal lattice and alat)

    Returns
    -------
        k_red: (..., 3) array of floats
    """
    B, _ = get_B_volume(save_dir)
    alat = float(_root(save_dir).find(".//output/atomic_structure").get("alat"))
    return (np.asarray(kcart) * (2 * np.pi / alat)) @ np.linalg.inv(B).T

def get_eigenvalues(save_dir, shift_Fermi=False):
    """
    Kohn-Sham eigenvalues in Hartree.

    Returns
    -------
        eigenvalues: (nband, nkpt) array of floats
    """
    root = _root(save_dir)
    bs = root.find(".//output/band_structure")
    fermi = float(bs.find("fermi_energy").text)

    eigs = []
    for ks in bs.findall("ks_energies"):
        e = _floats(ks.find("eigenvalues").text)
        eigs.append(e)

    eigs = np.asarray(eigs, dtype=np.float64).T  # (nband, nkpt)
    if shift_Fermi:
        eigs = eigs - fermi

    return eigs

# -----------------------------------------------------------------------------
# Plane-wave coefficients and G-vectors (wfc<ik>.hdf5)
# -----------------------------------------------------------------------------

def _read_all_wfc(save_dir):
    """
    Read all wfc<ik>.hdf5 files and assemble padded arrays.

    QE stores one file per k-point, each with its own number of plane waves (igwx) and no padding.
    We pad to nG_max = max_k igwx so the returned arrays match the (.., nkpt, nG_max, ..) shapes the
    downstream code expects, together with the per-k active count nG (number of plane waves per k).

    Returns
    -------
        C_nkg: (nband, nkpt, nG_max) complex
        nG:    (nkpt,) int
        G_red: (nkpt, nG_max, 3) int
    """
    files = _wfc_files(save_dir)
    nkpt = len(files)

    # First pass: shapes
    nbnd = None
    nG = np.zeros(nkpt, dtype=np.int64)
    for ik, path in enumerate(files):
        with h5py.File(path, "r") as f:
            nG[ik] = int(f.attrs["igwx"])
            if int(f.attrs.get("npol", 1)) != 1:
                raise NotImplementedError("qe_io currently supports npol=1 (no spinors) only.")
            nb = int(f.attrs["nbnd"])
            nbnd = nb if nbnd is None else nbnd
            if nb != nbnd:
                raise ValueError("Inconsistent nbnd across wfc files.")

    nG_max = int(nG.max())
    C_nkg = np.zeros((nbnd, nkpt, nG_max), dtype=np.complex128)
    G_red = np.zeros((nkpt, nG_max, 3), dtype=np.int64)

    # Second pass: fill
    for ik, path in enumerate(files):
        with h5py.File(path, "r") as f:
            nG_k = int(f.attrs["igwx"])
            scale = float(f.attrs.get("scale_factor", 1.0))

            mill = np.asarray(f["MillerIndices"][:], dtype=np.int64)  # (nG_k, 3)
            evc = np.asarray(f["evc"][:], dtype=np.float64)           # (nbnd, 2*nG_k)

        # interleaved real/imag -> complex
        C = (evc[:, 0::2] + 1j * evc[:, 1::2]) * scale               # (nbnd, nG_k)

        C_nkg[:, ik, :nG_k] = C
        G_red[ik, :nG_k, :] = mill

    return C_nkg, nG, G_red

def get_pot(filepath, subtract_mean=True, to_hartree=True):
    """
    Load the local Kohn-Sham potential V(r) from a Quantum ESPRESSO `pp.x` plot_file (filplot),
    e.g. produced with plot_num=1 (V_bare + V_Hartree + V_xc).

    QE's pp.x writes its plot data in *Rydberg*. We convert to Hartree by default so the
    resulting M^L is in the Hartree atomic units used throughout the pipeline.

    The returned array is laid out so that the existing caller pattern
        V = get_pot(...)[0].transpose(2, 1, 0)
    yields V[ix, iy, iz] aligned with psi[ix, iy, iz] from compute_psi_nk (x = a1 fastest).

    Inputs
    ------
        filepath: str
            Path to the pp.x plot_file (filplot), iflag=3 / 3D dump.
        subtract_mean: bool
            If True remove the constant (G=0) offset before returning (recommended).
        to_hartree: bool
            If True (default) convert from Rydberg to Hartree (factor 1/2).

    Returns
    -------
        V: (nr3, nr2, nr1) array of floats
            Local KS potential. Caller should `.transpose(2,1,0)` to get [ix,iy,iz] (see above).
        ngfft: tuple (nr1, nr2, nr3)
            FFT grid size.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    # line 0: title ; line 1: nrx1 nrx2 nrx3 nr1 nr2 nr3 nat ntyp
    hdr = lines[1].split()
    nr1, nr2, nr3 = int(hdr[3]), int(hdr[4]), int(hdr[5])
    nat, ntyp = int(hdr[6]), int(hdr[7])

    # line 2: ibrav celldm(1:6)
    ibrav = int(lines[2].split()[0])

    cur = 3
    cur += 3 if ibrav == 0 else 0   # explicit lattice vectors only when ibrav==0
    cur += 1                        # gcutm, dual, ecut, plot_num
    cur += ntyp                     # atom-type table
    cur += nat                      # atomic positions

    # remaining tokens are the data, written x-fastest (Fortran ir = ix + iy*nr1 + iz*nr1*nr2)
    data = np.fromstring(" ".join(lines[cur:]), sep=" ", dtype=np.float64)
    expected = nr1 * nr2 * nr3
    if data.size != expected:
        raise ValueError(f"plot_file data size {data.size} != nr1*nr2*nr3 = {expected}")

    # x-fastest flat -> [iz, iy, ix] via C-order reshape (nr3, nr2, nr1)
    V = data.reshape(nr3, nr2, nr1)

    if to_hartree:
        V = 0.5 * V                 # Rydberg -> Hartree

    if subtract_mean:
        V = V - V.mean()

    return V, (nr1, nr2, nr3)

def get_C_nk(save_dir):
    """
    Plane-wave coefficients C_nk(G).

    Returns
    -------
        C_nk: (nband, nkpt, nG_max) complex
        npw:  (nkpt,) int -- number of active plane waves per k-point
    """
    C_nkg, nG, _ = _read_all_wfc(save_dir)
    return C_nkg, nG

def get_G_red(save_dir):
    """
    Reciprocal lattice vectors G in reduced (signed integer) coordinates.

    Returns
    -------
        G_red: (nkpt, nG_max, 3) int
    """
    _, _, G_red = _read_all_wfc(save_dir)
    return G_red
