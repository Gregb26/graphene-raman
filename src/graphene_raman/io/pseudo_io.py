"""
pseudo_io.py
    Python module containing functions to extract the relevant information from .psp8 pseudopotential files. Currently extracts
    the nonlocal part of the pseudopotential, KB energies and KB radial projectors in real space, and transforms the radial
    projectors to reciprocal space via a Hankel transformation. To be used to evaluate non-local matrix elements.

    TODO
        - Generalize get_psps to handle different different kind of .psp8 files (currentely tailored to carbon only)
        - Could maybe add a module to get_psps to read the local part of the pseudopotential as well
"""

import numpy as np

def read_upf(path, to_hartree=True):
    """
    Read a UPF (v2) norm-conserving pseudopotential and extract the non-local Kleinman-Bylander data
    in the layout expected by fq_from_fr / compute_M_NL (the layout of the former ABINIT psp8 reader).

    UPF stores:
        - PP_MESH/PP_R   : radial grid r (size = mesh_size)
        - PP_BETA.i      : the radial KB projectors, tabulated as r*beta_i(r) (r*beta convention of the radial arrays)
        - PP_DIJ         : the D_ij matrix; for ONCV norm-conserving it is diagonal, D_ii = KB energy of projector i
        - PP_LOCAL       : the local potential
    Each PP_BETA.i carries its angular_momentum l; projectors are grouped by l into channels i.

    Units: UPF energies (PP_DIJ, PP_LOCAL) are in Rydberg; the code expects Hartree, so we
    convert by 1/2 when to_hartree=True (default). The projectors r*beta(r) are dimensionful but the same
    convention as .psp8 (verified by overlaying the Hankel form factors of C.psp8 and C.upf).

    Returns:
        ekb_li: (lmax+1, imax) array of floats   -- KB energies (Hartree)
        fr_li:  (lmax+1, imax, mmax) array       -- r*beta radial projectors on the real-space grid
        rgrid:  (mmax,) array                    -- radial grid
        lmax:   int
        imax:   int
        V_L:    (mmax,) array                    -- local potential (Hartree)
    """
    import xml.etree.ElementTree as ET

    root = ET.parse(path).getroot()

    def arr(tag_path):
        el = root.find(tag_path)
        return np.array(el.text.split(), dtype=np.float64)

    rgrid = arr(".//PP_MESH/PP_R")
    mmax = rgrid.size
    V_L = arr(".//PP_LOCAL")

    nl = root.find(".//PP_NONLOCAL")

    # Collect beta projectors with their angular momentum, in file order
    betas = []  # list of (l, values)
    for el in nl:
        if el.tag.startswith("PP_BETA"):
            l = int(el.get("angular_momentum"))
            vals = np.array(el.text.split(), dtype=np.float64)
            if vals.size < mmax:  # pad to full grid if stored shorter
                vals = np.concatenate([vals, np.zeros(mmax - vals.size)])
            betas.append((l, vals[:mmax]))

    # D_ij matrix (flat -> square), diagonal entries are the per-projector KB energies
    nproj = len(betas)
    Dij = arr(".//PP_DIJ").reshape(nproj, nproj)
    D = np.diag(Dij)

    lmax = max(l for l, _ in betas)
    # number of channels per l
    nproj_l = np.zeros(lmax + 1, dtype=int)
    for l, _ in betas:
        nproj_l[l] += 1
    imax = int(nproj_l.max())

    ekb_li = np.zeros((lmax + 1, imax))
    fr_li = np.zeros((lmax + 1, imax, mmax))

    ch = np.zeros(lmax + 1, dtype=int)  # current channel index per l
    for i, (l, vals) in enumerate(betas):
        c = ch[l]
        fr_li[l, c, :] = vals
        ekb_li[l, c] = D[i]
        ch[l] += 1

    if to_hartree:
        ekb_li = 0.5 * ekb_li
        V_L = 0.5 * V_L

    return ekb_li, fr_li, rgrid, lmax, imax, V_L

def fq_from_fr(r, fr_li, q):
    """
    Transforms the radial form factors fr_li in real space to reciprocal space on a q grid. This is done via a Hankel transofrmation:
    fq_li = int dr r fr_il(r) j_l(qr), where j_l(qr) is the spherical Bessel function of order l.

    Inputs:
        r: (mmax, ) array floats
            Grid in real space
        fr_li: (lmax+1, imax, mmax) array of floats
            Radial form factors on a 1d grid of size mmax in real space.
        q: (mqff) array of floats:
            1d grid in recirprocal space. qmax is chosen such that qmax = 2*sqrt(2*ecut) to make sure to resolve all K=|k+G| vectors
        
    Returns:
        fq_li: (lmax+1, imax, mqff) array of floats
            Radial form factors on a 1d grid of size mqff in reciprocal space. To be interpolated and evaluated at K=|k+G| vectors.
    """
    from scipy.special import spherical_jn
    from scipy.integrate import simpson

    lmax = fr_li.shape[0] - 1
    imax = fr_li.shape[1] 
    mqff = q.size

    qr = r[:, np.newaxis] * q[np.newaxis, :] # (mmax, mqff)
    fq_li = np.zeros((lmax+1, imax, mqff))

    for l in range(lmax + 1):
       jl = spherical_jn(l, qr) # (mmax, mqff)
       fr_i = fr_li[l, :, :] # (imax, mmax)
       integrand = (fr_i * r[np.newaxis, :])[:, :, np.newaxis] *jl[np.newaxis, :, :] # (imax, mmax, mqff)

       fq_li[l, :, :] =  simpson(integrand, x=r, axis=1) # integrate over r
    
    return fq_li # (lmax+1, imax, mqff)

