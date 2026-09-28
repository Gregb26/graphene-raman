"""
wannier_io.py
    Python module containing functions that extracts information from Wannier90 output files.
"""

import numpy as np

def read_w90_mat(w90_path):
    """
    Extracts the unitary matrix U, from a Wanier90 .mat outputfile, needed to go from the Bloch basis to the Wannier basis. 

    Inputs:
        w90_path: str,
            Path to Wannier90 .mat output file.
    Returns:
        mat: (nkpt, nband, nwann) array of complex
            Unitary U transformation matrix to go from Bloch basis to Wannier basis. If U is not a disentanglement matrix: nwann=nband. If 
            U is a disentanglement matrix: nwann < nband.
        k_red: (nkpt, 3) array of floats
            kpoints in reduced coordinates. Should (must) match the DFT k-points (e.g. qe_io.get_k_red).
    """

    from pathlib import Path

    with Path(w90_path).open() as f:

        date = f.readline()

        nkpt, nwann, nband = map(int, f.readline().split() )# number of kpoints, number of wannier functions, number of bloch bands. For U, Nw = Nb

        f.readline()

        mat = np.zeros((nkpt, nband, nwann), dtype=complex)
        k_red = np.zeros((nkpt, 3), dtype=float)

        # nkpt blocks with first row being the kpt the next N*J lines being real and imag part of the matrix in column major order
        for ik in range(nkpt):

            k_red[ik,:] = np.array([float(kpt) for kpt in f.readline().split()])

            block = np.zeros((nband*nwann), dtype=complex)
            for i in range(nwann*nband):
                Re, Im = map(float, f.readline().split())
                block[i] = Re + 1j*Im
            
            mat[ik] = block.reshape((nband,nwann), order='F')
            f.readline()


    # testing unitarity
    for ik in range(mat.shape[0]):
       assert np.allclose(mat[ik].conj().T @ mat[ik], np.eye(mat[ik].shape[1]), atol=1e-8), 'U matrix must be unitary'

    return mat, k_red

def check_hermicity_HR(HR, R):
    """
    Checks the Hermicity of the Hamiltonian written in the Wannier tight binding basis. The condition is
    H(R) = H^\\dag(-R).
    Inputs
        HR: (nrpts, nw, nw) array of complex
            Hamiltonien written in Wannier tight binding basis
        R: (nrpts, 3) array of floats
            R vectors computed by Wannier90 used in the tight binding basis
    """

    # Map R = (n1, n2, n3) to index in array R 
    R_tuples = [tuple(r) for r in R]
    R_index  = {r: i for i, r in enumerate(R_tuples)}

    for i, Rs in enumerate(R):
        R_ = tuple(-Rs) # -R

        # some R may be unpaired, skip those
        j = R_index.get(R_)
        if j is None:
            continue
        H = HR[i]
        H_h = HR[j].conj().T

        if not np.allclose(H, H_h, atol=1e-10):
            print(f"Hermiticity fails for R = {R_}, -R index = {j}")
            raise TypeError('Hamiltonian not Hermitian')
    
    return True

def read_w90_hr(w90_path):
    """
    Read a Wannier90 `seedname_hr.dat` file (the tight-binding Hamiltonian H(R) in the Wannier basis).

    Format: header line, nw, nrpts, the Wigner-Seitz degeneracy list (nrpts ints, 15 per line), then
    nrpts blocks of nw*nw lines "Rx Ry Rz m n Re Im".

    Returns
    -------
        HR: (nrpts, nw, nw) array of complex   -- H(R) in the Wannier basis
        R:  (nrpts, 3) array of ints           -- R vectors (reduced coords)
        ndegen: (nrpts,) array of ints         -- Wigner-Seitz degeneracies; H(k)=sum_R e^{ik.R} H(R)/ndegen(R)
    """
    from pathlib import Path
    with Path(w90_path).open() as f:
        f.readline()                       # header (date)
        nw = int(f.readline())
        nrpts = int(f.readline())

        ndeg = []
        while len(ndeg) < nrpts:           # degeneracy list, 15 per line
            ndeg += [int(x) for x in f.readline().split()]
        ndegen = np.array(ndeg[:nrpts], dtype=int)

        HR = np.zeros((nrpts, nw, nw), dtype=complex)
        R = np.zeros((nrpts, 3), dtype=int)
        for ir in range(nrpts):
            for _ in range(nw * nw):
                t = f.readline().split()
                rx, ry, rz, m, n = (int(x) for x in t[:5])
                HR[ir, m - 1, n - 1] = float(t[5]) + 1j * float(t[6])
            R[ir] = (rx, ry, rz)

    check_hermicity_HR(HR, R)
    return HR, R, ndegen

def _read_tb_blocks(f, nrpts, nw, ncomp):
    """
    Read one section of a `_tb.dat` (H: ncomp = 1, r: ncomp = 3) from the open file `f`: per R, a blank
    line, "R1 R2 R3", then nw*nw lines "m n Re Im ..." holding <0m|X|Rn> (m = row). Shared by H and r,
    hence one index convention.

    Returns
    -------
        X_R: (nrpts, ncomp, nw, nw) array of complex
        R:   (nrpts, 3) array of ints
    """
    X_R = np.zeros((nrpts, ncomp, nw, nw), dtype=complex)
    R = np.zeros((nrpts, 3), dtype=int)
    for ir in range(nrpts):
        f.readline()                       # blank line
        R[ir] = [int(n) for n in f.readline().split()]
        for _ in range(nw * nw):
            t = f.readline().split()
            m, n = int(t[0]) - 1, int(t[1]) - 1
            v = np.array(t[2:], dtype=float)                # Re, Im of each component, interleaved
            X_R[ir, :, m, n] = v[0::2] + 1j * v[1::2]       # shape error if not 2*ncomp values
    return X_R, R

def read_w90_tb(w90_path):
    """
    Read a Wannier90 `seedname_tb.dat` (format: memoire/EM/EM.md, section 2). Values are raw: no
    division by ndegen, r not hermitized; the Hermiticity of H is checked.

    Returns
    -------
        HR: (nrpts, nw, nw) array of complex     -- <0m|H|Rn>, eV
        R:  (nrpts, 3) array of ints             -- R vectors (reduced coords)
        ndegen: (nrpts,) array of ints           -- Wigner-Seitz degeneracies
        rR: (nrpts, 3, nw, nw) array of complex  -- <0m|r_alpha|Rn>, Angstrom; diagonal of r(0) = centres
        lattice: (3, 3) array of floats          -- Angstrom, lattice[:, i] = a_i
    """
    from pathlib import Path
    with Path(w90_path).open() as f:

        f.readline()                       # header (date)

        # lattice vectors in Angstrom, one per line, stored in columns
        a = [[float(x) for x in f.readline().split()] for _ in range(3)]
        lattice = np.column_stack(a) # lattice[:,i] = a_i

        nw = int(f.readline())             # number of wannier functions
        nrpts = int(f.readline())          # number of inequivalent R points used by Wannier90

        ndeg = []                          # Wigner-Seitz degeneracy list, 15 per line
        while len(ndeg) < nrpts:
            ndeg += [int(x) for x in f.readline().split()]
        ndegen = np.array(ndeg[:nrpts], dtype=int)

        HR, R = _read_tb_blocks(f, nrpts, nw, ncomp=1)
        rR, R_r = _read_tb_blocks(f, nrpts, nw, ncomp=3)

    assert np.array_equal(R_r, R), 'H and r sections of the _tb.dat have different R lists'
    HR = HR[:, 0]                          # (nrpts, nw, nw)
    check_hermicity_HR(HR, R)
    return HR, R, ndegen, rR, lattice

def read_w90_HR(w90_path):
    """Former name of `read_w90_tb`, returning (HR, R, ndegen). Kept for the frozen campaign drivers."""
    HR, R, ndegen, _, _ = read_w90_tb(w90_path)
    return HR, R, ndegen
