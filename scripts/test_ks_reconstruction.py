"""
test_ks_reconstruction.py
    Standalone validation of the QE I/O + matrix-element machinery (no external reference needed).

    The unit-cell Kohn-Sham wavefunctions diagonalise the unit-cell KS Hamiltonian, so

        H_mn(k) = T_mn(k) + <psi_mk|V_loc|psi_nk> + V^NL_mn(k)  ==  eps_nk * delta_mn

    where eps_nk are the KS eigenvalues written by QE in data-file-schema.xml. Reconstructing
    H from the pipeline pieces and recovering eps (diagonal, no constant offset) validates, in
    one shot and against a known answer:
        - wavefunction reading C_nk(G) and Miller indices         (qe_io.get_C_nk / get_G_red)
        - kinetic energy matrix T_mn(k) and the k+G vectors
        - the local KS potential                                  (qe_io.get_pot: Ry->Ha, G=0 ref, axes)
        - the UPF non-local projectors                            (read_upf, Hankel F_il, KB energies,
                                                                   spherical harmonics, phases, 4pi/sqrt(Omega))

    This is the same machinery used by compute_M_L and compute_M_NL for the electron-defect matrix.
"""

import numpy as np
from scipy.interpolate import CubicSpline

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.pseudo_io import read_upf, fq_from_fr
from electron_defect_interaction.utils.planewaves import mask_invalid_G
from electron_defect_interaction.utils.lattice import red_to_cart
from electron_defect_interaction.defects.non_local import (
    build_K_vectors, compute_phase, compute_angular_part, compute_M_NL,
)
from electron_defect_interaction.defects.local_R import compute_ML_R
from electron_defect_interaction.wavefunctions.wfk import compute_psi_nk


def sc_pot_on_uc_grid(uc_save, sc_save, pot_sc_file):
    """
    Restrict a *pristine* supercell local KS potential onto the unit-cell FFT grid, with no
    averaging/folding -- it is fed directly.

    A defect-free N1xN2xN3 supercell is a periodic tiling of the unit cell, so its local KS
    potential V_p is already periodic with the unit-cell period. The supercell FFT grid here has the
    same spacing as the unit-cell grid (QE 5x5 graphene: 150 = 5*30 in-plane, 192 = 1*192 out of
    plane), so the unit-cell grid points coincide with the first N1xN2xN3 block of the supercell
    grid. We therefore take that block of V_p directly -- no fold, no average. Because V_p is
    unit-cell-periodic, <psi_mk^uc | V_p | psi_nk^uc> is diagonal and equals <psi|V_uc|psi> up to a
    constant offset.

    Returns
    -------
        V_uc: (nr1, nr2, nr3) array, [ix,iy,iz], Hartree -- V_p restricted to the unit-cell grid.
        spread: float -- max|V_p(image) - V_p(image 0)| over the N1xN2xN3 unit-cell images, the
                deviation from perfect unit-cell periodicity (a diagnostic that the premise holds and
                that the two QE runs are mutually consistent; NOT used in the reconstruction).
    """
    A_uc, _ = qe_io.get_A_volume(uc_save)
    A_sc, _ = qe_io.get_A_volume(sc_save)
    N = np.rint(np.linalg.inv(A_uc) @ A_sc).astype(int)  # supercell multiplicity
    if not np.allclose(N, np.diag(np.diag(N))):
        raise ValueError(f"supercell is not a diagonal repeat of the unit cell:\n{N}")
    N1, N2, N3 = np.diag(N)

    n1, n2, n3 = qe_io.get_ngfft(uc_save)
    ng_sc = qe_io.get_ngfft(sc_save)
    for a, (nsc, nuc, Na) in enumerate(zip(ng_sc, (n1, n2, n3), (N1, N2, N3))):
        if nsc != Na * nuc:
            raise ValueError(f"axis {a}: sc grid {nsc} != {Na} x uc grid {nuc} "
                             f"(grids share spacing only if commensurate)")

    V_sc, _ = qe_io.get_pot(pot_sc_file, subtract_mean=False)
    V_sc = V_sc.transpose(2, 1, 0)              # [ix,iy,iz] on the supercell grid
    V_uc = V_sc[:n1, :n2, :n3]                  # image 0 == unit-cell potential (fed directly)

    blocks = V_sc.reshape(N1, n1, N2, n2, N3, n3)
    spread = np.max(np.abs(blocks - blocks[0, :, 0, :, 0, :][None, :, None, :, None, :]))

    return V_uc, spread


def reconstruct_ks_hamiltonian(uc_save, pot_file, upf_file, V_override=None):
    """
    Return H_mn(k) (nk, nb, nb) and the QE eigenvalues eps (nb, nk), all in Hartree.

    If V_override is given (a unit-cell-grid potential [ix,iy,iz] in Hartree, e.g. the pristine
    supercell potential restricted by sc_pot_on_uc_grid), it is used for the local term instead of
    reading pot_file.
    """
    C_nkg, nG = qe_io.get_C_nk(uc_save)
    G_red = qe_io.get_G_red(uc_save)
    k_red = qe_io.get_k_red(uc_save)
    B_uc, _ = qe_io.get_B_volume(uc_save)
    A_uc, Om = qe_io.get_A_volume(uc_save)
    ecut = float(qe_io.get_ecut(uc_save))
    tau = red_to_cart(qe_io.get_x_red(uc_save), A_uc)
    eps = qe_io.get_eigenvalues(uc_save)
    ngfft = qe_io.get_ngfft(uc_save)
    nb, nk, _ = C_nkg.shape

    keep = mask_invalid_G(nG)
    C = np.where(keep, C_nkg, 0.0)

    # Kinetic: T_mn(k) = 0.5 sum_G C*_mk(G) |k+G|^2 C_nk(G)
    K = red_to_cart(k_red[:, None, :] + G_red, B_uc)
    K2 = np.where(keep, np.sum(K**2, axis=2), 0.0)
    T = np.zeros((nk, nb, nb), complex)
    for ik in range(nk):
        Ck = C[:, ik, :]
        T[ik] = 0.5 * (Ck.conj() * K2[ik]) @ Ck.T

    # Non-local: V^NL_mn(k) = (4pi)^2/Omega sum_{a,l,i,ml} E_li B*_mk B_nk
    ekb_li, fr_li, rgrid, lmax, imax, _ = read_upf(upf_file)
    Kc, Kn, Kh = build_K_vectors(k_red, G_red, keep, B_uc)
    q = np.linspace(0, 2 * np.sqrt(2 * ecut), 2000)
    Fq = CubicSpline(q, fq_from_fr(rgrid, fr_li, q), axis=-1, extrapolate=False)
    F = Fq(Kn)
    phase = compute_phase(Kc, tau)
    Y = compute_angular_part(Kh, lmax)
    pref = 4 * np.pi / np.sqrt(Om)
    Bp = pref * np.einsum("nkg,likg,kglm,ksg->nkslim", np.conj(C), F, Y, phase, optimize=True)
    NL = np.einsum("li,pkslia,qkslia->kpq", ekb_li, Bp, np.conj(Bp), optimize=True)

    # Local: L_mn(k) = <psi_mk|V_loc|psi_nk>  (keep the full mean, no subtraction)
    if V_override is None:
        V, _ = qe_io.get_pot(pot_file, subtract_mean=False)
        V = V.transpose(2, 1, 0)
    else:
        V = V_override
    psi, _ = compute_psi_nk(C_nkg, nG, G_red, k_red, Om, ngfft=ngfft)
    Nr = np.prod(ngfft); dV = Om / Nr
    Vr = V.reshape(Nr)
    L = np.zeros((nk, nb, nb), complex)
    for ik in range(nk):
        P = psi[:, ik].reshape(nb, Nr)
        L[ik] = dV * (P.conj() * Vr) @ P.T

    return T + L + NL, eps


def null_test(uc_save, sc_p_save, pot_p, upf_file):
    """
    Identity/null check of the defect construction: if the 'defect' equals the pristine reference,
    V_ed = V_p - V_p = 0 and M_d - M_p = 0, so both M^L and M^NL must vanish exactly. This tests the
    difference construction itself (independently of the KS reconstruction above).
    """
    M_L = compute_ML_R(uc_save, sc_p_save, pot_p, pot_p,
                       subtract_mean=False, bands=[0, 1, 2, 3], io=qe_io)
    M_NL = compute_M_NL(uc_save, sc_p_save, sc_p_save, upf_file,
                        io=qe_io, pseudo_reader=read_upf)
    return np.max(np.abs(M_L)), np.max(np.abs(M_NL))


def main():
    import _paths
    grid = "5x5"
    uc_save = _paths.uc(grid)
    sc_p_save = _paths.sc_p(grid)
    pot_file = _paths.uc_pot(grid)
    pot_p = _paths.pot_p(grid)
    upf_file = _paths.upf(grid)

    tol = 1e-5

    # Test A: unit-cell KS Hamiltonian reconstruction == QE eigenvalues
    H, eps = reconstruct_ks_hamiltonian(uc_save, pot_file, upf_file)
    nk, nb, _ = H.shape
    eye = np.einsum("ki,ij->kij", np.einsum("kii->ki", H), np.eye(nb))
    offdiag_max = np.max(np.abs(H - eye))
    diag = np.real(np.einsum("kii->ki", H)).T
    shift = np.mean(diag - eps)
    diag_dev = np.max(np.abs(diag - eps - shift))
    okA = offdiag_max < tol and abs(shift) < tol and diag_dev < tol

    print(f"=== Test A: KS Hamiltonian reconstruction (bands={nb}, k={nk}) ===")
    print(f"  off-diagonal max      : {offdiag_max:.2e} Ha")
    print(f"  diag(H)-eps offset    : {shift:.2e} Ha")
    print(f"  max|diag(H)-eps-shift|: {diag_dev:.2e} Ha   -> {'PASS' if okA else 'FAIL'}")

    # Test B: null defect (defect = pristine) -> M^L = M^NL = 0
    mL, mNL = null_test(uc_save, sc_p_save, pot_p, upf_file)
    okB = mL < tol and mNL < tol
    print(f"\n=== Test B: null defect (defect = pristine) ===")
    print(f"  max|M^L|={mL:.2e}  max|M^NL|={mNL:.2e}   -> {'PASS' if okB else 'FAIL'}")

    # Test C: feed the PRISTINE SUPERCELL potential V_p (unchanged, not folded) into the unit-cell
    # reconstruction. V_p is periodic with the unit-cell period, so <psi^uc|V_p|psi^uc> is diagonal
    # and T + L + NL must still recover the unit-cell eigenvalues -- validating get_pot on the
    # supercell and the consistency of the pristine-supercell and unit-cell QE runs.
    V_uc, spread = sc_pot_on_uc_grid(uc_save, sc_p_save, pot_p)
    Hc, epsc = reconstruct_ks_hamiltonian(uc_save, pot_file, upf_file, V_override=V_uc)
    eyec = np.einsum("ki,ij->kij", np.einsum("kii->ki", Hc), np.eye(nb))
    offdiag_maxC = np.max(np.abs(Hc - eyec))
    diagC = np.real(np.einsum("kii->ki", Hc)).T
    shiftC = np.mean(diagC - epsc)
    diag_devC = np.max(np.abs(diagC - epsc - shiftC))
    okC = offdiag_maxC < tol and diag_devC < tol
    print(f"\n=== Test C: pristine supercell V_p -> unit-cell eigenvalues ===")
    print(f"  V_p unit-cell periodicity spread: {spread:.2e} Ha")
    print(f"  off-diagonal max      : {offdiag_maxC:.2e} Ha")
    print(f"  diag(H)-eps offset    : {shiftC:.2e} Ha")
    print(f"  max|diag(H)-eps-shift|: {diag_devC:.2e} Ha   -> {'PASS' if okC else 'FAIL'}")

    ok = okA and okB and okC
    print("\nRESULT:", "PASS" if ok else "FAIL", f"(tol={tol:.0e} Ha)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
