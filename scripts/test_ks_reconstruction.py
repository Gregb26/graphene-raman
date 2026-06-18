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


def reconstruct_ks_hamiltonian(uc_save, pot_file, upf_file):
    """Return H_mn(k) (nk, nb, nb) and the QE eigenvalues eps (nb, nk), all in Hartree."""
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
    V, _ = qe_io.get_pot(pot_file, subtract_mean=False)
    V = V.transpose(2, 1, 0)
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
    uc_save = "data/graphene/unit_cell/qe/defect_5x5.save"
    sc_p_save = "data/graphene/supercell/qe/defect_5x5_p.save"
    pot_file = f"{uc_save}/Vks_uc"
    pot_p = f"{sc_p_save}/Vks_5x5_p"
    upf_file = f"{uc_save}/C.upf"

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

    ok = okA and okB
    print("\nRESULT:", "PASS" if ok else "FAIL", f"(tol={tol:.0e} Ha)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
