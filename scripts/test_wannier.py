"""
test_wannier.py
    Standalone validation of the Wannier-interpolation pipeline for the electron-defect matrix M,
    using the QE-derived Wannier90 files in data/graphene/unit_cell/qe/defect_5x5.save/
    (wannier_u.mat, wannier_u_dis.mat, wannier_tb.dat, wannier_hr.dat). graphene 5x5: nb=16, nw=5.

    Tests, in increasing scope:
      (1) PARSERS
          - read_w90_mat(u.mat)     : U is unitary           U^dag U = U U^dag = I_nw
          - read_w90_mat(u_dis.mat) : U_dis is an isometry    U_dis^dag U_dis = I_nw, and
                                       P = U_dis U_dis^dag is a rank-nw projector (P^2=P, P^dag=P, trP=nw)
          - V = U_dis @ U is an isometry (V^dag V = I_nw)
          - read_w90_HR(tb.dat)     : Hermiticity H(R)=H(-R)^dag (reader asserts it), ndegen>0,
                                       R contains the origin; H(k) Hermitian with real eigenvalues
          - read_w90_hr(hr.dat) agrees with read_w90_HR(tb.dat) on H(R) and ndegen
      (2) ROUND-TRIP (Wannier-gauge Fourier identity, EXACT)
          - M_wk -> Mwk_to_Mwr -> Mwr_to_Mwk == M_wk to machine precision
      (3) FULL PIPELINE (gauge invariant, EXACT on coarse->coarse)
          - the Bloch back-rotation is a per-k unitary similarity and the FT round-trip is the identity,
            so eigvals(wannier_interpolate(M, k_coarse, k_coarse)) == eigvals(Mbk_to_Mwk(M, U, U_dis))
      (4) FINE GRID (smoke test)
          - interpolating to a denser grid yields a Hermitian operator and real, finite band energies
"""

import numpy as np

from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR, read_w90_hr
from electron_defect_interaction.wannier.wannier_interpolation import (
    Mbk_to_Mwk, Mwk_to_Mwr, Mwr_to_Mwk, wannier_interpolate,
)
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk

DATA = "data/graphene/unit_cell/qe/defect_5x5.save"
U_PATH = f"{DATA}/wannier_u.mat"
UDIS_PATH = f"{DATA}/wannier_u_dis.mat"
TB_PATH = f"{DATA}/wannier_tb.dat"
HR_PATH = f"{DATA}/wannier_hr.dat"


def random_hermitian_M(nb, nk, seed):
    """Random M[b,k,B,K] that is Hermitian as an operator on the (band,k) index: M[b,k,B,K]=conj(M[B,K,b,k])."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((nb * nk, nb * nk)) + 1j * rng.standard_normal((nb * nk, nb * nk))
    H = 0.5 * (A + A.conj().T)
    return H.reshape(nb, nk, nb, nk)


def mp_grid(N1, N2, N3):
    """Unshifted Gamma-centered Monkhorst-Pack grid, reduced coords, in the i->j->l nested order."""
    return np.array([[i / N1, j / N2, l / N3]
                     for i in range(N1) for j in range(N2) for l in range(N3)], dtype=float)


def test_parsers():
    print("=== (1) parser / matrix-property tests ===")
    ok = True

    # --- U matrix: unitary ---
    U, k_U = read_w90_mat(U_PATH)                       # reader already asserts U^dag U = I
    nk, nw, _ = U.shape
    e1 = max(np.max(np.abs(U[ik].conj().T @ U[ik] - np.eye(nw))) for ik in range(nk))
    e2 = max(np.max(np.abs(U[ik] @ U[ik].conj().T - np.eye(nw))) for ik in range(nk))
    okU = e1 < 1e-7 and e2 < 1e-7
    ok &= okU
    print(f"  U   (nk={nk}, nw={nw}): max|U^dU-I|={e1:.1e}  max|UU^d-I|={e2:.1e}  -> {'PASS' if okU else 'FAIL'}")

    # --- U_dis matrix: isometry + projector ---
    Ud, k_Ud = read_w90_mat(UDIS_PATH)                  # reader already asserts U_dis^dag U_dis = I_nw
    _, nb, nwd = Ud.shape
    iso = max(np.max(np.abs(Ud[ik].conj().T @ Ud[ik] - np.eye(nwd))) for ik in range(nk))
    pdev, trdev = 0.0, 0.0
    for ik in range(nk):
        P = Ud[ik] @ Ud[ik].conj().T                    # (nb, nb)
        pdev = max(pdev, np.max(np.abs(P @ P - P)), np.max(np.abs(P - P.conj().T)))
        trdev = max(trdev, abs(np.trace(P).real - nwd))
    okUd = iso < 1e-7 and pdev < 1e-6 and trdev < 1e-7
    ok &= okUd
    print(f"  Udis(nb={nb}, nw={nwd}): max|U^dU-I|={iso:.1e}  max|P^2-P|,|P-P^d|={pdev:.1e}  |trP-nw|={trdev:.1e}  -> {'PASS' if okUd else 'FAIL'}")

    # --- V = U_dis @ U : isometry ---
    perm = np.array([np.argmin(np.sum(np.abs(np.mod(k_U - k, 1.0) - np.round(np.mod(k_U - k, 1.0))), axis=1))
                     for k in k_Ud])                    # align U to U_dis k-order (identity if same order)
    viso = 0.0
    for ik in range(nk):
        V = Ud[ik] @ U[perm[ik]]                        # (nb, nw)
        viso = max(viso, np.max(np.abs(V.conj().T @ V - np.eye(nw))))
    okV = viso < 1e-7
    ok &= okV
    print(f"  V=Udis@U          : max|V^dV-I|={viso:.1e}  -> {'PASS' if okV else 'FAIL'}")

    # --- H(R) tight-binding from tb.dat ---
    Hwr, Rw, ndegen = read_w90_HR(TB_PATH)              # reader asserts H(R)=H(-R)^dag
    nrpts = Hwr.shape[0]
    has_origin = np.any(np.all(Rw == 0, axis=1))
    okH = (ndegen.min() > 0) and has_origin and (len(ndegen) == nrpts == len(Rw))
    # H(k) Hermitian, eigenvalues real
    Hwk, Ewk, _ = Hwr_to_Hwk(Hwr, Rw, k_U, ndegen=ndegen)
    herm = max(np.max(np.abs(Hwk[ik] - Hwk[ik].conj().T)) for ik in range(nk))
    okH = okH and herm < 1e-10 and np.max(np.abs(Ewk.imag)) < 1e-12
    ok &= okH
    print(f"  H(R) tb.dat (nrpts={nrpts}): ndegen>0 & origin & lengths={okH or 'see'};  max|H(k)-H(k)^d|={herm:.1e}  -> {'PASS' if okH else 'FAIL'}")

    # --- hr.dat agrees with tb.dat ---
    HR2, R2, nd2 = read_w90_hr(HR_PATH)
    idx = {tuple(r): i for i, r in enumerate(R2)}
    dH, dN = 0.0, 0
    for i, r in enumerate(Rw):
        j = idx[tuple(r)]
        dH = max(dH, np.max(np.abs(Hwr[i] - HR2[j])))
        dN = max(dN, abs(int(ndegen[i]) - int(nd2[j])))
    okHR = dH < 1e-5 and dN == 0                        # hr.dat is written with ~6 decimals
    ok &= okHR
    print(f"  hr.dat vs tb.dat  : max|dH(R)|={dH:.1e}  max|d ndegen|={dN}  -> {'PASS' if okHR else 'FAIL'}")

    return ok, U, Ud, k_U


def test_roundtrip(k_coarse):
    print("\n=== (2) round-trip Wannier-gauge Fourier identity (M_wk -> M_wr -> M_wk) ===")
    nw, nk = 5, len(k_coarse)
    from electron_defect_interaction.wannier.wannier_interpolation import _infer_mp_grid
    MP = _infer_mp_grid(k_coarse)
    Mwk = random_hermitian_M(nw, nk, seed=1)            # random Hermitian object in the Wannier gauge
    Mwr, R = Mwk_to_Mwr(Mwk, k_coarse, MP)
    Mwk2 = Mwr_to_Mwk(Mwr, R, k_coarse)
    err = np.max(np.abs(Mwk2 - Mwk))
    ok = err < 1e-10
    print(f"  MP grid={MP}, nk={nk}, nr={len(R)}: max|M_wk' - M_wk|={err:.2e}  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_pipeline(k_coarse, U, Ud):
    print("\n=== (3) full pipeline, gauge-invariant spectrum (coarse -> coarse) ===")
    nk = len(k_coarse)
    nb = Ud.shape[1]; nw = Ud.shape[2]
    Mbk = random_hermitian_M(nb, nk, seed=2)            # random Hermitian M in the ab-initio Bloch gauge

    # reference: spectrum of M rotated into the Wannier gauge (the disentangled subspace)
    Mwk_ref = Mbk_to_Mwk(Mbk, U, Ud)
    ev_ref = np.linalg.eigvalsh(Mwk_ref.reshape(nw * nk, nw * nk))

    # pipeline interpolated back onto the same coarse grid, in the smooth Bloch gauge
    Mbk_fine, E_fine = wannier_interpolate(Mbk, k_coarse, k_coarse, TB_PATH, U_PATH, UDIS_PATH)
    op = Mbk_fine.reshape(nw * nk, nw * nk)
    herm = np.max(np.abs(op - op.conj().T))
    ev_pipe = np.linalg.eigvalsh(op)

    spec_err = np.max(np.abs(np.sort(ev_pipe) - np.sort(ev_ref)))
    ok = spec_err < 1e-8 and herm < 1e-8
    print(f"  nb={nb}->nw={nw}, nk={nk}: max|M_bk^fine - h.c.|={herm:.2e}")
    print(f"  spectrum: max|eig(pipeline) - eig(Wannier gauge)|={spec_err:.2e}  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_fine_grid(k_coarse, Ud):
    print("\n=== (4) interpolation to a finer grid (smoke test) ===")
    nb, nw = Ud.shape[1], Ud.shape[2]
    Mbk = random_hermitian_M(nb, len(k_coarse), seed=3)
    k_fine = mp_grid(10, 10, 1)                         # denser Gamma-centered grid
    Mbk_fine, E_fine = wannier_interpolate(Mbk, k_coarse, k_fine, TB_PATH, U_PATH, UDIS_PATH)
    op = Mbk_fine.reshape(nw * len(k_fine), nw * len(k_fine))
    herm = np.max(np.abs(op - op.conj().T))
    ok = (Mbk_fine.shape == (nw, len(k_fine), nw, len(k_fine))
          and np.isfinite(Mbk_fine).all() and np.isfinite(E_fine).all()
          and herm < 1e-8)
    print(f"  k_fine: {len(k_fine)} pts -> M_bk^fine {Mbk_fine.shape}, E_fine {E_fine.shape}")
    print(f"  Hermitian operator: max|M-M^d|={herm:.2e}, finite={np.isfinite(Mbk_fine).all()}  -> {'PASS' if ok else 'FAIL'}")
    return ok


def test_real_matrix(matrix_path="M_ed.npy"):
    """Round-trip on the actual electron-defect matrix M = M^L + M^NL (M_ed.npy), if present."""
    print("\n=== (5) round-trip on the REAL electron-defect matrix (M_ed.npy) ===")
    import os
    from electron_defect_interaction.io import qe_io
    from electron_defect_interaction.wannier.wannier_interpolation import _match_kpoint_order, _infer_mp_grid

    if not os.path.exists(matrix_path):
        print(f"  {matrix_path} not found -> SKIP")
        return True

    M = np.load(matrix_path)                                # (nb, nk, nb, nk), M was built in get_k_red order
    nb, nk = M.shape[0], M.shape[1]
    U, kU = read_w90_mat(U_PATH)
    Ud, kUd = read_w90_mat(UDIS_PATH)
    nw = U.shape[1]
    k_coarse = qe_io.get_k_red(DATA)                        # the k-order in which M was computed
    if (nb, nk) != (Ud.shape[1], len(k_coarse)):
        print(f"  shapes incompatible (M {M.shape}, U_dis {Ud.shape}, nk={len(k_coarse)}) -> SKIP")
        return True
    Uo = U[_match_kpoint_order(kU, k_coarse)]
    Udo = Ud[_match_kpoint_order(kUd, k_coarse)]

    # Hermiticity of the physical operator over (band, k)
    op = M.reshape(nb * nk, nb * nk)
    herm = np.max(np.abs(op - op.conj().T))

    # (a) exact FT round-trip on the *physical* Wannier-gauge matrix (element-wise)
    Mwk = Mbk_to_Mwk(M, Uo, Udo)
    Mwr, R = Mwk_to_Mwr(Mwk, k_coarse, _infer_mp_grid(k_coarse))
    Mwk2 = Mwr_to_Mwk(Mwr, R, k_coarse)
    scale = max(np.max(np.abs(Mwk)), 1e-30)
    ft_rel = np.max(np.abs(Mwk2 - Mwk)) / scale

    # (b) full pipeline coarse->coarse, gauge-invariant spectrum
    ev_ref = np.linalg.eigvalsh(Mwk.reshape(nw * nk, nw * nk))
    Mbk_fine, _ = wannier_interpolate(M, k_coarse, k_coarse, TB_PATH, U_PATH, UDIS_PATH)
    ev_pipe = np.linalg.eigvalsh(Mbk_fine.reshape(nw * nk, nw * nk))
    spec_rel = np.max(np.abs(np.sort(ev_pipe) - np.sort(ev_ref))) / max(np.max(np.abs(ev_ref)), 1e-30)

    ok = ft_rel < 1e-10 and spec_rel < 1e-8
    print(f"  M_ed shape={M.shape}: max|M|={np.max(np.abs(M)):.3e} Ha, max|M-M^d|={herm:.2e} Ha")
    print(f"  (a) FT round-trip on physical M_wk : max|M_wk'-M_wk|/scale = {ft_rel:.2e}")
    print(f"  (b) full pipeline spectrum          : max|eig(pipe)-eig(M_wk)|/scale = {spec_rel:.2e}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ok1, U, Ud, k_coarse = test_parsers()
    ok2 = test_roundtrip(k_coarse)
    ok3 = test_pipeline(k_coarse, U, Ud)
    ok4 = test_fine_grid(k_coarse, Ud)
    ok5 = test_real_matrix()
    ok = ok1 and ok2 and ok3 and ok4 and ok5
    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
