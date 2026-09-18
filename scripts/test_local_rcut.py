"""
Synthetic validation of the R_cut machinery: (i) with R_local covering the defect support the local
rate is EXACT (matches dense compute_T); (ii) truncating R_local changes it (so R_cut is a real
convergence parameter); (iii) extract_V_loc / mwr_locality guardrails behave.
"""
import numpy as np

from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.defects.many_body.single_defect import compute_T
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from test_local_tmatrix import random_H   # reuse the synthetic Hamiltonian


def main():
    nw, Nc, eta = 3, 5, 0.10
    Hwr, Rw, ndeg = random_H(Hwr_nw := nw)
    R_full = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])   # defect on 5 sites
    nL = len(R_full)
    rng = np.random.default_rng(2)
    v = rng.normal(size=(nL * nw, nL * nw)) + 1j * rng.normal(size=(nL * nw, nL * nw))
    v = 0.5 * (v + v.conj().T)
    # make it physically LOCALIZED: on-site block dominant, neighbor blocks decaying (guardrail a)
    s = np.exp(-0.6 * np.linalg.norm(R_full, axis=1))          # (nL,) site weight, 1 at R0
    scale = np.outer(s, s)                                     # symmetric -> preserves Hermiticity
    V_full = (v.reshape(nL, nw, nL, nw) * scale[:, None, :, None]).reshape(nL * nw, nL * nw)

    # build a synthetic Mwr (nw, nr, nw, nr) whose R_full x R_full block equals V_full (index L*nw+w)
    R_mwr = R_full.copy()
    Mwr = np.zeros((nw, nL, nw, nL), complex)
    Vb = V_full.reshape(nL, nw, nL, nw)              # [L,w,L',w']
    Mwr[:, :, :, :] = np.transpose(Vb, (1, 0, 3, 2))  # (w,L,w',L')
    # guardrail helpers
    dist, wt = lt.mwr_locality(Mwr, R_mwr)
    print("[locality] ||Mwr(R,R0)|| vs |R-R0|:", dict(zip(np.round(dist, 3), np.round(wt, 3))))
    V_ex, res = lt.extract_V_loc(Mwr, R_mwr, R_full)
    print(f"[extract] V_loc matches V_full ? {np.allclose(V_ex, V_full, atol=1e-12)}  herm_res={res:.1e}")

    # dense reference from V_full (multi-R projection)
    kc = lt.mp_grid(Nc, Nc, 1)
    _, Ec, Uc = Hwr_to_Hwk(Hwr, Rw, kc, ndegen=ndeg)
    Nk = len(kc)
    ph = lt._phase(kc, R_full)                       # (Nk, nL)
    phi = np.einsum("kL,kwn->knLw", ph, Uc, optimize=True).reshape(Nk, nw, nL * nw)
    M_bk = np.einsum("kna,ab,KMb->nkMK", np.conj(phi), V_full, phi, optimize=True) / Nk
    eigs = Ec.T
    Gd = np.zeros((nw, Nk))
    for ik in range(Nk):
        for n in range(nw):
            T = compute_T(np.array([Ec[ik, n]]), eigs, M_bk, eta)
            Gd[n, ik] = -2.0 * T[0, n, ik, n, ik].imag
    Gd_perdef = Gd * Nk

    # local: full R_cut (exact) vs truncated to {0}
    Gl_full = lt.scattering_rate(Hwr, Rw, ndeg, V_full, R_full, kc, eta, k_int=kc)
    V0, _ = lt.extract_V_loc(Mwr, R_mwr, np.array([[0, 0, 0]]))
    Gl_trunc = lt.scattering_rate(Hwr, Rw, ndeg, V0, np.array([[0, 0, 0]]), kc, eta, k_int=kc)

    rel_full = np.max(np.abs(Gl_full - Gd_perdef)) / max(1e-30, np.max(np.abs(Gd_perdef)))
    rel_trunc = np.max(np.abs(Gl_trunc - Gd_perdef)) / max(1e-30, np.max(np.abs(Gd_perdef)))
    print(f"[R_cut full]  rel err vs dense = {rel_full:.2e} (seuil 1e-8, support complet = exact) : {'PASS' if rel_full < 1e-8 else 'FAIL'}")
    print(f"[R_cut {{0}}]   rel diff vs dense = {rel_trunc:.2e}   (should be O(1): truncation matters)")
    ok = rel_full < 1e-8 and rel_trunc > 1e-3 and np.allclose(V_ex, V_full, atol=1e-12)
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
