"""
Golden test (synthetic, no Wannier files needed): the local t-matrix scattering rate must reproduce
the dense single_defect.compute_T rate on the SAME coarse grid with full R_cut, to noise. Also checks
positivity (Gamma>=0) and Hermiticity of V_loc. Blocking before production.
"""
import numpy as np

from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.defects.many_body.single_defect import compute_T
from electron_defect_interaction.defects.many_body import local_tmatrix as lt


def random_H(nw, seed=0):
    """A random Hermitian tight-binding H(R) on R in {0,+-x,+-y}, with H(-R)=H(R)^dag."""
    rng = np.random.default_rng(seed)
    Rpos = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    Rw, Hwr = [], []
    for R in Rpos:
        h = rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))
        if R == (0, 0, 0):
            h = 0.5 * (h + h.conj().T)          # on-site Hermitian
            Rw.append(R); Hwr.append(h)
        else:
            Rw.append(R); Hwr.append(h)
            Rw.append(tuple(-np.array(R))); Hwr.append(h.conj().T)   # ensure H(k) Hermitian
    return np.array(Hwr), np.array(Rw), np.ones(len(Rw))


def main():
    nw, Nc, eta = 3, 5, 0.10
    Hwr, Rw, ndeg = random_H(nw)
    R_local = np.array([[0, 0, 0]])                # on-site defect, full support -> exact
    rng = np.random.default_rng(1)
    v = rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))
    V_loc = 0.5 * (v + v.conj().T)                 # (nw, nw) Hermitian, at R0=0
    assert np.allclose(V_loc, V_loc.conj().T, atol=1e-12), "V_loc not Hermitian"

    kc = lt.mp_grid(Nc, Nc, 1)                      # coarse grid = internal = output (golden)
    _, Ec, Uc = Hwr_to_Hwk(Hwr, Rw, kc, ndegen=ndeg)   # (Nk,nw),(Nk,nw,nw)
    Nk = len(kc)

    # dense reference M_bk[n,k,n',k'] = (1/Nk) phi_nk^dag V_loc phi_n'k'
    ph = lt._phase(kc, R_local)                     # (Nk, 1)
    phi = np.einsum("kL,kwn->knLw", ph, Uc, optimize=True).reshape(Nk, nw, len(R_local) * nw)
    M_bk = np.einsum("kna,ab,KMb->nkMK", np.conj(phi), V_loc, phi, optimize=True) / Nk

    # dense on-shell Gamma via compute_T (eigs = Ec.T is (nw,Nk))
    eigs = Ec.T
    Gd = np.zeros((nw, Nk))
    for ik in range(Nk):
        for n in range(nw):
            T = compute_T(np.array([Ec[ik, n]]), eigs, M_bk, eta)   # (1,nw,Nk,nw,Nk)
            Gd[n, ik] = -2.0 * T[0, n, ik, n, ik].imag
    Gd_perdef = Gd * Nk                              # dense is ~1/Nk; per-defect = *Nk

    # local scheme (same grids)
    Gl = lt.scattering_rate(Hwr, Rw, ndeg, V_loc, R_local, kc, eta, k_int=kc, perdef=True)

    err = np.max(np.abs(Gl - Gd_perdef))
    rel = err / max(1e-30, np.max(np.abs(Gd_perdef)))
    print(f"[golden] max|Gamma_local - Gamma_dense*Nk| = {err:.3e}  rel = {rel:.3e}")
    print(f"[golden] Gamma range (local): {Gl.min():.3e} .. {Gl.max():.3e}")
    print(f"[positivity] min Gamma = {Gl.min():.3e}  -> {'OK' if Gl.min() >= -1e-8 else 'FAIL'}")
    ok = rel < 1e-8 and Gl.min() >= -1e-8
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
