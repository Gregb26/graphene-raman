"""local_green_batch must equal local_green to roundoff (exact restructuring), on real Wannier data
(wannier/27x27, tracked); then scattering_rate_fast (batched, energy grid) vs scattering_rate (per state).
Run from the repo root. Prints PASS/FAIL and exits 0/1; wrapped by tests/test_scripts_tmatrix.py."""
import sys, time
import numpy as np
from graphene_raman.io.wannier_io import read_w90_tb
from graphene_raman.defects.many_body import local_tmatrix as lt


def main(tb_path="wannier/27x27/wannier_tb.dat"):
    Hwr, Rw, nd, _, _ = read_w90_tb(tb_path)
    k_int = lt.mp_grid(90, 90, 1)
    Hwk, _, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=nd)
    R_all = np.array([(i, j, 0) for i in range(-3, 4) for j in range(-3, 4)])
    ok = True
    for rc in (0, 1, 2, 3):
        Rloc = R_all[np.linalg.norm(R_all, axis=1) <= rc + 1e-9]
        egrid = np.linspace(-6.0, -2.0, 7)
        t0 = time.time(); ref = np.array([lt.local_green(Hwk, k_int, Rloc, e, 0.02) for e in egrid]); t1 = time.time()
        new = lt.local_green_batch(Hwk, k_int, Rloc, egrid, 0.02, k_chunk=3000, e_chunk=4); t2 = time.time()
        err = np.abs(new - ref).max() / np.abs(ref).max()
        ok &= err < 1e-12
        print(f"Rcut={rc} nL={len(Rloc)} nD={len(lt._diff_table(Rloc)[0])}: max rel |batch - ref| = {err:.2e}  (old {t1-t0:.2f}s, new {t2-t1:.2f}s)   {'OK' if err < 1e-12 else 'FAIL'}")
    # end-to-end: scattering_rate_fast (batched) vs scattering_rate (per-state exact) on a tiny case
    V = np.zeros((5, 5)); V[3, 3] = 7.0 / 27.211386; V[4, 4] = 0.05 / 27.211386; V[3, 4] = V[4, 3] = -0.5 / 27.211386
    k_out = lt.mp_grid(12, 12, 1)
    g_ref = lt.scattering_rate(Hwr, Rw, nd, V, [(0, 0, 0)], k_out, 0.05, k_int=k_int)
    g_new = lt.scattering_rate_fast(Hwr, Rw, nd, V, [(0, 0, 0)], k_out, 0.05, k_int=k_int, ne_per_eta=8)
    m = np.isfinite(g_new)
    rel = np.abs(g_new[m] - g_ref[m]).max() / np.abs(g_ref[m]).max()
    ok &= rel < 1e-2   # energy-grid interpolation (ne_per_eta = 8), not roundoff: 5.4e-3 measured (C8 of NOTES_TGAMMA)
    print(f"fast(batched, ne_per_eta=8) vs exact per-state: max rel diff = {rel:.2e} (energy-grid interpolation, not roundoff; seuil 1e-2)   {'OK' if rel < 1e-2 else 'FAIL'}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
