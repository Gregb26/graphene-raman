#!/usr/bin/env python
"""
nkint_check_post.py -- P13: convergence of the INTERNAL k-grid N_k^int of g0 (reference size, R_cut 3, frozen grid/eta/window).
Reads results/M/resigma_<size>_rc3_nk<N>.npz (scripts/rcut_resigma.py --nk-int, job submit_nkint_check.sh) and reports, per nk_int:
  * STATE median of Gamma (= Gamma N_cells, intensive V_loc) over the on-shell states of the +-e_window (NOT the Lorentzian curve),
  * median Re Sigma, E_res = argmax Gamma over states with |eps - E_D| <= 1.5 eV, Gamma_T at E_D (mean of the states with eps = E_D, i.e. K),
  * the same restricted to |eps - E_D| <= 0.3 eV,
with relative deviations to the densest nk_int. Writes results/M/nkint_check_<size>.csv and prints markdown tables.
Usage: python scripts/nkint_check_post.py --size 9x9 --nk 150,300,450,600 [--pattern results/M/resigma_{S}_rc3_nk{N}.npz]
"""
import argparse, csv, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--size", default="9x9"); ap.add_argument("--nk", default="150,300,450,600")
ap.add_argument("--pattern", default="results/M/resigma_{S}_rc3_nk{N}.npz"); ap.add_argument("--zoom", type=float, default=0.3); a = ap.parse_args()
nks = [int(x) for x in a.nk.split(",")]; rows = []
for N in nks:
    f = a.pattern.format(S=a.size, N=N); z = np.load(f); S = z["Sigma_rc3"]; E = z["E_out"]; ED = float(z["E_D"]); ew = float(z["e_window"])
    assert int(z["nk_int"]) == N, (f, z["nk_int"]); m = np.isfinite(S.real); G = -2 * S.imag; R = S.real; d = np.abs(E - ED)
    def stats(msk):
        mm = msk & (d <= 1.5); e_res = float(E[mm][np.argmax(G[mm])] - ED) if mm.any() else np.nan
        return dict(n=int(msk.sum()), medG=float(np.median(G[msk])) * 1e3, medRe=float(np.median(R[msk])) * 1e3, E_res=e_res)
    full = stats(m); zoom = stats(m & (d <= a.zoom)); atED = m & (d < 1e-6)
    rows.append(dict(nk_int=N, nk_tot=N * N, file=f, n_states=full["n"], medG=full["medG"], medRe=full["medRe"], E_res=full["E_res"],
                     G_ED=float(G[atED].mean()) * 1e3, n_ED=int(atED.sum()), z_n=zoom["n"], z_medG=zoom["medG"], z_medRe=zoom["medRe"], z_E_res=zoom["E_res"]))
ref = rows[-1]
def rel(r, k): return 100.0 * (r[k] - ref[k]) / ref[k]
print(f"\n### N_k^int, {a.size}, R_cut 3 — états de la fenêtre ±{ew:.0f} eV (référence : nk_int = {ref['nk_int']})\n")
print("| nk_int | N_k^int | n états | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) | Γ_T(E_D) (meV, états à K) | écart |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['nk_int']} | {r['nk_tot']} | {r['n_states']} | {r['medG']:.2f} | {rel(r,'medG'):+.2f} % | {r['medRe']:.2f} | {rel(r,'medRe'):+.2f} % | {r['E_res']:+.3f} | {r['G_ED']:.2f} ({r['n_ED']} états) | {rel(r,'G_ED'):+.2f} % |")
print(f"\n### idem, états à |ε − E_D| ≤ {a.zoom} eV\n")
print("| nk_int | n états | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) |")
print("|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['nk_int']} | {r['z_n']} | {r['z_medG']:.2f} | {rel(r,'z_medG'):+.2f} % | {r['z_medRe']:.2f} | {rel(r,'z_medRe'):+.2f} % | {r['z_E_res']:+.3f} |")
out = f"results/M/nkint_check_{a.size}.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"\nécrit {out}")
