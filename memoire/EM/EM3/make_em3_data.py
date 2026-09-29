"""
EM3: data of the §2.5 figure and table, from the package functions and the EM2 / M4 outputs (no number copied by hand
except the hermiticity budget, EM.md section 5):
  em2_postw90_sigma_ti.npz : postw90 kubo of the run k1201_ti (transl_inv = true, reference curve, EM.md decision 7),
                             read by read_postw90 of EM2/em2_B_compare.py, same format as EM2/em2_postw90_sigma.npz;
  em3_ring_2p33.npz        : in-plane |hbar v_cv| on the 2.33 eV ring (ring, ntheta = 720, mu = E_D), three variants;
  em_table.md, em_table.csv: numbers of the section (hbar v_F, |hbar v_cv| against the DFT, ring averages, sigma at the
                             lasers with postw90, van Hove, hbar omega_froz, error budgets).
Stops (AssertionError) if a control fails. Run: .venv/bin/python memoire/EM/EM3/make_em3_data.py
"""

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EM2, M4 = REPO / "memoire" / "EM" / "EM2", REPO / "memoire" / "EM" / "M4_sigma"
os.environ.setdefault("GRAPHENE_RAMAN", str(REPO))           # read by em2_B_compare at import (PROJECTS too: its
os.environ.setdefault("PROJECTS", str(REPO.parent))          # os.environ.get default is evaluated eagerly)
sys.path.insert(0, str(EM2))
from em2_B_compare import read_postw90, at, G_S, SIGMA0, C_BOHR, BOHR_W90   # noqa: E402
from electron_defect_interaction.electron_photon import (centres_only, compute_velocity, make_grid_tb,  # noqa: E402
                                                         make_wannier_tb, ring)

E_D = -4.238895                          # eV (EM.md section 2)
HW_RING, NTHETA = 2.33, 720              # eV, angles (0.5 degree; the 48 DFT angles are every 15th)
LASERS = {1.96: 633, 2.33: 532, 2.54: 488}
HBAR = 6.582119569e-16                   # eV s
HERM_BUDGET = 1.3e-3                     # max |delta hbar v_cv| / hbar v_F dropped by hermitize, 2.33 eV ring (EM.md section 5)
VARIANTS = ("full", "centres_only", "no_berry")
log = []


def say(s=""):
    print(s); log.append(s)


# ---------------- postw90, transl_inv
win = (EM2 / "postw90" / "k1201_ti" / "wannier.win").read_text()
assert "transl_inv = true" in win and "berry_kmesh = 1201 1201 1" in win and "kubo_smr_fixed_en_width = 0.0565685425" in win
hw_p, sig_p = read_postw90("k1201_ti")
np.savez(HERE / "em2_postw90_sigma_ti.npz", hw=hw_p, sigma=sig_p, berry_kmesh="k1201", eta=0.04, width_w90=0.0565685425,
         mu=E_D, g_s=G_S, sigma0_SI=SIGMA0, c_A=C_BOHR * BOHR_W90, source="postw90/k1201_ti/wannier-kubo_S_*.dat",
         transl_inv=True)

# ---------------- ring at 2.33 eV, three variants
tb = make_wannier_tb(REPO / "wannier" / "27x27" / "wannier_tb.dat")
K = make_grid_tb(tb, 3).K
models = {"full": (tb, "berry"), "centres_only": (centres_only(tb), "berry"), "no_berry": (tb, "no_berry")}


def hv_inplane(k):
    """{variant: in-plane |hbar v_cv| (nk,)}, c and v by energy around E_D (same pair at every k, asserted)."""
    out = {}
    for name, (model, mode) in models.items():
        _, eps, _, hv = compute_velocity(model, k, mode)
        n = np.sum(eps < E_D, axis=1); assert np.all(n == n[0])
        out[name] = np.sqrt(np.abs(hv[:, 0, n[0], n[0] - 1])**2 + np.abs(hv[:, 1, n[0], n[0] - 1])**2)
    return out


k_ring, theta, q, _ = ring(tb, K, HW_RING, mu=E_D, ntheta=NTHETA)
hv_ring = hv_inplane(k_ring)
scal = json.loads((M4 / "em_scalars.json").read_text())
hv_F = scal["hv_F_eV_A"]["inter"]
np.savez(HERE / "em3_ring_2p33.npz", theta=theta, q=q, **{f"hv_{v}": hv_ring[v] for v in VARIANTS}, hv_F=hv_F,
         hw=HW_RING, mu=E_D, ntheta=NTHETA)

# ---------------- controls
A = np.load(EM2 / "em2_A.npz")
m = (A["sets"] == "ring") & np.isclose(A["hw"], HW_RING)
B = 2 * np.pi * np.linalg.inv(tb.lattice).T
same_k = hv_inplane(A["k_red"][m] @ B.T)
d1 = max(np.abs(same_k[v] - A[f"hv_{v}"][m]).max() for v in VARIANTS)
d_theta = np.abs(theta[::15] - A["theta"][m]).max()
d2 = max(np.abs(hv_ring[v][::15] - A[f"hv_{v}"][m]).max() for v in VARIANTS)
say(f"control 1: at the {m.sum()} DFT k, |hv_cv| recomputed - em2_A {d1:.1e}; ring grid [::15]: theta {d_theta:.1e} rad, "
    f"|hv_cv| {d2:.1e} eV A")
assert m.sum() == 48 and d1 < 1e-10 and d_theta < 1e-12 and d2 < 1e-10

SIG = {v: np.load(M4 / f"em_sigma_{v}_N1200_eta0.04.npz") for v in VARIANTS}
hw = SIG["full"]["hw"]; assert np.abs(hw - hw_p).max() < 1e-6
iL = {e: at(hw, e) for e in LASERS}
sxx = np.array([SIG["full"]["sigma"][i, 0, 0] for i in iL.values()])
say(f"control 2: M4 full sigma_xx/sigma_0 at the lasers {np.round(sxx, 5)}")
assert np.abs(sxx - [1.2637, 1.4056, 1.5115]).max() < 1e-4
d3 = sig_p[iL[2.33], 0, 0] - SIG["full"]["sigma"][iL[2.33], 0, 0]
say(f"control 3: postw90 ti - M4 at 2.33 eV (xx) {d3:+.2e}")
assert abs(d3 + 4.4e-4) < 5e-6
s02 = sig_p[at(hw, 0.2), 0, 0]
say(f"control 4: postw90 ti sigma_xx/sigma_0 at 0.2 eV {s02:.5f} (g_s = {G_S} applied)")
assert abs(s02 - 1.002) < 1e-3

# ---------------- numbers of the table
rows = []   # (block, quantity, context, variant, value, unit)
hvF_dft = float(A["hvF_dft"])
rows += [("vF", "hbar v_F", "Wannier, q = 1e-3 1/A", "full", hv_F, "eV A"),
         ("vF", "hbar v_F", "DFT, slopes on q = 0.005 1/A", "dft", hvF_dft, "eV A"),
         ("vF", "v_F", "Wannier", "full", hv_F * 1e-10 / HBAR, "m/s"),
         ("vF", "v_F", "DFT", "dft", hvF_dft * 1e-10 / HBAR, "m/s")]
for e in LASERS:
    mm = (A["sets"] == "ring") & np.isclose(A["hw"], e)
    for v in VARIANTS:
        r = 100 * (A[f"hv_{v}"][mm] / A["hv_qe"][mm] - 1)
        rows += [("dft", f"|hv_cv|/DFT - 1 {s}", f"{e} eV ({mm.sum()} k)", v, f(r), "%")
                 for s, f in (("min", np.min), ("mean", np.mean), ("max", np.max))]
stats = json.loads((M4 / "em_ring_stats.json").read_text())
for e in LASERS:
    st = stats[str(e)]
    for v in VARIANTS:
        rows += [("avg", f"<|e_{a}.hv_cv|^2>", f"{e} eV", v, st[v]["avg"][i], "(hbar v_F)^2/2") for i, a in enumerate("xy")]
    rows += [("avg", f"<|e_{a}.hv_cv|^2> centres/full", f"{e} eV", "centres_only",
              st["centres_only"]["avg"][i] / st["full"]["avg"][i], "") for i, a in enumerate("xy")]
for e, i in iL.items():
    for v in VARIANTS:
        rows.append(("sigma", "sigma_xx/sigma_0", f"{e} eV", v, SIG[v]["sigma"][i, 0, 0], ""))
    rows.append(("sigma", "sigma_yy/sigma_0", f"{e} eV", "no_berry", SIG["no_berry"]["sigma"][i, 1, 1], ""))
    rows.append(("sigma", "sigma_xx/sigma_0", f"{e} eV", "postw90_ti", sig_p[i, 0, 0], ""))
    rows.append(("sigma", "postw90 ti - full (xx)", f"{e} eV", "postw90_ti", sig_p[i, 0, 0] - SIG["full"]["sigma"][i, 0, 0], ""))
ivh = int(np.argmax(SIG["full"]["sigma"][:, 0, 0]))
rows += [("misc", "van Hove peak (position)", "full, sigma_xx", "full", hw[ivh], "eV"),
         ("misc", "hbar omega_froz", f"N = {scal['N_froz']}", "full", scal["hw_froz_eV"], "eV")]

# budgets at 2.33 eV (sigma_xx/sigma_0 of the full mode unless stated)
pil = np.load(M4 / "pilote" / "em_sigma_full_N1200.npz")
ie = {float(x): j for j, x in enumerate(pil["eta"])}
s02_, s04_ = pil["sigma"][ie[0.02], iL[2.33], 0, 0], pil["sigma"][ie[0.04], iL[2.33], 0, 0]
eta_bias = (s04_ - s02_) * 4 / 3                     # sigma(0.04) - sigma(0), quadratic in eta (Richardson on 0.02, 0.04)
ours = np.load(EM2 / "em2_B_ours.npz")
prefactor = ours["shift_delta"][iL[2.33], 0, 0] - SIG["full"]["sigma"][iL[2.33], 0, 0]
ws = sig_p[iL[2.33], 0, 0] - read_postw90("k1201_ti_nows")[1][iL[2.33], 0, 0]
dft_max = np.abs(A["hv_full"][A["sets"] == "ring"] / A["hv_qe"][A["sets"] == "ring"] - 1).max()
rows += [("budget", "hermiticity of r (dropped by hermitize), max |d hv_cv|/hbar v_F", "2.33 eV ring", "full", HERM_BUDGET, ""),
         ("budget", "|hv_cv| full vs DFT, max relative", "72 ring k", "full", dft_max, ""),
         ("budget", "eta bias sigma(0.04) - sigma(0)", "2.33 eV", "full", eta_bias, ""),
         ("budget", "prefactor 1/Delta - 1/hw", "2.33 eV", "full", prefactor, ""),
         ("budget", "use_ws_distance (postw90 ti - ti_nows)", "2.33 eV", "postw90_ti", ws, "")]

with open(HERE / "em_table.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["block", "quantity", "context", "variant", "value", "unit"])
    w.writerows([r[:4] + (f"{r[4]:.6g}",) + r[5:] for r in rows])

# ---------------- markdown
val = {(r[1], r[2], r[3]): r[4] for r in rows}
V = {"full": "complet", "centres_only": "centres seuls", "no_berry": "sans Berry"}
md = ["# Chiffres du §2.5 (EM3)", "",
      f"Généré par `make_em3_data.py` (ne pas éditer à la main). Wannierisation 27×27, μ = E_D = {E_D} eV ; σ : N = 1200, η = 0.04 eV.",
      "", "## Vitesse de Fermi", "", "| | ħv_F (eV·Å) | v_F (10⁵ m/s) |", "|---|---|---|",
      f"| Wannier (complet, q = 10⁻³ Å⁻¹) | {hv_F:.4f} | {hv_F * 1e-15 / HBAR:.3f} |",
      f"| DFT directe (pentes, q = 0.005 Å⁻¹) | {hvF_dft:.5f} | {hvF_dft * 1e-15 / HBAR:.3f} |",
      "", "## |ħv_cv| dans le plan contre la DFT directe (écart relatif k par k, %, min / moy. / max)", "",
      "| Anneau | k | complet | centres seuls | sans Berry |", "|---|---|---|---|---|"]
for e in LASERS:
    n = int(((A["sets"] == "ring") & np.isclose(A["hw"], e)).sum())
    md.append(f"| {e} eV ({LASERS[e]} nm) | {n} | " + " | ".join(
        " / ".join(f"{val[(f'|hv_cv|/DFT - 1 {s}', f'{e} eV ({n} k)', v)]:+.3f}" for s in ("min", "mean", "max")) for v in VARIANTS) + " |")
md += ["", "## Moyennes sur les anneaux (F15), unités (ħv_F)²/2", "",
       "| Anneau | ⟨\\|e_x·ħv_cv\\|²⟩ complet / centres / sans Berry | ⟨\\|e_y·ħv_cv\\|²⟩ complet / centres / sans Berry | centres/complet x ; y |",
       "|---|---|---|---|"]
for e in LASERS:
    cells = [" / ".join(f"{val[(f'<|e_{a}.hv_cv|^2>', f'{e} eV', v)]:.4f}" for v in VARIANTS) for a in "xy"]
    ratio = " ; ".join(f"{val[(f'<|e_{a}.hv_cv|^2> centres/full', f'{e} eV', 'centres_only')]:.4f}" for a in "xy")
    md.append(f"| {e} eV | {cells[0]} | {cells[1]} | {ratio} |")
md += ["", "## σ(ω)/σ₀ aux énergies laser", "",
       "| λ (nm) | ħω (eV) | complet xx | centres xx | sans Berry xx | sans Berry yy | postw90 (`transl_inv`) xx | postw90 − complet |",
       "|---|---|---|---|---|---|---|---|"]
for e in LASERS:
    c = f"{e} eV"
    md.append(f"| {LASERS[e]} | {e} | " + " | ".join(f"{val[('sigma_xx/sigma_0', c, v)]:.4f}" for v in VARIANTS)
              + f" | {val[('sigma_yy/sigma_0', c, 'no_berry')]:.4f} | {val[('sigma_xx/sigma_0', c, 'postw90_ti')]:.4f}"
              f" | {val[('postw90 ti - full (xx)', c, 'postw90_ti')]:+.1e} |")
md += ["", f"Pic de van Hove (position seulement) : {hw[ivh]:.2f} eV ; ħω_froz = {scal['hw_froz_eV']:.4f} eV (N = {scal['N_froz']}).",
       "", "## Budgets d'erreur (à 2.33 eV)", "", "| Source | Grandeur | Valeur |", "|---|---|---|",
       f"| hermiticité de r (partie jetée par `hermitize`) | max \\|δħv_cv\\|/ħv_F sur l'anneau | {HERM_BUDGET:.1e} |",
       f"| DFT directe (contrôle du budget précédent) | max \\|écart relatif\\| de \\|ħv_cv\\| complet, 72 k | {dft_max:.1e} |",
       f"| élargissement η = 0.04 eV | σ_xx(η) − σ_xx(0) | {eta_bias:+.1e} |",
       f"| préfacteur 1/(ε_c − ε_v) de postw90 contre 1/ħω | Δσ_xx | {prefactor:+.1e} |",
       f"| `use_ws_distance` (non appliqué, décision 4) | Δσ_xx (postw90) | {ws:+.1e} |"]
commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
dirty = bool(subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "src"], capture_output=True, text=True).stdout.strip())
md += ["", f"Commit {commit}{' (src modifié)' if dirty else ''}. Sources : `M4_sigma/` (σ, anneaux, scalaires, pilote), "
       "`EM2/em2_A.npz` (DFT), `EM2/postw90/k1201_ti*`, `EM2/em2_B_ours.npz` ; hermiticité : EM.md §5 (F10)."]
(HERE / "em_table.md").write_text("\n".join(md) + "\n")
say("\n".join(md))
(HERE / "make_em3_data.log").write_text("\n".join(log) + "\n")
