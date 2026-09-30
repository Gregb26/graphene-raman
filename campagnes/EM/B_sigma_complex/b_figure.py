"""
Figure of perspective B (exploratory, not a thesis figure): b_sigma_complex.{pdf,png} next to this script.
  (a) undoped, Re and Im sigma_xx/sigma_0 at 300 K, Im at T = 0: Kramers-Kronig signature around the van Hove peak;
  (b) sigma(mu) - sigma(E_D) at 300 K for n and p doping (mu = E_D +/- 0.3 eV) against the Dirac cone at 300 K with the
      same kernel, (i/pi) int_0^inf dD [f(-D/2) - f(D/2)] [K(hw - D) + K(hw + D)], by quadrature (the lattice
      background cancels in the difference).
Run: .venv/bin/python campagnes/EM/B_sigma_complex/b_figure.py (after b_prod.py)
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad
from electron_defect_interaction.electron_photon import gaussian_complex

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
plt.style.use(ROOT / "figures" / "memoire.mplstyle")
sys.path.insert(0, str(ROOT / "scripts"))
from _palette import NAVY, ORANGE, SKY, REF, INK, MUTED   # noqa: E402

Z = np.load(HERE / "b_sigma_complex_N1200_eta0.04.npz")
hw, eta, E_D = Z["hw"], float(Z["eta"]), float(Z["E_D"])
S = {str(c): Z["sigma"][i][:, 0, 0] for i, c in enumerate(Z["cases"])}
MU = {str(c): float(m) - E_D for c, m in zip(Z["cases"], Z["mu"])}
KT = 0.025


def cone(h, mu, kT):
    """Dirac cone, gaussian_complex, occupations f(-D/2) - f(D/2) at chemical potential mu, temperature kT."""
    f = lambda e: 0.5 * (1 - np.tanh((e - mu) / (2 * kT)))
    g = lambda D: (f(-D / 2) - f(D / 2)) * (gaussian_complex(h - D, eta) + gaussian_complex(h + D, eta))
    parts = [(0, h), (h, np.inf)]
    return 1j / np.pi * sum(quad(lambda D: g(D).real, *p, limit=400)[0] + 1j * quad(lambda D: g(D).imag, *p, limit=400)[0]
                            for p in parts)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.3))

# (a) undoped
a1.plot(hw, S["undoped_300K"].real, color=NAVY, lw=1.4, label=r"Re, 300 K")
a1.plot(hw, S["undoped_300K"].imag, color=NAVY, lw=1.2, ls="--", label=r"Im, 300 K")
a1.plot(hw, S["undoped_T0"].imag, color=SKY, lw=1.0, ls=":", label=r"Im, $T = 0$")
a1.axhline(0, color=MUTED, lw=0.6, zorder=0); a1.axhline(1, color=MUTED, lw=0.6, ls=":", zorder=0)
a1.set_xlim(0, 6); a1.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)"); a1.set_ylabel(r"$\sigma_{xx}(\omega)/\sigma_0$")
a1.legend(fontsize=7, frameon=False, loc="upper left")
a1.set_title("(a) non dopé", loc="left", fontsize=10)

# (b) doping difference against the cone
m = hw <= 1.6
h = hw[m]
c = np.array([cone(x, 0.3, KT) - cone(x, 0.0, KT) for x in h])
for case, col, lab in (("n_300K", NAVY, "dopage n"), ("p_300K", ORANGE, "dopage p")):
    d = (S[case] - S["undoped_300K"])[m]
    a2.plot(h, d.real, color=col, lw=1.3, label=rf"Re, {lab}")
    a2.plot(h, d.imag, color=col, lw=1.1, ls="--", label=rf"Im, {lab}")
a2.plot(h, c.real, color=REF, lw=0.8, ls=(0, (1, 1.5)), label="cône de Dirac")
a2.plot(h, c.imag, color=REF, lw=0.8, ls=(0, (1, 1.5)))
a2.axhline(0, color=MUTED, lw=0.6, zorder=0)
a2.set_xlim(0.2, 1.6); a2.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)")
a2.set_ylabel(r"$\sigma_{xx}(\mu) - \sigma_{xx}(E_D)$ (unités de $\sigma_0$)")
a2.legend(fontsize=7, frameon=False, loc="lower right")
a2.set_title(r"(b) $\mu = E_D \pm 0.3$ eV, 300 K", loc="left", fontsize=10)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"b_sigma_complex.{ext}")
dn, dp = (S["n_300K"] - S["undoped_300K"])[m], (S["p_300K"] - S["undoped_300K"])[m]
print(f"max |(n) - cone| {np.abs(dn - c).max():.1e}, max |(p) - cone| {np.abs(dp - c).max():.1e}, max |n - p| {np.abs(dn - dp).max():.1e}")
print("écrit", HERE / "b_sigma_complex.pdf/.png")
