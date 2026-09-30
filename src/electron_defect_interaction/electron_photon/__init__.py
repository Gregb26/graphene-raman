"""
Electron-photon coupling (EM series): Wannier-interpolated velocity operator and optical response.

From a Wannier tight-binding model (H(R), r(R), layout of Wannier90's `_tb.dat`), computes

    hbar v(k) = V(k)^dagger [ dH(k)/dk + i [H(k), A(k)] ] V(k),

with V(k) diagonalizing H(k) and A(k) the Berry connection in the Wannier gauge. Modules, in order
of dependence:
    tb_model           WannierTB, analytic graphene model (M0), make_wannier_tb from a _tb.dat,
                       centres_only (M1), pz_block (M3)
    kgrid              reciprocal lattice, k grids and paths, GridTB
    velocity_operator  fourier (H, dH, A) -> hermitize (A) -> velocity, chained by compute_velocity
    ring               resonant k points around K (eps_c - eps_v = hbar omega), fermi_velocity,
                       ring_stats and ring_kpoints_crystal (M3), map_around_K (M4)
    kubo               sigma(omega)/sigma_0, block by block, driver `sigma_on_grid`, doping and temperature
                       (kT, perspective A) with the reference kubo_doped_finite_T_analytical, complex sigma
                       with a complex kernel (gaussian_complex, lorentzian_complex; perspective B)
    diagnostics        reports on a model: hermiticity_report, symmetry_report (M1),
                       frozen_window_limit (M4)
Plan, derivations and reference values: memoire/EM/EM.md.

Conventions:
    - Units: Angstrom, eV, 1/Angstrom; hbar v in eV*Angstrom.
    - Vectors are always 3D (z = 0), as in `_tb.dat`.
    - X_mn(R) = <0m| X |Rn> and X(k) = sum_R e^{+ik.R} X(R) / ndegen(R) (lattice gauge, no tau in
      the phases), computed by the single routine `fourier` for H, dH/dk and A.
    - `lattice` and `B` hold their basis vectors in COLUMNS; lists of vectors (R, k) are ROWS,
      shape (N, 3), so that X_cart = X_red @ lattice.T.
    - Occupations are decided by energy (eps < mu), never by band index.
"""

from .tb_model import WannierTB, make_graphene_tb, make_wannier_tb, centres_only, pz_block
from .kgrid import GridTB, reciprocal, k_grid, make_grid_tb, kpath
from .velocity_operator import dagger, hermitize, fourier, velocity, compute_velocity
from .ring import ring, fermi_velocity, ring_stats, ring_kpoints_crystal, map_around_K
from .kubo import (gaussian_eta, lorentzian_complex, gaussian_complex, kubo_accumulate, kubo_normalize,
                   sigma_on_grid, kubo_doped_finite_T_analytical, kubo_doped_complex_analytical)
from .diagnostics import hermiticity_report, symmetry_report, frozen_window_limit

__all__ = [
    "WannierTB", "make_graphene_tb", "make_wannier_tb", "centres_only", "pz_block",
    "GridTB", "reciprocal", "k_grid", "make_grid_tb", "kpath",
    "dagger", "hermitize", "fourier", "velocity", "compute_velocity",
    "ring", "fermi_velocity", "ring_stats", "ring_kpoints_crystal", "map_around_K",
    "gaussian_eta", "lorentzian_complex", "gaussian_complex", "kubo_accumulate", "kubo_normalize",
    "sigma_on_grid", "kubo_doped_finite_T_analytical", "kubo_doped_complex_analytical",
    "hermiticity_report", "symmetry_report", "frozen_window_limit",
]
