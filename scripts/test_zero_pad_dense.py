"""
test_zero_pad_dense.py
    Validation of the exact zero-padding densification of M^L (defects/local_G.py).

    Test 1 (standalone, no data): zero_pad_potential on a random real potential
        - coincident reciprocal coefficients are preserved exactly: Ved_G_pad[::p,::p,:] == fftn(Ved)/N
        - the padded potential is real: max|ifftn(Ved_G_pad).imag| ~ machine eps
        - real-space round trip: ifftn(Ved_G_pad).real[:Nx,:Ny,:]/p^2 == Ved, padding region == 0

    Test 2 (non-regression, needs the 5x5 .save data): feeding the ORIGINAL coarse k-points back
        through the dense code path (Ndiag*p, ngfft*p, Ved_G_pad) must reproduce compute_ML_G to
        machine precision -- the coarse grid is the subset of the dense grid where delta_k_sc lands
        exactly on the preserved coefficients.

    Prints PASS/FAIL and exits 0/1 (no pytest), like the other scripts/ validation tests.
"""
import sys
import numpy as np

import _paths
from electron_defect_interaction.defects.local_G import (
    zero_pad_potential, prep_reciprocal_inputs, compute_ML_G, compute_ML_G_dense,
)

# ---- Test 2 data (see scripts/_paths.py) ----
UC    = _paths.uc("5x5")
SC_P  = _paths.sc_p("5x5")
SC_D  = _paths.sc_d("5x5")
POT_P = _paths.pot_p("5x5")
POT_D = _paths.pot_d("5x5")
BANDS = range(4)          # restrict to a few bands to keep the O(nk^2) test fast


def test_zero_pad(p=3):
    rng = np.random.default_rng(0)
    Nx, Ny, Nz = 8, 10, 3
    Ved = rng.standard_normal((Nx, Ny, Nz))          # real defect potential

    Ved_G_pad = zero_pad_potential(Ved, p)
    No = Ved.size                                     # original Nx*Ny*Nz

    ok_shape = Ved_G_pad.shape == (p * Nx, p * Ny, Nz)

    # (a) coincident reciprocal coefficients preserved exactly (drop-in for prep["Ved_G"])
    Ved_G = np.fft.fftn(Ved) / No
    err_coeff = np.max(np.abs(Ved_G_pad[::p, ::p, :] - Ved_G))

    # (b)/(c) real-space round trip, reconstructed with the SAME convention as prep (Ved = ifftn*No)
    rec = np.fft.ifftn(Ved_G_pad)
    err_imag = np.max(np.abs(rec.imag))
    Ved_rec = rec.real * No                           # == real-space zero-padded potential
    expected = np.zeros_like(Ved_rec)
    expected[:Nx, :Ny, :] = Ved
    err_real = np.max(np.abs(Ved_rec[:Nx, :Ny, :] - Ved))
    err_pad  = np.max(np.abs(Ved_rec - expected))      # also catches the zero padding region

    ok = ok_shape and err_coeff < 1e-13 and err_imag < 1e-12 and err_real < 1e-12 and err_pad < 1e-12
    print("=== Test 1: zero_pad_potential ===")
    print(f"  shape {Ved_G_pad.shape} (expected {(p*Nx, p*Ny, Nz)})  ok={ok_shape}")
    print(f"  max|coincident coeff - Ved_G|        = {err_coeff:.2e}")
    print(f"  max|imag(ifftn(Ved_G_pad))|          = {err_imag:.2e}")
    print(f"  max|ifftn*No [:Nx,:Ny] - Ved|        = {err_real:.2e}")
    print(f"  max|ifftn*No - zero_padded(Ved)|     = {err_pad:.2e}")
    print(f"  -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_non_regression(p=3):
    from electron_defect_interaction.io import qe_io

    print("=== Test 2: dense path reproduces compute_ML_G at the coarse k-points ===")
    try:
        prep = prep_reciprocal_inputs(UC, SC_P, POT_P, POT_D, bands=BANDS)
    except (FileNotFoundError, OSError) as e:
        print(f"  SKIP (data not available: {e})\n")
        return None

    M_orig = compute_ML_G(prep, block_size=512)

    # raw primitive quantities of the SAME (coarse) grid: it is the subset of the dense grid whose
    # delta_k_sc * p lands exactly on the preserved coefficients of Ved_G_pad.
    C_nkg, nG = qe_io.get_C_nk(UC)
    C_dense = C_nkg[list(BANDS), ...]
    G_dense = qe_io.get_G_red(UC)
    k_dense = qe_io.get_k_red(UC)

    M_dense = compute_ML_G_dense(prep, p, k_dense, C_dense, G_dense, nG, block_size=512)

    err = np.max(np.abs(M_dense - M_orig))
    scale = np.max(np.abs(M_orig))
    ok = err < 1e-10 * max(scale, 1.0)
    print(f"  p={p}  shapes {M_orig.shape} vs {M_dense.shape}")
    print(f"  max|M_dense - M_orig| = {err:.2e}   (|M|max = {scale:.2e})")
    print(f"  -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    r1 = test_zero_pad(p=3)
    r2 = test_non_regression(p=3)

    results = [r for r in (r1, r2) if r is not None]
    allok = all(results)
    print("RESULT:", "PASS" if allok else "FAIL", "" if r2 is not None else "(Test 2 skipped)")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
