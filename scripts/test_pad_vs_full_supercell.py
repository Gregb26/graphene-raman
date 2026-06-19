"""
test_pad_vs_full_supercell.py
    PHYSICAL validation of the zero-padding densification: does a SMALL supercell defect potential,
    zero-padded by p, reproduce the M^L computed from the p-times LARGER supercell?

    This is NOT the machine-precision non-regression check (scripts/test_zero_pad_dense.py); here the
    two potentials are genuinely different physical objects and the residual measures the finite-size
    (defect-image) error of the small supercell. It is exact only in the limit where the small
    supercell already contains the full defect potential (compact support).

    Both matrices are evaluated on the SAME dense primitive k-grid with the SAME dense unit-cell
    wavefunctions, so the only difference is the source of the local potential:

        M_full = compute_ML_G(prep_large)                          # large supercell, native
        M_pad  = compute_ML_G_dense(prep_small, p, k_dense, ...)   # small supercell, zero-padded

    Only the LOCAL part M^L is compared (M^NL depends on the relaxed atomic positions, which differ
    between the two supercells and are unrelated to the padding mechanism).

    Example (6x6 -> 12x12, p=2; uc on the 12x12 grid):
        .venv/bin/python scripts/test_pad_vs_full_supercell.py \\
            --uc-dense   data/graphene/unit_cell/qe/defect_unit_cell_12x12.save \\
            --sc-small-p data/graphene/supercell/qe/defect_6x6_p.save \\
            --pot-small-p data/graphene/supercell/qe/Vks_6x6_p \\
            --pot-small-d data/graphene/supercell/qe/Vks_6x6_d \\
            --sc-large-p data/graphene/supercell/qe/defect_12x12_p.save \\
            --pot-large-p data/graphene/supercell/qe/Vks_12x12_p \\
            --pot-large-d data/graphene/supercell/qe/Vks_12x12_d \\
            --p 2 --bands 0,1,2,3
"""
import sys
import argparse
import numpy as np

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.defects.local_G import (
    prep_reciprocal_inputs, compute_ML_G, compute_ML_G_dense,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uc-dense", required=True,
                    help="unit-cell .save on the DENSE (p*N x p*N) k-grid; provides the wavefunctions "
                         "used by BOTH matrices")
    ap.add_argument("--sc-small-p", required=True, help="small pristine supercell .save (for A_sc geometry)")
    ap.add_argument("--pot-small-p", required=True, help="small pristine pp.x potential (Vks)")
    ap.add_argument("--pot-small-d", required=True, help="small defective pp.x potential (Vks)")
    ap.add_argument("--sc-large-p", required=True, help="large pristine supercell .save")
    ap.add_argument("--pot-large-p", required=True, help="large pristine pp.x potential (Vks)")
    ap.add_argument("--pot-large-d", required=True, help="large defective pp.x potential (Vks)")
    ap.add_argument("--p", type=int, required=True, help="densification factor (large = p * small)")
    ap.add_argument("--bands", default="0,1,2,3", help="comma-separated band indices, or 'all'")
    ap.add_argument("--tol", type=float, default=1e-3, help="relative tolerance for PASS (physical, not machine)")
    args = ap.parse_args()

    bands = None if args.bands.strip() == "all" else [int(b) for b in args.bands.split(",")]

    # --- ground truth: M^L from the LARGE supercell on the dense uc grid ---
    prep_large = prep_reciprocal_inputs(args.uc_dense, args.sc_large_p,
                                        args.pot_large_p, args.pot_large_d, bands=bands)
    # sanity: the large supercell must fold onto exactly the dense uc grid (Ndiag = p * Ndiag_small)
    print(f"Ndiag(large) = {prep_large['Ndiag']}  ngfft(large) = {prep_large['ngfft']}")
    M_full = compute_ML_G(prep_large)

    # --- padded: M^L from the SMALL supercell potential, zero-padded by p, on the SAME dense grid ---
    prep_small = prep_reciprocal_inputs(args.uc_dense, args.sc_small_p,
                                        args.pot_small_p, args.pot_small_d, bands=bands)
    print(f"Ndiag(small) = {prep_small['Ndiag']}  ngfft(small) = {prep_small['ngfft']}  "
          f"-> effective {prep_small['Ndiag'] * np.array([args.p, args.p, 1])}")
    assert np.array_equal(prep_small["Ndiag"] * np.array([args.p, args.p, 1]), prep_large["Ndiag"]), \
        "p * Ndiag(small) must equal Ndiag(large): supercell sizes are not commensurate with p"

    C_dense, nG_dense = qe_io.get_C_nk(args.uc_dense)
    if bands is not None:
        C_dense = C_dense[bands, ...]
    G_dense = qe_io.get_G_red(args.uc_dense)
    k_dense = qe_io.get_k_red(args.uc_dense)

    M_pad = compute_ML_G_dense(prep_small, args.p, k_dense, C_dense, G_dense, nG_dense)

    # --- compare (physical finite-size residual) ---
    diff = np.abs(M_pad - M_full)
    scale = np.max(np.abs(M_full))
    abs_err = float(np.max(diff))
    rel_err = abs_err / max(scale, 1e-30)
    rel_fro = float(np.linalg.norm(M_pad - M_full) / max(np.linalg.norm(M_full), 1e-30))
    ok = rel_err < args.tol

    print(f"\n=== M^L: small (padded x{args.p}) vs large supercell, dense grid ===")
    print(f"  shapes        : {M_full.shape}")
    print(f"  |M|max (full) : {scale:.3e} Ha")
    print(f"  max abs diff  : {abs_err:.3e} Ha")
    print(f"  max rel diff  : {rel_err:.3e}")
    print(f"  Frobenius rel : {rel_fro:.3e}")
    print(f"  -> {'PASS' if ok else 'FAIL'}  (tol={args.tol:g}); residual = finite-size error of the small supercell")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
