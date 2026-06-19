"""
run_test_A_batch.py
    Batch version of Test A (Kohn-Sham Hamiltonian reconstruction) for every graphene
    unit cell. For each cell size it reconstructs H(k) = T + V_loc + V^NL and compares
    its diagonal to the QE eigenvalues eps_nk.

    Unlike the pass/fail check in test_ks_reconstruction.py, this saves the full
    per-state deviation  dev = diag(H) - eps  (shape nb x nk) for each cell, so the
    agreement can be plotted. Output goes to results/test_A/.

    Run from the project root (paths are relative). Intended for a SLURM batch job
    (see scripts/submit_test_A.sh), but also runs interactively.

    Usage:  python -u scripts/run_test_A_batch.py [size ...]
            (default: all sizes 5x5 .. 12x12)
"""
import os
import sys

import numpy as np

sys.path.insert(0, "scripts")
from test_ks_reconstruction import reconstruct_ks_hamiltonian
from electron_defect_interaction.io import qe_io

SIZES = sys.argv[1:] or ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]
OUTDIR = "results/test_A"
TOL = 1e-5


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    summary = []
    for size in SIZES:
        uc = f"data/graphene/unit_cell/qe/defect_{size}.save"
        pot_uc = f"{uc}/Vks_uc"
        upf = f"{uc}/C.upf"

        print(f"[{size}] reconstructing H(k) ...", flush=True)
        H, eps = reconstruct_ks_hamiltonian(uc, pot_uc, upf)   # H (nk, nb, nb), eps (nb, nk)
        nk, nb, _ = H.shape

        diag = np.real(np.einsum("kii->ki", H)).T              # (nb, nk)
        eye = np.einsum("ki,ij->kij", np.einsum("kii->ki", H), np.eye(nb))
        offdiag_max = float(np.max(np.abs(H - eye)))           # how diagonal H really is
        dev = diag - eps                                       # (nb, nk)  H_reconstruit - eps
        shift = float(np.mean(dev))                            # constant offset (gauge of V_loc mean)
        dev_shifted = dev - shift                              # deviation after removing the offset
        max_dev = float(np.max(np.abs(dev_shifted)))
        k_red = qe_io.get_k_red(uc)                            # (nk, 3) for plotting context

        np.savez(
            os.path.join(OUTDIR, f"ks_recon_{size}.npz"),
            size=size, nb=nb, nk=nk,
            diag=diag, eps=eps,
            dev=dev, dev_shifted=dev_shifted,
            shift=shift, offdiag_max=offdiag_max,
            k_red=k_red, tol=TOL,
        )

        ok = offdiag_max < TOL and abs(shift) < TOL and max_dev < TOL
        summary.append((size, nb, nk, offdiag_max, shift, max_dev, "PASS" if ok else "FAIL"))
        print(f"[{size}] offdiag={offdiag_max:.2e}  shift={shift:.2e}  "
              f"max|dev|={max_dev:.2e}  -> {'PASS' if ok else 'FAIL'}", flush=True)

    header = f"{'size':>6} {'nb':>4} {'nk':>5} {'offdiag':>11} {'shift':>11} {'max|dev|':>11}  result"
    lines = [header]
    for s in summary:
        lines.append(f"{s[0]:>6} {s[1]:>4} {s[2]:>5} {s[3]:>11.2e} {s[4]:>11.2e} {s[5]:>11.2e}  {s[6]}")
    with open(os.path.join(OUTDIR, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    all_ok = all(s[6] == "PASS" for s in summary)
    print("\n".join(lines))
    print(f"RESULT: {'PASS' if all_ok else 'FAIL'} (tol={TOL:.0e} Ha)")
    print(f"Saved per-cell .npz + summary.txt to {OUTDIR}/")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
