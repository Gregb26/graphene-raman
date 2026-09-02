#!/usr/bin/env python
"""
migrate_M_norm.py
    One-time migration of the pre-existing M_ed_<N>.npy (built before the 1/N_cells Bloch
    normalization was added to compute_M.py). For each raw matrix it:
        * tags the raw file as bloch_norm='unit_cell' (so it can never be mistaken for normalized),
        * writes M_ed_<N>_norm.npy = M / N_cells tagged bloch_norm='supercell'.
    Refuses to touch anything already tagged 'supercell' (double-application guard).

    Usage: python scripts/migrate_M_norm.py
"""
import glob
import os
import re

import numpy as np

from electron_defect_interaction.io import matrix_io

SIZE_RE = re.compile(r"M_ed_(\d+x\d+)\.npy$")


def main():
    for f in sorted(glob.glob("results/M/M_ed_*x*.npy")):
        if f.endswith("_norm.npy"):
            continue
        m = SIZE_RE.search(os.path.basename(f))
        if not m:
            continue
        size = m.group(1)
        meta = matrix_io.read_manifest(f)
        if meta and meta.get("bloch_norm") == matrix_io.SUPERCELL:
            print(f"SKIP {f}: already bloch_norm=supercell (guard against double 1/N_cells)")
            continue

        M = np.load(f)
        N_cells = M.shape[1]
        # tag the raw file so future loads know it is un-normalized
        if meta is None:
            matrix_io.save_M(f, M, matrix_io.UNIT_CELL, note="raw, pre-migration")
        out = f"results/M/M_ed_{size}_norm.npy"
        matrix_io.save_M(out, M / N_cells, matrix_io.SUPERCELL,
                         N_cells=int(N_cells), source=os.path.basename(f),
                         note="migrated: applied 1/N_cells to a unit-cell-normalized M")
        print(f"OK   {out}  N_cells={N_cells}  max|M|:{np.max(np.abs(M)):.3e}->{np.max(np.abs(M))/N_cells:.3e}")


if __name__ == "__main__":
    main()
