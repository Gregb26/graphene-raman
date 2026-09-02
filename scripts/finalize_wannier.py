#!/usr/bin/env python
"""
finalize_wannier.py
    Post-run finalize for the unit-cell wannierizations (run once the wann_<N> SLURM jobs finish).
    For each size it: reports the .wout spreads and the interpolated-vs-DFT band error (raw numbers,
    no pass/fail judgement), copies _tb.dat/_u.mat/_u_dis.mat/.wout into the project (git-trackable),
    and writes the gauge-consistency manifest. Then prints the compute_spectral_wannier command.

    Usage: python scripts/finalize_wannier.py [5x5 7x7 8x8]
"""
import os
import re
import shutil
import sys

import numpy as np

from electron_defect_interaction.io.wannier_io import read_w90_HR, read_w90_mat
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.io.wannier_provenance import write_wannier_manifest

GRAPHENE = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/unit_cell"
DEST = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects/wannier"


def parse_spreads(wout):
    txt = open(wout).read()
    tot = re.findall(r"Omega Total\s*=\s*([-\d.]+)", txt)
    spr = re.findall(r"WF centre and spread\s+\d+\s+\([^)]*\)\s+([-\d.]+)", txt)
    conv = "Wannierisation convergence criteria satisfied" in txt or "Convergence criteria satisfied" in txt
    return (float(tot[-1]) if tot else None), [float(x) for x in spr[-5:]], conv


def read_eig(path):
    d = np.loadtxt(path)                       # columns: band, kidx, energy(eV)
    nb, nk = int(d[:, 0].max()), int(d[:, 1].max())
    E = np.zeros((nk, nb))
    for row in d:
        E[int(row[1]) - 1, int(row[0]) - 1] = row[2]
    return E


def finalize(N):
    wd = f"{GRAPHENE}/{N}"
    tb, u, ud = f"{wd}/wannier_tb.dat", f"{wd}/wannier_u.mat", f"{wd}/wannier_u_dis.mat"
    wout, eig = f"{wd}/wannier.wout", f"{wd}/wannier.eig"
    for f in (tb, u, wout):
        if not os.path.exists(f):
            print(f"[{N}] MISSING {os.path.basename(f)} -- wannier job not finished/failed; skipping.")
            return
    tot, spr, conv = parse_spreads(wout)
    print(f"[{N}] spreads: Omega_total = {tot} Ang^2 | per-WF = {spr} | converged flag = {conv}")

    # interpolated bands vs DFT (wannier.eig), on the u.mat k-grid (same win order)
    if os.path.exists(eig):
        Hwr, Rw, nd = read_w90_HR(tb)
        _, kU = read_w90_mat(u)
        _, Ew, _ = Hwr_to_Hwk(Hwr, Rw, kU, ndegen=nd)      # (nk, nw) eV
        Edft = read_eig(eig)                               # (nk, nb) eV
        errs = np.array([np.min(np.abs(Edft[ik] - Ew[ik, n]))
                         for ik in range(len(kU)) for n in range(Ew.shape[1])])
        # pi region: wannier bands within +-2.5 eV of the median crossing (raw, no judgement)
        mid = np.median(Ew)
        near = np.array([np.min(np.abs(Edft[ik] - Ew[ik, n]))
                         for ik in range(len(kU)) for n in range(Ew.shape[1])
                         if abs(Ew[ik, n] - mid) < 2.5])
        print(f"[{N}] band error interp-vs-DFT: global max={errs.max()*1e3:.2f} meV median={np.median(errs)*1e3:.2f} meV; "
              f"pi-region(|E-med|<2.5eV) max={near.max()*1e3:.2f} meV" if near.size else f"[{N}] band error global max={errs.max()*1e3:.2f} meV")
    else:
        print(f"[{N}] no wannier.eig -- skipping band comparison.")

    # rapatriate + manifest (git-trackable project space)
    d = f"{DEST}/{N}"
    os.makedirs(d, exist_ok=True)
    for f in (tb, u, ud, wout, eig):
        if os.path.exists(f):
            shutil.copy(f, d)
    man = f"{d}/wannier_manifest.json"
    write_wannier_manifest(man, f"{d}/wannier_tb.dat", f"{d}/wannier_u.mat",
                           f"{d}/wannier_u_dis.mat" if os.path.exists(f"{d}/wannier_u_dis.mat") else None)
    print(f"[{N}] copied to {d} + manifest {os.path.basename(man)}")
    print(f"[{N}] driver: python scripts/compute_spectral_wannier.py --size {N} "
          f"--manifest {man} --grids 60,120,240 --etas 0.05,0.02,0.01 --rcut 0,1,2 --nk-int 300")


def main():
    sizes = sys.argv[1:] or ["5x5", "7x7", "8x8"]
    for N in sizes:
        finalize(N)


if __name__ == "__main__":
    main()
