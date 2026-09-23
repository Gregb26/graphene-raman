#!/usr/bin/env python3
"""R1 phase 2 — extraction des contrôles d'une relaxation pw.x (relax.out) de la lacune 9x9.
Usage : analyze_relax.py nspin1 [nspin2]   (répertoires contenant relax.in / relax.out)
Sorties : chiffres bruts (pas d'interprétation). spglib (ASE) sur positions initiales perturbées et finales."""
import re, sys, os, numpy as np
from ase import Atoms
import spglib

BOHR = 0.529177210903; RY_EV = 13.605693122994
SRC = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/super_cell/9x9/defective/scf.in"
NEIGH = (64, 80, 81); PAIR = (64, 80); THIRD = 81
VAC = np.array([13/27, 13/27, 0.0])
SYMPREC_A = 1e-3

def read_input_positions(path):
    txt = open(path).read()
    cell = np.array([[float(v) for v in l.split()] for l in txt.split("CELL_PARAMETERS bohr")[1].split("ATOMIC_POSITIONS")[0].strip().splitlines()])
    pos = txt.split("ATOMIC_POSITIONS crystal")[1].split("K_POINTS")[0].strip().splitlines()
    X = np.array([[float(v) for v in l.split()[1:4]] for l in pos]); sp = [l.split()[0] for l in pos]
    return cell, X, sp

def parse_out(path):
    txt = open(path).read()
    E = [float(m) for m in re.findall(r"^!\s+total energy\s+=\s+([-\d.]+) Ry", txt, re.M)]
    Ftot = [float(m) for m in re.findall(r"Total force =\s+([\d.Ee+-]+)", txt)]
    it = [int(m) for m in re.findall(r"convergence has been achieved in\s+(\d+) iterations", txt)]
    # max composante de force par pas (bloc "Forces acting on atoms")
    fmax = []
    for blk in re.findall(r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n((?:\s+atom.*\n)+)", txt):
        F = np.array([[float(v) for v in l.split("=")[1].split()] for l in blk.strip().splitlines()])
        fmax.append(np.abs(F).max())
    # positions après chaque pas
    steps = []
    for blk in re.findall(r"ATOMIC_POSITIONS \(crystal\)\n((?:\s*C1?\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+.*\n)+)", txt):
        steps.append(np.array([[float(v) for v in l.split()[1:4]] for l in blk.strip().splitlines()]))
    magt = [float(m) for m in re.findall(r"total magnetization\s+=\s+([-\d.]+) Bohr mag/cell", txt)]
    maga = [float(m) for m in re.findall(r"absolute magnetization\s+=\s+([-\d.]+) Bohr mag/cell", txt)]
    # magnétisation au dernier scf de chaque pas : celle qui précède chaque "!"
    mag_steps = []
    for m in re.finditer(r"total magnetization\s+=\s+([-\d.]+) Bohr mag/cell\n\s+absolute magnetization\s+=\s+([-\d.]+) Bohr mag/cell\n\n\s+convergence has been achieved", txt):
        mag_steps.append((float(m.group(1)), float(m.group(2))))
    fin = re.search(r"Begin final coordinates(.*?)End final coordinates", txt, re.S)
    final = None
    if fin:
        final = np.array([[float(v) for v in l.split()[1:4]] for l in fin.group(1).split("ATOMIC_POSITIONS (crystal)")[1].strip().splitlines() if re.match(r"\s*C", l)])
    conv = re.search(r"bfgs converged in\s+(\d+) scf cycles and\s+(\d+) bfgs steps", txt)
    nbfgs = re.findall(r"number of bfgs steps\s+=\s+(\d+)", txt)
    wall = re.search(r"PWSCF\s+:.*?([\d.]+h?[\d. ]*m?[\d.]*s) WALL", txt)
    done = "JOB DONE" in txt
    return dict(E=E, Ftot=Ftot, it=it, fmax=fmax, steps=steps, magt=magt, maga=maga, mag_steps=mag_steps,
                final=final, conv=conv, nbfgs=nbfgs, done=done, txt=txt)

def mic(d, cell):
    f = d @ np.linalg.inv(cell); f -= np.round(f); return f @ cell

def dist(X, cell, i, j):
    return np.linalg.norm(mic((X[j-1]-X[i-1]) @ cell, cell)) * BOHR

def spg(cell, X, symprec_A):
    at = Atoms("C"*len(X), scaled_positions=X, cell=cell*BOHR, pbc=True)
    ds = spglib.get_symmetry_dataset((at.cell[:], at.get_scaled_positions(), at.numbers), symprec=symprec_A)
    # site de la lacune : atome factice (Z=2) à la position idéale
    at2 = Atoms("C"*len(X)+"He", scaled_positions=np.vstack([X, VAC]), cell=cell*BOHR, pbc=True)
    ds2 = spglib.get_symmetry_dataset((at2.cell[:], at2.get_scaled_positions(), at2.numbers), symprec=symprec_A)
    return ds, ds2

def report(d):
    cell, X0, sp = read_input_positions(os.path.join(d, "relax.in"))
    _, Xideal, _ = read_input_positions(SRC)
    o = parse_out(os.path.join(d, "relax.out"))
    n = len(o["E"])
    print(f"\n===== {d} : JOB DONE={o['done']}, pas scf/BFGS = {n}, bfgs converged = {bool(o['conv'])}"
          + (f" ({o['conv'].group(1)} scf, {o['conv'].group(2)} bfgs)" if o['conv'] else ""))
    Xf = o["final"] if o["final"] is not None else (o["steps"][-1] if o["steps"] else X0)
    tag = "finales" if o["final"] is not None else "DERNIER PAS DISPONIBLE (non final)"
    print(f"E finale = {o['E'][-1]:.8f} Ry = {o['E'][-1]*RY_EV:.6f} eV ; Ftot final = {o['Ftot'][-1] if o['Ftot'] else 'nan'} Ry/bohr ; "
          f"|F|max composante final = {o['fmax'][-1] if o['fmax'] else 'nan'} Ry/bohr ; scf iter max/par pas = {max(o['it']) if o['it'] else 'nan'}")
    print("--- distances C-C entre voisins sous-coordonnés (A)")
    for tagX, X in (("idéal", Xideal), ("initial perturbé", X0), (tag, Xf)):
        print(f"  {tagX:26s}: d(64-80)={dist(X,cell,64,80):.5f}  d(64-81)={dist(X,cell,64,81):.5f}  d(80-81)={dist(X,cell,80,81):.5f}")
    D = np.array([mic((Xf[i]-Xideal[i]) @ cell, cell) for i in range(len(Xf))]) * BOHR
    Dp = np.array([mic((Xf[i]-X0[i]) @ cell, cell) for i in range(len(Xf))]) * BOHR
    im = np.argmax(np.linalg.norm(D, axis=1))
    print(f"--- déplacement max vs idéal = {np.linalg.norm(D,axis=1).max():.5f} A (atome {im+1}) ; vs initial perturbé = {np.linalg.norm(Dp,axis=1).max():.5f} A ; hors plan max |dz| = {np.abs(D[:,2]).max():.2e} A")
    for a in NEIGH:
        print(f"  atome {a}: |d|={np.linalg.norm(D[a-1]):.5f} A, d_inplane={np.linalg.norm(D[a-1][:2]):.5f} A, dz={D[a-1][2]:+.2e} A")
    print(f"--- atome 81 : d(81-64)={dist(Xf,cell,81,64):.5f} A, d(81-80)={dist(Xf,cell,81,80):.5f} A, dz(81)={D[80][2]:+.2e} A, d_inplane(81)={np.linalg.norm(D[80][:2]):.5f} A")
    print("--- suivi par pas : k, E (Ry), dE (meV), Ftot (Ry/bohr), Fmax comp., scf iter, d(64-80) (A)" + (", m_tot, m_abs" if o["mag_steps"] else ""))
    Xs = [X0] + o["steps"]
    for k in range(n):
        Xk = Xs[k] if k < len(Xs) else Xs[-1]
        dE = (o["E"][k]-o["E"][k-1])*RY_EV*1e3 if k else 0.0
        line = f"  {k:3d} {o['E'][k]:.8f} {dE:+9.3f} {o['Ftot'][k] if k < len(o['Ftot']) else float('nan'):.6f} {o['fmax'][k] if k < len(o['fmax']) else float('nan'):.6f} {o['it'][k] if k < len(o['it']) else -1:4d} {dist(Xk,cell,64,80):.5f}"
        if o["mag_steps"] and k < len(o["mag_steps"]): line += f" {o['mag_steps'][k][0]:+.4f} {o['mag_steps'][k][1]:.4f}"
        print(line)
    if o["magt"]:
        print(f"--- magnétisation finale : totale = {o['magt'][-1]:+.4f} µB, absolue = {o['maga'][-1]:.4f} µB")
    print(f"--- spglib (symprec {SYMPREC_A} A)")
    for tagX, X in (("initial perturbé", X0), (tag, Xf)):
        ds, ds2 = spg(cell, X, SYMPREC_A)
        i_vac = len(X)
        print(f"  {tagX:26s}: supercellule {ds.international} (n° {ds.number}, {len(ds.rotations)} op.) ; site lacune (atome factice) : {ds2.site_symmetry_symbols[i_vac]} ; groupe avec factice {ds2.international} ({len(ds2.rotations)} op.)")
    return o

if __name__ == "__main__":
    res = {}
    for d in sys.argv[1:]:
        res[os.path.basename(os.path.normpath(d))] = report(d)
    if "nspin1" in res and "nspin2" in res:
        dE = (res["nspin2"]["E"][-1]-res["nspin1"]["E"][-1])*RY_EV*1e3
        print(f"\n===== E(nspin2) - E(nspin1) = {dE:+.3f} meV")
