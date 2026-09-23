#!/usr/bin/env python3
"""R1 — génère nspin1/relax.in et nspin2/relax.in à partir du scf.in 9x9 défaut du ch. 4.
Modifs : calculation='relax' (cellule fixe, bfgs, nstep 200, forc 1e-4, etot 1e-5), nosym/noinv,
perturbation initiale : les deux voisins de la lacune 64 et 80 rapprochés de 0.03 Å chacun le long
de leur droite (dans le plan) ; nspin2 : espèce C1 (même pseudo) pour les voisins 64, 80, 81,
starting_magnetization(C1)=0.5. Tout le reste (cutoffs, Γ, mv 0.01 Ry, conv_thr 1e-10, pseudo) inchangé."""
import re, numpy as np, sys, os
SRC = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/super_cell/9x9/defective/scf.in"
SCRATCH = "/home/gregb26/links/scratch/qe_tmp/vacancy_relaxed"
BOHR = 0.529177210903
DISP_A = 0.03                       # Å, par atome
PAIR = (64, 80)                     # indices 1-based dans l'input défaut (voisins de la lacune)
THIRD = 81
NEIGH = (64, 80, 81)

txt = open(SRC).read()
head, rest = txt.split("CELL_PARAMETERS bohr")
cell_txt, pos_txt = rest.split("ATOMIC_POSITIONS crystal")
pos_txt, kp = pos_txt.split("K_POINTS")
A = np.array([[float(v) for v in l.split()] for l in cell_txt.strip().splitlines()])
pos = [l.split() for l in pos_txt.strip().splitlines()]
X = np.array([[float(v) for v in p[1:4]] for p in pos])
assert len(X) == 161
# --- vérification des voisins (distance minimale à la lacune 13/27,13/27)
rv = np.array([13/27, 13/27, 0]) @ A
d = X @ A - rv; d -= np.round(d @ np.linalg.inv(A)) @ A
nn = tuple(sorted(np.argsort(np.linalg.norm(d, axis=1))[:3] + 1))
assert nn == NEIGH, nn
# --- perturbation
i, j = PAIR[0]-1, PAIR[1]-1
ri, rj = X[i] @ A, X[j] @ A
u = (rj - ri) / np.linalg.norm(rj - ri); u[2] = 0.0
dd = DISP_A / BOHR
Xp = X.copy()
Xp[i] = (ri + dd*u) @ np.linalg.inv(A)
Xp[j] = (rj - dd*u) @ np.linalg.inv(A)
log = []
log.append(f"d(64-80) idéal = {np.linalg.norm(rj-ri):.6f} bohr = {np.linalg.norm(rj-ri)*BOHR:.6f} A")
log.append(f"u (64->80) = {u}")
for k, name in ((i, 64), (j, 80)):
    dc = (Xp[k]-X[k]) @ A
    log.append(f"atome {name}: cart. {X[k]@A} -> {Xp[k]@A} bohr ; delta cart = {dc} bohr = {dc*BOHR} A (|d|={np.linalg.norm(dc)*BOHR:.5f} A) ; delta cryst = {Xp[k]-X[k]}")
log.append(f"d(64-80) perturbé = {np.linalg.norm((Xp[j]-Xp[i])@A)*BOHR:.6f} A")
open(os.path.join(os.path.dirname(__file__) or '.', "perturbation.log"), "w").write("\n".join(log)+"\n")
print("\n".join(log))

def build(nspin):
    h = head
    h = h.replace("! ground state SCF of graphene 9x9 supercell", f"! R1 — relax (cellule fixe) de la supercellule 9x9 avec lacune, nspin={nspin}")
    h = h.replace("calculation = 'scf'       ! SCF ground state calculation", "calculation = 'relax'     ! ions relaxés, cellule fixe")
    h = h.replace("prefix      = 'defect_9x9_d'", f"prefix      = 'vac_9x9_relax_nspin{nspin}'")
    h = h.replace("outdir      = '/home/gregb26/links/scratch/qe_tmp/defect_9x9_d'", f"outdir      = '{SCRATCH}/nspin{nspin}'")
    h = h.replace("  pseudo_dir  =", "  forc_conv_thr = 1.0d-4\n  etot_conv_thr = 1.0d-5\n  nstep       = 200\n  tprnfor     = .true.\n  pseudo_dir  =")
    h = h.replace("  assume_isolated = '2D' ! helps for 2d systems\n", "  assume_isolated = '2D' ! helps for 2d systems\n  nosym       = .true.\n  noinv       = .true.\n")
    if nspin == 2:
        h = h.replace("  ntyp        = 1        ! number of atom types", "  ntyp        = 2        ! C (158 atomes) + C1 (3 voisins de la lacune, même pseudo)")
        h = h.replace("  noinv       = .true.\n", "  noinv       = .true.\n  nspin       = 2\n  starting_magnetization(1) = 0.0   ! C\n  starting_magnetization(2) = 0.5   ! C1 (atomes 64, 80, 81)\n")
    h = h.replace("  mixing_beta      =  0.3\n/", "  mixing_beta      =  0.3\n/\n\n&IONS\n  ion_dynamics = 'bfgs'\n/")
    h = h.replace("  C 12.011 C.upf ! atom type and pseudo potential", "  C  12.011 C.upf ! atom type and pseudo potential" + ("\n  C1 12.011 C.upf ! voisins de la lacune (64, 80, 81), même pseudo" if nspin == 2 else ""))
    lines = []
    for k, p in enumerate(pos):
        sp = "C1" if (nspin == 2 and k+1 in NEIGH) else "C"
        lines.append(f"  {sp:<2} {Xp[k,0]:.10f}  {Xp[k,1]:.10f}  {Xp[k,2]:.10f}")
    return h + "CELL_PARAMETERS bohr" + cell_txt + "ATOMIC_POSITIONS crystal\n" + "\n".join(lines) + "\n\nK_POINTS" + kp

for n in (1, 2):
    out = os.path.join(os.path.dirname(__file__) or '.', f"nspin{n}", "relax.in")
    open(out, "w").write(build(n))
    print("écrit", out)
