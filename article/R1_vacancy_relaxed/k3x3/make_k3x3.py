#!/usr/bin/env python3
"""R1b — quatre scf 3x3x1 (Γ-centré, nosym/noinv conservés) sur les géométries idéale et relaxées de R1.
Positions relaxées extraites du bloc 'Begin/End final coordinates' de relax.out (aucune recopie manuelle)."""
import re, os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R1 = os.path.dirname(D)
SRC = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/super_cell/9x9/defective/scf.in"
SCRATCH = "/home/gregb26/links/scratch/qe_tmp/vacancy_relaxed/k3x3"
BOHR = 0.529177210903; NEIGH = (64, 80, 81)
POSRE = r"^\s*(C1?)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$"

def split_input(txt):
    head, rest = txt.split("CELL_PARAMETERS bohr"); cell_txt, pos_txt = rest.split("ATOMIC_POSITIONS crystal")
    pos_txt, kp = pos_txt.split("K_POINTS")
    cell = np.array([[float(v) for v in l.split()] for l in cell_txt.strip().splitlines()])
    pos = [(m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in (re.match(POSRE, l) for l in pos_txt.strip().splitlines()) if m]
    return head, cell_txt, cell, pos

def final_positions(out):
    txt = open(out).read()
    blk = re.search(r"Begin final coordinates.*?ATOMIC_POSITIONS \(crystal\)\n(.*?)End final coordinates", txt, re.S).group(1)
    return [(m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in (re.match(POSRE, l) for l in blk.strip().splitlines()) if m]

def dist(pos, cell, i, j):
    d = (np.array(pos[j-1][1:]) - np.array(pos[i-1][1:])); d -= np.round(d); return np.linalg.norm(d @ cell) * BOHR

def build(nspin, pos, name):
    head, cell_txt, cell, _ = split_input(open(f"{R1}/nspin{nspin}/relax.in").read())
    h = head.replace(f"! R1 — relax (cellule fixe) de la supercellule 9x9 avec lacune, nspin={nspin}", f"! R1b — scf 3x3x1 ({name}), nspin={nspin}")
    h = h.replace("calculation = 'relax'     ! ions relaxés, cellule fixe", "calculation = 'scf'")
    h = h.replace(f"prefix      = 'vac_9x9_relax_nspin{nspin}'", f"prefix      = 'k3x3_{name}'")
    h = h.replace(f"outdir      = '{'/home/gregb26/links/scratch/qe_tmp/vacancy_relaxed'}/nspin{nspin}'", f"outdir      = '{SCRATCH}/{name}'")
    h = h.replace("  forc_conv_thr = 1.0d-4\n  etot_conv_thr = 1.0d-5\n  nstep       = 200\n", "")
    h = h.replace("&IONS\n  ion_dynamics = 'bfgs'\n/\n\n", "")
    assert "calculation = 'scf'" in h and f"{SCRATCH}/{name}" in h and "nstep" not in h and "&IONS" not in h and "tprnfor     = .true." in h and "nosym       = .true." in h
    if nspin == 2:
        assert "nspin       = 2" in h and "C1 12.011 C.upf" in h
        pos = [("C1" if k+1 in NEIGH else "C", x, y, z) for k, (_, x, y, z) in enumerate(pos)]
    else:
        pos = [("C", x, y, z) for (_, x, y, z) in pos]
    assert len(pos) == 161
    lines = [f"  {s:<2} {x:.10f}  {y:.10f}  {z:.10f}" for s, x, y, z in pos]
    txt = h + "CELL_PARAMETERS bohr" + cell_txt + "ATOMIC_POSITIONS crystal\n" + "\n".join(lines) + "\n\nK_POINTS automatic\n  3 3 1 0 0 0\n"
    os.makedirs(f"{D}/{name}", exist_ok=True); open(f"{D}/{name}/scf.in", "w").write(txt)
    return cell, pos

_, _, cell, ideal = split_input(open(SRC).read()); assert len(ideal) == 161
rel1 = final_positions(f"{R1}/nspin1/relax.out"); rel2 = final_positions(f"{R1}/nspin2/relax.out")
assert len(rel1) == 161 and len(rel2) == 161
assert [s for s, *_ in rel2].count("C1") == 3 and [rel2[i-1][0] for i in NEIGH] == ["C1"]*3
for name, nspin, pos, dexp in (("ideal_nspin1", 1, ideal, 2.46585), ("ideal_nspin2", 2, ideal, 2.46585), ("relax1_nspin1", 1, rel1, 2.19194), ("relax2_nspin2", 2, rel2, 1.98798)):
    c, p = build(nspin, pos, name)
    d = dist(p, c, 64, 80); assert abs(d - dexp) < 5e-5, (name, d)
    print(f"{name}: nat={len(p)} d(64-80)={d:.5f} A  C1={[s for s,*_ in p].count('C1')}")
