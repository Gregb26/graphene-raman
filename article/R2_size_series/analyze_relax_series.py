#!/usr/bin/env python3
"""R2 phase 2 — généralisation de 9x9/analyze_relax.py (R1) à la série 5x5→12x12.
Usage : analyze_relax_series.py            -> analyse brute complète (comme analyse_2026-09-22.txt de R1), une section par (N, spin)
        analyze_relax_series.py --md       -> tableaux Markdown pour R2_rapport.md
Sorties : chiffres bruts, aucune interprétation. Le 9x9 est relu depuis ../9x9/ (R1) avec le même code.
Références : F = « ! total energy » (Ry) ; E idéal = F du scf.out du ch. 4 (super_cell/NxN/defective/scf.out, nspin1) ;
distances par image minimale ; déplacements vs positions idéales du scf.in du ch. 4 ; spglib via ASE, symprec 1e-3 Å."""
import re, sys, os, subprocess, numpy as np
from ase import Atoms
import spglib

BOHR = 0.529177210903; RY_EV = 13.605693122994
HERE = os.path.dirname(os.path.abspath(__file__))
QE = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects"
SYMPREC_A = 1e-3
SIZES = [5, 6, 7, 8, 9, 10, 11, 12]

def params():
    """Lit R2_phase0_table.txt : N -> (sous-réseau, s_red lacune, voisins, paire, troisième)."""
    P = {}
    for l in open(os.path.join(HERE, "R2_phase0_table.txt")):
        m = re.match(r"\s*(\d+)\s+([AB])\s+([\d.]+)\s+([\d.]+)\s+\(([\d, ]+)\)\s+\(([\d, ]+)\)\s+(\d+)\s+(\d+)", l)
        if m:
            N = int(m.group(1)); neigh = tuple(int(v) for v in m.group(5).split(","))
            pair = tuple(int(v) for v in m.group(6).split(",")); third = int(m.group(7))
            P[N] = dict(sub=m.group(2), vac=np.array([float(m.group(3)), float(m.group(4)), 0.0]), neigh=neigh, pair=pair, third=third, nbnd=int(m.group(8)))
    return P

def workdir(N, spin):
    return os.path.join(HERE, f"{N}x{N}", spin) if N != 9 else os.path.join(HERE, "..", "9x9", spin)

def projwfc_out(N):
    return os.path.join(HERE, f"{N}x{N}", "nspin2", "projwfc", "projwfc.out") if N != 9 else os.path.join(HERE, "..", "9x9", "projwfc", "gamma_nspin2", "projwfc.out")

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
    fmax = []
    for blk in re.findall(r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n((?:\s+atom.*\n)+)", txt):
        F = np.array([[float(v) for v in l.split("=")[1].split()] for l in blk.strip().splitlines()])
        fmax.append(np.abs(F).max())
    steps = []
    for blk in re.findall(r"ATOMIC_POSITIONS \(crystal\)\n((?:\s*C1?\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+.*\n)+)", txt):
        steps.append(np.array([[float(v) for v in l.split()[1:4]] for l in blk.strip().splitlines()]))
    magt = [float(m) for m in re.findall(r"total magnetization\s+=\s+([-\d.]+) Bohr mag/cell", txt)]
    maga = [float(m) for m in re.findall(r"absolute magnetization\s+=\s+([-\d.]+) Bohr mag/cell", txt)]
    mag_steps = []
    for m in re.finditer(r"total magnetization\s+=\s+([-\d.]+) Bohr mag/cell\n\s+absolute magnetization\s+=\s+([-\d.]+) Bohr mag/cell\n\n\s+convergence has been achieved", txt):
        mag_steps.append((float(m.group(1)), float(m.group(2))))
    fin = re.search(r"Begin final coordinates(.*?)End final coordinates", txt, re.S)
    final = None
    if fin:
        final = np.array([[float(v) for v in l.split()[1:4]] for l in fin.group(1).split("ATOMIC_POSITIONS (crystal)")[1].strip().splitlines() if re.match(r"\s*C", l)])
    conv = re.search(r"bfgs converged in\s+(\d+) scf cycles and\s+(\d+) bfgs steps", txt)
    wall = re.findall(r"PWSCF\s+:\s+(.*?)\s+CPU\s+(.*?)\s+WALL", txt)
    notconv = len(re.findall(r"convergence NOT achieved", txt))
    ef = re.findall(r"the Fermi energy is\s+([-\d.]+) ev", txt)
    TS = re.findall(r"smearing contrib\. \(-TS\)\s+=\s+([-\d.]+) Ry", txt)
    return dict(E=E, Ftot=Ftot, it=it, fmax=fmax, steps=steps, magt=magt, maga=maga, mag_steps=mag_steps, final=final, conv=conv,
                wall=(wall[-1][1] if wall else None), done=("JOB DONE" in txt), notconv=notconv, ef=(float(ef[-1]) if ef else None),
                TS=(float(TS[-1]) if TS else None))

def parse_lowdin(path):
    """{atome: dict(tot, up, dn, pol, pol_s, pol_p, pol_pz, pol_px, pol_py)} + (charge totale, spilling)."""
    t = open(path).read(); res = {}
    pat = re.compile(r"Atom #\s+(\d+): total charge =\s+([\d.]+), s =\s+([\d.]+), p =\s+([\d.]+),\s*\n"
                     r"\s+spin up\s+=\s+([\d.]+), s =\s+([\d.]+),\s*\n"
                     r"\s+spin up\s+=\s+[\d.]+, p =\s+([\d.]+), pz=\s+([\d.]+), px=\s+([\d.]+), py=\s+([\d.]+),\s*\n"
                     r"\s+spin down\s+=\s+([\d.]+), s =\s+([\d.]+),\s*\n"
                     r"\s+spin down\s+=\s+[\d.]+, p =\s+([\d.]+), pz=\s+([\d.]+), px=\s+([\d.]+), py=\s+([\d.]+),\s*\n"
                     r"\s+polarization =\s+([-\d.]+), s =\s+([-\d.]+), p =\s+([-\d.]+)")
    for m in pat.finditer(t):
        i = int(m.group(1)); g = [float(v) for v in m.groups()[1:]]
        if i in res: continue
        res[i] = dict(tot=g[0], up=g[3], dn=g[9], pol=g[15], pol_s=g[16], pol_p=g[17], pol_pz=g[6]-g[12], pol_px=g[7]-g[13], pol_py=g[8]-g[14])
    sp = re.search(r"Spilling Parameter:\s+([\d.]+)", t)
    return res, (float(sp.group(1)) if sp else None)

def sublattices(N, Xideal):
    """A = atomes dont round(3N x) mod 3 == 1 (règle s_uc = mod(N s, 1) : A ≈ 1/3, B ≈ 2/3), même règle que R1b/tag_vacancy_sublattice."""
    A, B = [], []
    for k, x in enumerate(Xideal[:, 0]):
        r = round(x * 3 * N) % 3
        assert r in (1, 2), (N, k, x, r)
        (A if r == 1 else B).append(k + 1)
    return A, B

def mic(d, cell):
    f = d @ np.linalg.inv(cell); f -= np.round(f); return f @ cell

def dist(X, cell, i, j):
    return np.linalg.norm(mic((X[j-1]-X[i-1]) @ cell, cell)) * BOHR

def spg(cell, X, vac, symprec_A):
    at = Atoms("C"*len(X), scaled_positions=X, cell=cell*BOHR, pbc=True)
    ds = spglib.get_symmetry_dataset((at.cell[:], at.get_scaled_positions(), at.numbers), symprec=symprec_A)
    at2 = Atoms("C"*len(X)+"He", scaled_positions=np.vstack([X, vac]), cell=cell*BOHR, pbc=True)
    ds2 = spglib.get_symmetry_dataset((at2.cell[:], at2.get_scaled_positions(), at2.numbers), symprec=symprec_A)
    return ds, ds2

def elapsed(d):
    try:
        jid = open(os.path.join(d, "JOBID")).read().strip()
        r = subprocess.run(["sacct", "-j", jid, "-n", "-X", "-o", "JobID,State,Elapsed"], capture_output=True, text=True).stdout.split()
        return jid, r[1], r[2]
    except Exception as e:
        return None, None, None

def analyze(N, spin, P, verbose=True):
    d = workdir(N, spin); p = P[N]; p1, p2 = p["pair"]; t3 = p["third"]
    cell, X0, _ = read_input_positions(os.path.join(d, "relax.in"))
    _, Xideal, _ = read_input_positions(f"{QE}/super_cell/{N}x{N}/defective/scf.in")
    o = parse_out(os.path.join(d, "relax.out"))
    Eideal = parse_out(f"{QE}/super_cell/{N}x{N}/defective/scf.out")["E"][-1]
    n = len(o["E"])
    Xf = o["final"] if o["final"] is not None else (o["steps"][-1] if o["steps"] else X0)
    tag = "finales" if o["final"] is not None else "DERNIER PAS DISPONIBLE (non final)"
    D = np.array([mic((Xf[i]-Xideal[i]) @ cell, cell) for i in range(len(Xf))]) * BOHR
    Dp = np.array([mic((Xf[i]-X0[i]) @ cell, cell) for i in range(len(Xf))]) * BOHR
    nrm = np.linalg.norm(D, axis=1); im = int(np.argmax(nrm)) + 1
    jid, state, elap = elapsed(d)
    ds0, ds0v = spg(cell, X0, p["vac"], SYMPREC_A); dsf, dsfv = spg(cell, Xf, p["vac"], SYMPREC_A)
    R = dict(N=N, spin=spin, nat=len(Xf), done=o["done"], conv=o["conv"], nscf=(int(o["conv"].group(1)) if o["conv"] else None),
             nbfgs=(int(o["conv"].group(2)) if o["conv"] else None), notconv=o["notconv"], E=o["E"][-1], Eideal=Eideal,
             dE_relax=(o["E"][-1]-Eideal), Ftot=o["Ftot"][-1], fmax=o["fmax"][-1], itmax=max(o["it"]), it_last=o["it"][-1],
             d_pair=(dist(Xideal,cell,p1,p2), dist(X0,cell,p1,p2), dist(Xf,cell,p1,p2)),
             d_13=dist(Xf,cell,p1,t3), d_23=dist(Xf,cell,p2,t3), d_13_ideal=dist(Xideal,cell,p1,t3),
             dmax=nrm.max(), imax=im, dmax_pert=np.linalg.norm(Dp,axis=1).max(), dz_max=np.abs(D[:,2]).max(),
             d3_inplane=np.linalg.norm(D[t3-1][:2]), dz3=D[t3-1][2], d_pair_inplane=(np.linalg.norm(D[p1-1][:2]), np.linalg.norm(D[p2-1][:2])),
             magt=(o["magt"][-1] if o["magt"] else None), maga=(o["maga"][-1] if o["maga"] else None),
             spg0=(ds0.international, ds0.number, len(ds0.rotations), ds0v.site_symmetry_symbols[len(X0)], ds0v.international, len(ds0v.rotations)),
             spgf=(dsf.international, dsf.number, len(dsf.rotations), dsfv.site_symmetry_symbols[len(Xf)], dsfv.international, len(dsfv.rotations)),
             wall=o["wall"], jid=jid, state=state, elapsed=elap, ef=o["ef"], TS=o["TS"], Xideal=Xideal, cell=cell, X0=X0, Xf=Xf, o=o, p=p, D=D)
    if verbose:
        print(f"\n===== {N}x{N} {spin} ({d}) : JOB DONE={o['done']}, sacct {jid} {state} {elap}, pas scf/BFGS = {n}, bfgs converged = {bool(o['conv'])}"
              + (f" ({R['nscf']} scf, {R['nbfgs']} bfgs)" if o['conv'] else "") + f", « convergence NOT achieved » = {o['notconv']}, WALL = {o['wall']}")
        print(f"E finale = {R['E']:.8f} Ry = {R['E']*RY_EV:.6f} eV ; E scf idéal ch. 4 = {Eideal:.8f} Ry ; ΔE_relax = {R['dE_relax']:.8f} Ry = {R['dE_relax']*RY_EV*1e3:+.3f} meV ; "
              f"Ftot final = {R['Ftot']} Ry/bohr ; |F|max composante final = {R['fmax']:.6g} Ry/bohr ; scf iter max/par pas = {R['itmax']} ; E_F = {R['ef']} eV ; -TS = {R['TS']} Ry")
        print(f"--- distances C-C entre voisins sous-coordonnés (A) ; lacune sous-réseau {p['sub']} à {p['vac'][:2]} ; voisins {p['neigh']} ; paire {p['pair']} ; troisième {t3}")
        for tagX, X in (("idéal", Xideal), ("initial perturbé", X0), (tag, Xf)):
            print(f"  {tagX:26s}: d({p1}-{p2})={dist(X,cell,p1,p2):.5f}  d({p1}-{t3})={dist(X,cell,p1,t3):.5f}  d({p2}-{t3})={dist(X,cell,p2,t3):.5f}")
        print(f"--- déplacement max vs idéal = {R['dmax']:.5f} A (atome {im}) ; vs initial perturbé = {R['dmax_pert']:.5f} A ; hors plan max |dz| = {R['dz_max']:.2e} A")
        for a in p["neigh"]:
            print(f"  atome {a}: |d|={np.linalg.norm(D[a-1]):.5f} A, d_inplane={np.linalg.norm(D[a-1][:2]):.5f} A, dz={D[a-1][2]:+.2e} A")
        print(f"--- atome {t3} : d({t3}-{p1})={R['d_13']:.5f} A, d({t3}-{p2})={R['d_23']:.5f} A, dz({t3})={R['dz3']:+.2e} A, d_inplane({t3})={R['d3_inplane']:.5f} A")
        print(f"--- suivi par pas : k, E (Ry), dE (meV), Ftot (Ry/bohr), Fmax comp., scf iter, d({p1}-{p2}) (A)" + (", m_tot, m_abs" if o["mag_steps"] else ""))
        Xs = [X0] + o["steps"]
        for k in range(n):
            Xk = Xs[k] if k < len(Xs) else Xs[-1]
            dE = (o["E"][k]-o["E"][k-1])*RY_EV*1e3 if k else 0.0
            line = f"  {k:3d} {o['E'][k]:.8f} {dE:+9.3f} {o['Ftot'][k] if k < len(o['Ftot']) else float('nan'):.6f} {o['fmax'][k] if k < len(o['fmax']) else float('nan'):.6f} {o['it'][k] if k < len(o['it']) else -1:4d} {dist(Xk,cell,p1,p2):.5f}"
            if o["mag_steps"] and k < len(o["mag_steps"]): line += f" {o['mag_steps'][k][0]:+.4f} {o['mag_steps'][k][1]:.4f}"
            print(line)
        if o["magt"]:
            print(f"--- magnétisation finale : totale = {R['magt']:+.4f} µB, absolue = {R['maga']:.4f} µB")
        print(f"--- spglib (symprec {SYMPREC_A} A)")
        for tagX, s in (("initial perturbé", R["spg0"]), (tag, R["spgf"])):
            print(f"  {tagX:26s}: supercellule {s[0]} (n° {s[1]}, {s[2]} op.) ; site lacune (atome factice) : {s[3]} ; groupe avec factice {s[4]} ({s[5]} op.)")
    return R

def lowdin(N, P, R, verbose=True):
    f = projwfc_out(N)
    if not os.path.exists(f): 
        if verbose: print(f"--- Löwdin {N}x{N} : {f} absent")
        return None
    L, spill = parse_lowdin(f); p = P[N]; A, B = sublattices(N, R["Xideal"])
    if not L or "JOB DONE" not in open(f).read():
        if verbose: print(f"--- Löwdin {N}x{N} : {f} incomplet (JOB DONE absent ou bloc Löwdin vide)")
        return None
    vac_sub, oth_sub = (A, B) if p["sub"] == "A" else (B, A)
    assert all(a in oth_sub for a in p["neigh"]), (N, p["neigh"])
    assert len(L) == R["nat"], (N, len(L), R["nat"])
    s = lambda idx: sum(L[i]["pol"] for i in idx)
    res = dict(m3=L[p["third"]]["pol"], m1=L[p["pair"][0]]["pol"], m2=L[p["pair"][1]]["pol"], sum_vac=s(vac_sub), sum_oth=s(oth_sub), total=s(A+B),
               n_vac=len(vac_sub), n_oth=len(oth_sub), charge=sum(L[i]["tot"] for i in A+B), spill=spill, comp3=L[p["third"]],
               sub_vac=p["sub"], sub_oth=("B" if p["sub"] == "A" else "A"))
    if verbose:
        print(f"--- Löwdin (projwfc.x, {f}) : {len(L)} atomes, charge totale {res['charge']:.3f} e, spilling {spill}")
        print(f"  moment atome {p['third']} (isolé) = {res['m3']:+.4f} µB ; atome {p['pair'][0]} = {res['m1']:+.4f} ; atome {p['pair'][1]} = {res['m2']:+.4f}")
        print(f"  somme sous-réseau {res['sub_vac']} (lacune, {res['n_vac']} atomes) = {res['sum_vac']:+.4f} ; somme sous-réseau {res['sub_oth']} ({res['n_oth']} atomes) = {res['sum_oth']:+.4f} ; total = {res['total']:+.4f} µB (pw.x : {R['magt']:+.2f})")
        c = res["comp3"]
        print(f"  composantes atome {p['third']} : s {c['pol_s']:.4f}, p {c['pol_p']:.4f} (pz {c['pol_pz']:.4f}, px {c['pol_px']:.4f}, py {c['pol_py']:.4f})")
    return res

def markdown(P, RES, LOW):
    ry = lambda x: f"{x:.8f}"
    print("## Tableau récapitulatif (une ligne par (N, spin) ; F = « ! total energy » ; 9x9 = R1 relu avec le même script)\n")
    print("| N | nat | K replié sur Γ | spin | job ID | état | durée (sacct) | pas BFGS / cycles SCF | SCF non conv. | F finale (Ry) | F scf idéal ch. 4 (Ry) | ΔE_relax (Ry) | ΔE_relax (meV) | force totale finale (Ry/bohr) | comp. max (Ry/bohr) | d(paire) idéal → perturbé → final (Å) | d(p1–3e) / d(p2–3e) final (Å) | dépl. max vs idéal (Å) (atome) | dépl. 3e dans le plan (Å) | hors plan max (Å) | m_tot / m_abs (µB) | spglib initial | spglib final |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for N in SIZES:
        for spin in ("nspin1", "nspin2"):
            R = RES[(N, spin)]; p = P[N]
            mag = f"{R['magt']:+.2f} / {R['maga']:.2f}" if R["magt"] is not None else "—"
            s0, sf = R["spg0"], R["spgf"]
            print(f"| {N} | {R['nat']} | {'oui' if N % 3 == 0 else 'non'} | {spin} | {R['jid']} | {R['state']} | {R['elapsed']} | {R['nbfgs']} / {R['nscf']} | {R['notconv']} | {ry(R['E'])} | {ry(R['Eideal'])} | {ry(R['dE_relax'])} | {R['dE_relax']*RY_EV*1e3:+.2f} | {R['Ftot']:.2e} | {R['fmax']:.2e} | "
                  f"{R['d_pair'][0]:.5f} → {R['d_pair'][1]:.5f} → {R['d_pair'][2]:.5f} | {R['d_13']:.5f} / {R['d_23']:.5f} | {R['dmax']:.5f} ({R['imax']}) | {R['d3_inplane']:.5f} | {R['dz_max']:.1e} | {mag} | "
                  f"{s0[0]} (n° {s0[1]}, {s0[2]} op.), site {s0[3]} | {sf[0]} (n° {sf[1]}, {sf[2]} op.), site {sf[3]} |")
    print("\n## ΔE_spin = F(nspin2) − F(nspin1) et moments de Löwdin (projwfc.x sur le .save nspin2, lwrite_overlaps = .false.)\n")
    print("| N | K replié sur Γ | lacune / paire / 3e | ΔE_spin (Ry) | ΔE_spin (meV) | m_tot / m_abs pw.x (µB) | Löwdin atome 3e (µB) | Löwdin paire (µB) | Löwdin somme sous-réseau de la lacune (µB, n) | Löwdin somme autre sous-réseau (µB, n) | Löwdin total (µB) | comp. 3e : s / p (pz, px, py) | charge Löwdin (e) | spilling |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for N in SIZES:
        R1_, R2_ = RES[(N, "nspin1")], RES[(N, "nspin2")]; p = P[N]; dE = R2_["E"] - R1_["E"]; L = LOW.get(N)
        lw = (f"{L['m3']:+.4f} | {L['m1']:+.4f} / {L['m2']:+.4f} | {L['sum_vac']:+.4f} ({L['sub_vac']}, {L['n_vac']}) | {L['sum_oth']:+.4f} ({L['sub_oth']}, {L['n_oth']}) | {L['total']:+.4f} | "
              f"{L['comp3']['pol_s']:.4f} / {L['comp3']['pol_p']:.4f} ({L['comp3']['pol_pz']:.4f}, {L['comp3']['pol_px']:.4f}, {L['comp3']['pol_py']:.4f}) | {L['charge']:.3f} | {L['spill']}") if L else "absent | | | | | | |"
        print(f"| {N} | {'oui' if N % 3 == 0 else 'non'} | {p['sub']} / ({p['pair'][0]}, {p['pair'][1]}) / {p['third']} | {ry(dE)} | {dE*RY_EV*1e3:+.2f} | {R2_['magt']:+.2f} / {R2_['maga']:.2f} | {lw} |")

if __name__ == "__main__":
    md = "--md" in sys.argv
    P = params(); RES = {}; LOW = {}
    for N in SIZES:
        for spin in ("nspin1", "nspin2"):
            RES[(N, spin)] = analyze(N, spin, P, verbose=not md)
        LOW[N] = lowdin(N, P, RES[(N, "nspin2")], verbose=not md)
        if not md:
            dE = (RES[(N,"nspin2")]["E"]-RES[(N,"nspin1")]["E"])
            print(f"\n===== {N}x{N} : E(nspin2) - E(nspin1) = {dE:.8f} Ry = {dE*RY_EV*1e3:+.3f} meV")
    if md: markdown(P, RES, LOW)
