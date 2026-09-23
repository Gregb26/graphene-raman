#!/usr/bin/env python3
"""R1b — extraction brute : énergies (F, E, -TS), forces, magnétisation, itérations, durée des scf 3x3x1 et de R1 (Γ) ;
différences demandées en Ry et meV avec F et avec E ; charges/moments de Löwdin (projwfc.x) par atome."""
import re, os, sys, glob, numpy as np
RY_MEV = 13605.693122994
D = os.path.dirname(os.path.abspath(__file__)); R1 = os.path.dirname(D)
NEIGH = (64, 80, 81)
SUBA = None  # sous-réseau A = atomes dont (x*27) mod 3 == 1 dans l'input idéal (lacune sur A)

def sublattices():
    src = "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/super_cell/9x9/defective/scf.in"
    pos = [l.split() for l in open(src) if re.match(r"^\s+C\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*$", l)]
    A, B = [], []
    for k, p in enumerate(pos):
        r = round(float(p[1]) * 27) % 3
        (A if r == 1 else B).append(k + 1)
    assert len(A) == 80 and len(B) == 81, (len(A), len(B))
    assert all(a in B for a in NEIGH)
    return A, B

def parse_scf(path):
    t = open(path).read()
    g = lambda pat, cast=float: (cast(re.findall(pat, t, re.M)[-1]) if re.findall(pat, t, re.M) else None)
    F = g(r"^!\s+total energy\s+=\s+([-\d.]+) Ry")
    E = g(r"internal energy E=F\+TS\s+=\s+([-\d.]+) Ry")
    TS = g(r"smearing contrib\. \(-TS\)\s+=\s+([-\d.]+) Ry")
    Ftot = g(r"Total force =\s+([\d.Ee+-]+)")
    blocks = re.findall(r"Forces acting on atoms \(cartesian axes, Ry/au\):\n\n((?:\s+atom.*\n)+)", t)
    fmax = None
    if blocks:
        Fm = np.array([[float(v) for v in l.split("=")[1].split()] for l in blocks[-1].strip().splitlines()])
        fmax = np.abs(Fm).max(); fnorm = np.linalg.norm(Fm, axis=1); imax = int(np.argmax(fnorm)) + 1
    mt = g(r"total magnetization\s+=\s+([-\d.]+) Bohr mag/cell"); ma = g(r"absolute magnetization\s+=\s+([-\d.]+) Bohr mag/cell")
    it = g(r"convergence has been achieved in\s+(\d+) iterations", int)
    wall = re.findall(r"PWSCF\s+:\s+(.*?)\s+CPU\s+(.*?)\s+WALL", t)
    nk = g(r"number of k points=\s+(\d+)", int)
    ef = g(r"the Fermi energy is\s+([-\d.]+) ev")
    done = "JOB DONE" in t
    return dict(F=F, E=E, TS=TS, Ftot=Ftot, fmax=fmax, fnorm_max=(fnorm.max() if blocks else None), iatom_fmax=(imax if blocks else None),
                mt=mt, ma=ma, it=it, wall=(wall[-1][1] if wall else None), nk=nk, ef=ef, done=done)

def parse_lowdin(path):
    """Retourne {atome: (charge totale, up, down, polarisation)} depuis le bloc 'Lowdin Charges' de projwfc.out."""
    t = open(path).read()
    res = {}
    for m in re.finditer(r"Atom #\s+(\d+): total charge =\s+([\d.]+),.*?spin up\s+=\s+([\d.]+),.*?spin down\s+=\s+([\d.]+),.*?polarization =\s+([-\d.]+)", t, re.S):
        i = int(m.group(1))
        if i not in res: res[i] = (float(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)))
    tot = re.search(r"Spilling Parameter:\s+([\d.]+)", t)
    return res, (float(tot.group(1)) if tot else None)

if __name__ == "__main__":
    runs = {"R1 nspin1 (Γ, relax)": f"{R1}/nspin1/relax.out", "R1 nspin2 (Γ, relax)": f"{R1}/nspin2/relax.out",
            "ideal_nspin1 (3x3)": f"{D}/ideal_nspin1/scf.out", "ideal_nspin2 (3x3)": f"{D}/ideal_nspin2/scf.out",
            "relax1_nspin1 (3x3)": f"{D}/relax1_nspin1/scf.out", "relax2_nspin2 (3x3)": f"{D}/relax2_nspin2/scf.out"}
    R = {k: (parse_scf(v) if os.path.exists(v) else None) for k, v in runs.items()}
    print("| calcul | JOB DONE | nk | F = ! total energy (Ry) | E = F+TS (Ry) | −TS (Ry) | force max comp. (Ry/bohr) | |F_atome| max (atome) | force totale | m_tot / m_abs (µB) | E_F (eV) | itér. SCF | WALL |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, r in R.items():
        if r is None: print(f"| {k} | absent |"); continue
        fm = f"{r['fmax']:.2e}" if r['fmax'] is not None else "—"; fn = f"{r['fnorm_max']:.2e} ({r['iatom_fmax']})" if r['fnorm_max'] is not None else "—"
        mag = f"{r['mt']:+.2f} / {r['ma']:.2f}" if r['mt'] is not None else "—"
        f8 = lambda x: f"{x:.8f}" if x is not None else "—"
        print(f"| {k} | {r['done']} | {r['nk']} | {f8(r['F'])} | {f8(r['E'])} | {f8(r['TS'])} | {fm} | {fn} | {f8(r['Ftot'])} | {mag} | {r['ef']} | {r['it']} | {r['wall']} |")
    def dE(a, b, key):
        if R[a] is None or R[b] is None or R[a][key] is None or R[b][key] is None: return "—"
        d = R[a][key] - R[b][key]; return f"{d:+.8f} Ry = {d*RY_MEV:+.2f} meV"
    ideal_gamma = parse_scf("/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/defects/super_cell/9x9/defective/scf.out")
    R["ch4 idéal (Γ, 4 sym)"] = ideal_gamma
    print("\n| différence | avec F | avec E |\n|---|---|---|")
    for lab, a, b in (("ΔE_spin(Γ) = R1 nspin2 − R1 nspin1", "R1 nspin2 (Γ, relax)", "R1 nspin1 (Γ, relax)"),
                      ("ΔE_spin(3×3) = relax2_nspin2 − relax1_nspin1", "relax2_nspin2 (3x3)", "relax1_nspin1 (3x3)"),
                      ("ΔE_spin,idéal(3×3) = ideal_nspin2 − ideal_nspin1", "ideal_nspin2 (3x3)", "ideal_nspin1 (3x3)"),
                      ("ΔE_relax,1(3×3) = relax1_nspin1 − ideal_nspin1", "relax1_nspin1 (3x3)", "ideal_nspin1 (3x3)"),
                      ("ΔE_relax,2(3×3) = relax2_nspin2 − ideal_nspin2", "relax2_nspin2 (3x3)", "ideal_nspin2 (3x3)"),
                      ("ΔE_relax(Γ) = R1 nspin1 − scf ch. 4 idéal", "R1 nspin1 (Γ, relax)", "ch4 idéal (Γ, 4 sym)")):
        print(f"| {lab} | {dE(a, b, 'F')} | {dE(a, b, 'E')} |")
    A, B = sublattices()
    P = {"Γ (R1 nspin2)": f"{R1}/projwfc/gamma_nspin2/projwfc.out", "3×3 (relax2_nspin2)": f"{R1}/projwfc/k3x3_relax2_nspin2/projwfc.out"}
    L = {k: (parse_lowdin(v) if os.path.exists(v) else (None, None)) for k, v in P.items()}
    print("\n| moment de Löwdin (µB = up − down) | " + " | ".join(P) + " |\n|---|" + "---|" * len(P))
    rows = [("atome 81", [81]), ("atome 64", [64]), ("atome 80", [80]), ("somme sous-réseau A (80 atomes, celui de la lacune)", A), ("somme sous-réseau B (81 atomes)", B), ("total (161 atomes)", A + B)]
    for lab, idx in rows:
        cells = []
        for k in P:
            d, _ = L[k]
            cells.append("—" if not d else f"{sum(d[i][3] for i in idx if i in d):+.4f}")
        print(f"| {lab} | " + " | ".join(cells) + " |")
    for k in P:
        d, sp = L[k]
        if d: print(f"| charge de Löwdin totale, {k} : {sum(v[0] for v in d.values()):.3f} e (161 atomes, {len(d)} entrées) ; spilling {sp} |")
