#!/usr/bin/env python
"""
make_inputs_r7_relax.py -- R7c (2026-09-26) : relaxation (cellule fixe, BFGS, nspin 1) des lacunes 15x15, 18x18, 21x21, 24x24, 27x27,
en partant de la geometrie relaxee de la 12x12 (R2, nspin1) transplantee.

Protocole R2 conserve pour les namelists (calculation='relax', forc_conv_thr 1e-4, etot_conv_thr 1e-5, nstep 200, tprnfor, nosym, noinv,
&IONS bfgs ; le reste = scf.in R7 de la meme taille : cutoffs, Gamma, mv 0.01 Ry, conv_thr 1e-10, nbnd R7, pseudo). Ce qui change par
rapport a R2 : la geometrie initiale n'est plus la maille ideale perturbee (paire rapprochee de 0,03 A) mais le champ de deplacements
u(r) = x_relax - x_ideal de la 12x12 R2 nspin1 (relax.out, coordonnees finales), transplante autour de la lacune de la cellule N x N.
La 12x12 a sa lacune sur le sous-reseau B, les cellules R7 sur A : le champ est transplante par inversion autour de la lacune
(r -> -r, u -> -u), qui echange A et B sur le nid d'abeille. Rayon de transplant R_CUT (< L12/2 = 14,795 A, unicite de l'image minimale) ;
au-dela, positions ideales.
Regressions (--check) : (1) auto-transplant 12x12 sans inversion = coordonnees finales R2 ; (2) namelists derives du scf.in 12x12 = relax.in
R2 12x12 nspin1 hors titre et positions ; (3) l'inversion envoie la maille ideale B-lacune sur la maille ideale A-lacune (12x12).
Usage : python make_inputs_r7_relax.py --check ; --write [--sizes 15,18,21,24,27] ; --profile (profil |u|(r) de la 12x12)
Sorties : super_cell_relaxed/series/NxN/nspin1/{relax.in, submit.relax}, super_cell_relaxed/series/NxN/transplant.log
"""
import os, re, sys, argparse
import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))
SC = os.path.join(GQ, "defects", "super_cell")
SER = os.path.join(GQ, "defects", "super_cell_relaxed", "series")
SCRATCH = "/home/gregb26/links/scratch/qe_tmp/vacancy_relaxed/series/{N}x{N}/nspin1"
BOHR = 0.529177210903
R_CUT = 14.5                                                                    # A ; L12/2 = 14.795 A
TOL_MATCH = 1e-4                                                                # A, appariement des sites ideaux
RES = {15: dict(nodes=1, ntasks_per_node=192, time="12:00:00"),
       18: dict(nodes=1, ntasks_per_node=192, time="1-00:00:00"),
       21: dict(nodes=2, ntasks_per_node=96, time="2-00:00:00"),
       24: dict(nodes=3, ntasks_per_node=64, time="3-00:00:00"),
       27: dict(nodes=6, ntasks_per_node=32, time="6-00:00:00")}
sys.path.insert(0, WORK)
import make_inputs_r7 as g                                                      # positions(N, vac), vacancy_A(N), read()


def parse_scf(path):
    txt = g.read(path)
    head, rest = txt.split("CELL_PARAMETERS bohr")
    cell_txt, pos_txt = rest.split("ATOMIC_POSITIONS crystal")
    pos_txt, kp = pos_txt.split("K_POINTS")
    A = np.array([[float(v) for v in l.split()] for l in cell_txt.strip().splitlines()])          # lignes = a_i (bohr)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in pos_txt.strip().splitlines()])
    return head, cell_txt, X, A, kp


def final_coords(relax_out):
    txt = g.read(relax_out)
    blk = re.search(r"Begin final coordinates(.*?)End final coordinates", txt, re.S).group(1)
    blk = blk.split("ATOMIC_POSITIONS (crystal)")[1]
    return np.array([[float(v) for v in l.split()[1:4]] for l in blk.strip().splitlines() if re.match(r"\s*C\b", l)])


def minimg_cart(dx_red, A_A):
    """dx_red (n,3) reduits -> image minimale, cartesien (A). A_A lignes = a_i en A."""
    s = dx_red - np.rint(dx_red)
    return s @ A_A


def field_12x12():
    """Champ de deplacements de la 12x12 R2 nspin1 : r_i (cart. A, relatif a la lacune, image minimale), u_i (cart. A)."""
    head, cell_txt, X0, A_b, kp = parse_scf(os.path.join(SC, "12x12", "defective", "scf.in"))
    Xf = final_coords(os.path.join(SER, "12x12", "nspin1", "relax.out"))
    assert len(Xf) == len(X0) == 287, (len(Xf), len(X0))
    A_A = A_b * BOHR
    vac = (5, 5, "B"); s12 = np.array([(5 + 2 / 3) / 12, (5 + 2 / 3) / 12, 0.0])
    r = minimg_cart(X0 - s12, A_A); u = minimg_cart(Xf - X0, A_A)
    # unicite de l'image minimale : distance a la 2e image
    ok_unique = True
    for i in range(len(X0)):
        d = X0[i] - s12; cands = [np.linalg.norm((d - np.rint(d) + np.array([a, b, 0])) @ A_A) for a in (-1, 0, 1) for b in (-1, 0, 1)]
        cands.sort()
        if cands[1] - cands[0] < 1e-6:
            ok_unique = False
    return dict(X0=X0, Xf=Xf, A_A=A_A, A_b=A_b, s12=s12, r=r, u=u, unique=ok_unique, head=head, cell_txt=cell_txt, kp=kp)


def profile(F, lines):
    rn = np.linalg.norm(F["r"], axis=1); un = np.linalg.norm(F["u"], axis=1); uz = np.abs(F["u"][:, 2])
    lines.append(f"champ 12x12 R2 nspin1 : 287 atomes ; |u| max {un.max():.5f} A (atome {un.argmax()+1}, |r| {rn[un.argmax()]:.3f} A) ; |u_z| max {uz.max():.1e} A ; image minimale unique : {F['unique']}")
    for lo, hi in ((0, 2), (2, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, R_CUT), (R_CUT, 20)):
        m = (rn >= lo) & (rn < hi)
        if m.any():
            lines.append(f"  {lo:5.1f} <= |r| < {hi:5.1f} A : {m.sum():3d} atomes, |u| max {un[m].max():.5f} A, moyen {un[m].mean():.5f} A")
    beyond = rn >= R_CUT
    lines.append(f"au-dela de R_CUT = {R_CUT} A : {beyond.sum()} atomes non transplantes, |u| max {un[beyond].max() if beyond.any() else 0:.5f} A (discontinuite du champ a la troncature)")
    # trois voisins de la lacune (indices 1-based) et leurs distances finales
    nn = np.argsort(rn)[:3]; lines.append(f"voisins de la lacune 12x12 (1-based) : {[int(i)+1 for i in nn]} ; |r| ideal {rn[nn].round(5)} A")
    Xc = (F["X0"] + 0) @ F["A_A"]; Xfc = F["Xf"] @ F["A_A"]
    def dfin(i, j):
        d = F["Xf"][i] - F["Xf"][j]; return np.linalg.norm((d - np.rint(d)) @ F["A_A"])
    lines.append(f"distances finales entre voisins : {dfin(nn[0], nn[1]):.5f}, {dfin(nn[0], nn[2]):.5f}, {dfin(nn[1], nn[2]):.5f} A")
    return rn, un


def field_relaxed(Ns):
    """Champ de deplacements d'une cellule Ns x Ns relaxee de la serie R7c (lacune A, i = j = Ns//2) : r_i, u_i cartesiens (A)."""
    head, cell_txt, X0, A_b, kp = parse_scf(os.path.join(SC, f"{Ns}x{Ns}", "defective", "scf.in"))
    Xf = final_coords(os.path.join(SER, f"{Ns}x{Ns}", "nspin1", "relax.out"))
    assert len(Xf) == len(X0) == 2 * Ns * Ns - 1, (len(Xf), len(X0))
    A_A = A_b * BOHR; vac = g.vacancy_A(Ns); s = np.array([(vac[0] + 1 / 3) / Ns] * 2 + [0.0])
    r = minimg_cart(X0 - s, A_A); u = minimg_cart(Xf - X0, A_A)
    return dict(X0=X0, Xf=Xf, A_A=A_A, A_b=A_b, s12=s, r=r, u=u, unique=True, head=head, cell_txt=cell_txt, kp=kp, Ns=Ns, vac=vac,
                rcut_max=float(np.linalg.norm(A_A[0]) / 2))                          # L/2 : unicite de l'image minimale


def transplant(F, N, invert=True, rcut=R_CUT):
    """Positions reduites de la cellule N x N (lacune A, i=j=N//2) avec le champ transplante. Retourne (X_red, info)."""
    vac = g.vacancy_A(N) if (invert or F.get("Ns")) else (5, 5, "B")
    i0 = vac[0]; f = 1 / 3 if vac[2] == "A" else 2 / 3
    sN = np.array([(i0 + f) / N, (i0 + f) / N, 0.0])
    XN = np.array(g.positions(N, vac)); assert len(XN) == 2 * N * N - 1
    A_A = F["A_A"] * (N / float(F.get("Ns", 12)))                                # cellule source x N/N_source (comme le generateur R7)
    rN = minimg_cart(XN - sN, A_A)
    sign = -1.0 if invert else 1.0
    src = sign * F["r"]; usrc = sign * F["u"]                                  # inversion : r -> -r, u -> -u
    rn12 = np.linalg.norm(F["r"], axis=1); sel = np.where(rn12 < rcut)[0]
    u_N = np.zeros_like(rN); matched = np.zeros(len(XN), bool); used = 0
    for i in sel:
        d = np.linalg.norm(rN - src[i], axis=1); j = int(np.argmin(d))
        if d[j] > TOL_MATCH:
            raise RuntimeError(f"N={N}: site 12x12 {i+1} (|r| {rn12[i]:.3f} A) sans partenaire (min {d[j]:.2e} A)")
        if matched[j]:
            raise RuntimeError(f"N={N}: site N x N {j+1} apparie deux fois")
        matched[j] = True; u_N[j] = usrc[i]; used += 1
    Xc = XN @ A_A + u_N
    Xred = np.mod(Xc @ np.linalg.inv(A_A), 1.0)
    # controles geometriques
    rn = np.linalg.norm(rN, axis=1); nn = np.argsort(rn)[:3]
    def dist(i, j):
        d = Xred[i] - Xred[j]; return np.linalg.norm((d - np.rint(d)) @ A_A)
    dmin = np.inf
    for i in range(len(Xred)):
        d = Xred - Xred[i]; d -= np.rint(d); dd = np.linalg.norm(d @ A_A, axis=1); dd[i] = np.inf; dmin = min(dmin, dd.min())
    info = dict(N=N, vac=vac, sN=sN.tolist(), n_transplanted=used, n_sel=len(sel), nn_1based=[int(i) + 1 for i in nn],
                d_nn=[dist(nn[0], nn[1]), dist(nn[0], nn[2]), dist(nn[1], nn[2])], dmin=dmin, umax=float(np.linalg.norm(u_N, axis=1).max()),
                n_moved=int((np.linalg.norm(u_N, axis=1) > 0).sum()))
    return Xred, info


def build_head(head, N, title):
    """Derivation R2 (build_head de make_inputs_series.py, nspin 1), titre libre."""
    out = []; lines = head.splitlines(); k = 0
    while k < len(lines) and lines[k].startswith("!"):
        k += 1
    out += title; block = None
    for l in lines[k:]:
        s = l.strip()
        if s.startswith("&"):
            block = s.upper()
        if re.match(r"\s*calculation\s*=", l):
            l = "  calculation = 'relax'     ! ions relaxés, cellule fixe"
        elif re.match(r"\s*prefix\s*=", l):
            l = f"  prefix      = 'vac_{N}x{N}_relax_nspin1'"
        elif re.match(r"\s*outdir\s*=", l):
            l = f"  outdir      = '{SCRATCH.format(N=N)}'"
        elif re.match(r"\s*pseudo_dir\s*=", l):
            out += ["  forc_conv_thr = 1.0d-4", "  etot_conv_thr = 1.0d-5", "  nstep       = 200", "  tprnfor     = .true."]
        out.append(l)
        if re.match(r"\s*assume_isolated\s*=", l):
            out += ["  nosym       = .true.", "  noinv       = .true."]
        if s == "/" and block == "&ELECTRONS":
            out += ["", "&IONS", "  ion_dynamics = 'bfgs'", "/"]; block = None
    return "\n".join(out) + "\n"


def build_relax_in(N, Xred, title):
    head, cell_txt, X, A, kp = parse_scf(os.path.join(SC, f"{N}x{N}", "defective", "scf.in"))
    assert len(X) == len(Xred)
    lines = [f"  C  {x:.10f}  {y:.10f}  {z:.10f}" for x, y, z in Xred]
    return build_head(head, N, title) + "\nCELL_PARAMETERS bohr" + cell_txt + "ATOMIC_POSITIONS crystal\n" + "\n".join(lines) + "\n\nK_POINTS" + kp


def build_submit(N):
    t = g.read(os.path.join(SER, "12x12", "nspin1", "submit.relax")); r = RES[N]
    t = t.replace("#SBATCH --job-name=R2_relax_12x12_nspin1", f"#SBATCH --job-name=R7c_relax_{N}x{N}_nspin1")
    t = t.replace("#SBATCH --ntasks=64\n", f"#SBATCH --nodes={r['nodes']}\n#SBATCH --ntasks={r['nodes'] * r['ntasks_per_node']}\n#SBATCH --ntasks-per-node={r['ntasks_per_node']}\n")
    t = re.sub(r"#SBATCH --mem-per-cpu=\S+", "#SBATCH --mem=0", t)
    t = re.sub(r"#SBATCH --time=\S+", f"#SBATCH --time={r['time']}", t)
    t = t.replace("mkdir -p /home/gregb26/links/scratch/qe_tmp/vacancy_relaxed/series/12x12/nspin1", "mkdir -p " + SCRATCH.format(N=N))
    return t


def title_for(N, info, F=None, rcut=R_CUT):
    if F is not None and F.get("Ns"):
        src = f"champ de déplacements de la {F['Ns']}x{F['Ns']} relaxée (R7c nspin1) transplanté sans inversion (lacune A -> A)"
    else:
        src = "champ de déplacements de la 12x12 relaxée (R2 nspin1) transplanté par inversion (lacune B -> A)"
    return [f"! R7c — relax (cellule fixe) de la supercellule {N}x{N} avec lacune, nspin=1",
            f"! dérivé de super_cell/{N}x{N}/defective/scf.in (R7) ; géométrie initiale = {src}, R_cut = {rcut} A, {info['n_transplanted']} atomes déplacés ; voisins {info['nn_1based']}"]


def check(F):
    ok = True
    # (1) auto-transplant 12x12 sans inversion, tous les atomes
    Xred, info = transplant(F, 12, invert=False, rcut=1e9)
    d = Xred - np.mod(F["Xf"], 1.0); d -= np.rint(d); dev = np.abs(d).max()
    print(f"[check 1] auto-transplant 12x12 (B, sans inversion, tous atomes) = coordonnees finales R2 : ecart max {dev:.1e} (reduit) ; {info['n_transplanted']} atomes"); ok &= dev < 1e-9
    # (2) namelists : derivation du scf.in 12x12 = relax.in R2 hors titre/positions
    r2 = g.read(os.path.join(SER, "12x12", "nspin1", "relax.in"))
    r2_title = [l for l in r2.splitlines()[:2]]
    mine = build_relax_in(12, np.mod(F["Xf"], 1.0), r2_title)
    strip = lambda t: "\n".join(l for l in t.splitlines() if not l.startswith("  C  ") and not l.startswith("!"))
    same = strip(mine) == strip(r2)
    print(f"[check 2] namelists/cellule derives du scf.in 12x12 == relax.in R2 12x12 nspin1 (hors titre et positions) : {same}"); ok &= same
    if not same:
        import difflib; sys.stdout.writelines(difflib.unified_diff(strip(r2).splitlines(True), strip(mine).splitlines(True), "R2", "R7c"))
    # (3) inversion : maille ideale B-lacune 12x12 -> maille ideale A-lacune 12x12 (ensembles de vecteurs relatifs)
    vacA = g.vacancy_A(12); XA = np.array(g.positions(12, vacA)); sA = np.array([(vacA[0] + 1 / 3) / 12] * 2 + [0.0])
    rA = minimg_cart(XA - sA, F["A_A"]); rB_inv = -F["r"]
    sel = np.linalg.norm(rB_inv, axis=1) < R_CUT; nmatch = 0
    for v in rB_inv[sel]:
        if np.linalg.norm(rA - v, axis=1).min() < TOL_MATCH: nmatch += 1
    print(f"[check 3] inversion B -> A sur la 12x12 ideale : {nmatch}/{sel.sum()} sites apparies a {TOL_MATCH} A"); ok &= nmatch == sel.sum()
    # (4) profil du champ
    lines = []; profile(F, lines); print("\n".join("[profil] " + l for l in lines))
    print("[check] REGRESSION", "PASS" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true"); ap.add_argument("--write", action="store_true"); ap.add_argument("--sizes", default="15,18,21,24,27")
    ap.add_argument("--source", type=int, default=12, help="12 (R2, inversion B->A) ou N relaxee de R7c (A->A, sans inversion)")
    ap.add_argument("--force", action="store_true", help="ecrase relax.in existant (sauvegarde .bak)")
    a = ap.parse_args()
    if a.source == 12:
        F = field_12x12(); rcut = R_CUT; invert = True
        if a.check and not check(F):
            sys.exit(1)
    else:
        F = field_relaxed(a.source); rcut = F["rcut_max"] - 0.3; invert = False
        Xs, info_s = transplant(F, a.source, invert=False, rcut=1e9)
        d = Xs - np.mod(F["Xf"], 1.0); d -= np.rint(d); dev = np.abs(d).max()
        rn = np.linalg.norm(F["r"], axis=1); un = np.linalg.norm(F["u"], axis=1); beyond = rn >= rcut
        print(f"[check source {a.source}x{a.source}] auto-transplant = coordonnees finales : ecart max {dev:.1e} ; R_cut {rcut:.2f} A ; "
              f"{beyond.sum()} atomes au-dela (|u| max {un[beyond].max() if beyond.any() else 0:.5f} A, moyen {un[beyond].mean() if beyond.any() else 0:.5f}) ; |u| max {un.max():.5f} A")
        lines = []; 
        for lo, hi in ((0, 3), (3, 8), (8, 13), (13, 18), (18, 23), (23, 28), (28, 35)):
            m = (rn >= lo) & (rn < hi)
            if m.any(): lines.append(f"{lo}-{hi} A : {m.sum()} at., |u| max {un[m].max():.4f}, moyen {un[m].mean():.4f}")
        print("[profil " + f"{a.source}x{a.source}] " + " ; ".join(lines))
        if a.check and dev > 1e-9:
            print("[check] FAIL"); sys.exit(1)
    for N in [int(s) for s in a.sizes.split(",")]:
        Xred, info = transplant(F, N, invert=invert, rcut=rcut)
        print(f"[{N}x{N}] lacune {info['vac']} s_red ({info['sN'][0]:.6f}, {info['sN'][1]:.6f}) ; transplantes {info['n_transplanted']}/{info['n_sel']} ; voisins {info['nn_1based']} d = "
              f"{[round(v, 5) for v in info['d_nn']]} A ; d_min {info['dmin']:.4f} A ; |u| max {info['umax']:.5f} A ; ressources {RES[N]}")
        if a.write:
            d = os.path.join(SER, f"{N}x{N}", "nspin1"); os.makedirs(d, exist_ok=True)
            for f, txt in (("relax.in", build_relax_in(N, Xred, title_for(N, info, F, rcut))), ("submit.relax", build_submit(N))):
                p = os.path.join(d, f)
                if os.path.exists(p):
                    if a.force and f == "relax.in" and g.read(p) != txt:
                        import shutil; bak = p + f".bak_source{'12' if a.source == 12 else a.source}"; shutil.copy2(p, bak); print(f"[write] sauvegarde {bak}")
                    else:
                        print(f"[write] {p} existe : NON ecrase"); continue
                open(p, "w").write(txt); print(f"[write] {p}")
            srcname = f"{F['Ns']}x{F['Ns']} R7c nspin1 (relax.out final), sans inversion (A -> A)" if F.get("Ns") else "12x12 R2 nspin1 (relax.out final) par inversion (B -> A)"
            lines = [f"R7c {N}x{N} nspin1 : geometrie initiale transplantee de la {srcname} ; R_cut {rcut} A ; TOL {TOL_MATCH} A",
                     f"lacune {info['vac']} s_red {info['sN']} ; {info['n_transplanted']} atomes deplaces sur {2*N*N-1} ({info['n_sel']} sites 12x12 sous R_cut, tous apparies)",
                     f"voisins (1-based) {info['nn_1based']} ; distances entre voisins {[round(v, 5) for v in info['d_nn']]} A (12x12 R2 final : 2.15980 / 2.56934 / 2.56934) ; d_min {info['dmin']:.5f} A ; |u| max {info['umax']:.5f} A"]
            if not F.get("Ns"):
                profile(F, lines)
            open(os.path.join(SER, f"{N}x{N}", "transplant.log"), "a" if F.get("Ns") else "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
