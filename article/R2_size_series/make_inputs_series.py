#!/usr/bin/env python3
"""R2 — série en taille de la lacune relaxée. Dérivé de R1 make_inputs.py (article/R1_vacancy_relaxed/), généralisé :
détection de la lacune et de ses trois voisins par distance à image minimale (quel que soit le sous-réseau),
paire perturbée = les deux voisins d'indice le plus bas (rapprochés de 0.03 Å chacun, dans le plan), le troisième reste.
Pour chaque N : series/NxN/{nspin1,nspin2}/relax.in + submit.relax, series/NxN/perturbation.log.
Protocole identique à R1 : calculation='relax' (cellule fixe, bfgs, nstep 200, forc 1e-4, etot 1e-5, tprnfor),
nosym/noinv, nspin2 = espèce C1 (même pseudo) sur les trois voisins, starting_magnetization(C1)=0.5.
Tout le reste du scf.in du ch. 4 (cutoffs, Γ, mv 0.01 Ry, conv_thr 1e-10, nbnd, pseudo) inchangé.
Usage : python3 make_inputs_series.py [N ...]   (défaut 5 6 7 8 10 11 12 ; 9 = détection seule, R1 déjà fait)"""
import os, re, sys, json
import numpy as np

ROOT = os.environ.get("PROJECTS", "/home/gregb26/links/projects/rrg-cotemich-ac/gregb26")
SRC = ROOT + "/graphene/qe/defects/super_cell/{N}x{N}/defective/scf.in"
SIDECAR = ROOT + "/graphene-raman/results/M/M_ed_{N}x{N}.json"      # convention A/B du dépôt (tag_vacancy_sublattice.py)
SCRATCH = "/home/gregb26/links/scratch/qe_tmp/vacancy_relaxed/series/{N}x{N}/nspin{n}"
HERE = os.path.dirname(os.path.abspath(__file__))
BOHR = 0.529177210903
DISP_A = 0.03            # Å par atome
A_CC_TOL = 1e-3          # Å
SIZES = [int(a) for a in sys.argv[1:]] or [5, 6, 7, 8, 10, 11, 12]

def minimg(d, A, Ainv):
    """vecteurs cartésiens d -> image minimale dans la cellule A (lignes = vecteurs)."""
    s = d @ Ainv
    s -= np.round(s)
    return s @ A

def parse(path):
    txt = open(path).read()
    head, rest = txt.split("CELL_PARAMETERS bohr")
    cell_txt, pos_txt = rest.split("ATOMIC_POSITIONS crystal")
    pos_txt, kp = pos_txt.split("K_POINTS")
    A = np.array([[float(v) for v in l.split()] for l in cell_txt.strip().splitlines()])
    pos = [l.split() for l in pos_txt.strip().splitlines()]
    X = np.array([[float(v) for v in p[1:4]] for p in pos])
    return head, cell_txt, pos, X, A, kp

def detect(N, X, A):
    """lacune + 3 voisins. Retourne dict (indices 1-based)."""
    Ainv = np.linalg.inv(A)
    R = X @ A
    nat = len(R)
    assert nat == 2 * N * N - 1, (N, nat)
    # distances à image minimale
    D = np.zeros((nat, nat))
    for i in range(nat):
        d = minimg(R - R[i], A, Ainv)
        D[i] = np.linalg.norm(d, axis=1)
    np.fill_diagonal(D, np.inf)
    a_cc = D.min() * BOHR                                  # Å
    tol = A_CC_TOL / BOHR
    nn = [np.where(np.abs(D[i] - a_cc / BOHR) < tol)[0] for i in range(nat)]
    two = [i for i in range(nat) if len(nn[i]) == 2]
    assert all(len(nn[i]) in (2, 3) for i in range(nat)), "coordination inattendue"
    assert len(two) == 3, f"{len(two)} atomes bicoordonnés (3 attendus)"
    # position de la lacune : pour chaque voisin i, liaison manquante = -(v_j + v_k)
    est = []
    for i in two:
        v = minimg(R[nn[i]] - R[i], A, Ainv)
        est.append(R[i] - v.sum(axis=0))
    est = np.array(est)
    est_rel = minimg(est - est[0], A, Ainv) + est[0]
    assert np.abs(est_rel - est_rel.mean(axis=0)).max() * BOHR < 1e-4, "estimations de la lacune incohérentes"
    r_vac = est_rel.mean(axis=0)
    d_vac = np.linalg.norm(minimg(R[two] - r_vac, A, Ainv), axis=1) * BOHR
    assert np.all(np.abs(d_vac - a_cc) < A_CC_TOL), (d_vac, a_cc)
    s_vac = np.mod(r_vac @ Ainv, 1.0)
    s_uc = np.mod(N * s_vac, 1.0)
    sub = "A" if np.allclose(s_uc[:2], 1/3, atol=1e-3) else ("B" if np.allclose(s_uc[:2], 2/3, atol=1e-3) else "?")
    neigh = tuple(sorted(i + 1 for i in two))
    return dict(a_cc=a_cc, r_vac=r_vac, s_vac=s_vac, sub=sub, neigh=neigh, pair=neigh[:2], third=neigh[2], d_vac=d_vac)

def perturb(X, A, pair, log):
    Ainv = np.linalg.inv(A)
    i, j = pair[0] - 1, pair[1] - 1
    ri, rj = X[i] @ A, X[j] @ A
    dij = minimg((rj - ri)[None, :], A, Ainv)[0]
    u = dij / np.linalg.norm(dij); u[2] = 0.0
    dd = DISP_A / BOHR
    Xp = X.copy()
    Xp[i] = (ri + dd * u) @ Ainv
    Xp[j] = (rj - dd * u) @ Ainv
    log.append(f"d({pair[0]}-{pair[1]}) idéal = {np.linalg.norm(dij):.6f} bohr = {np.linalg.norm(dij)*BOHR:.6f} A")
    log.append(f"u ({pair[0]}->{pair[1]}) = {u}")
    for k, name in ((i, pair[0]), (j, pair[1])):
        dc = (Xp[k] - X[k]) @ A
        log.append(f"atome {name}: cart. {X[k]@A} -> {Xp[k]@A} bohr ; delta cart = {dc} bohr = {dc*BOHR} A (|d|={np.linalg.norm(dc)*BOHR:.5f} A) ; delta cryst = {Xp[k]-X[k]}")
    dp = minimg(((Xp[j] - Xp[i]) @ A)[None, :], A, Ainv)[0]
    log.append(f"d({pair[0]}-{pair[1]}) perturbé = {np.linalg.norm(dp)*BOHR:.6f} A")
    return Xp

def build_head(head, N, nspin, neigh):
    out = []
    lines = head.splitlines()
    # titre : 2 premières lignes de commentaire
    title = [f"! R2 — relax (cellule fixe) de la supercellule {N}x{N} avec lacune, nspin={nspin}",
             f"! dérivé de super_cell/{N}x{N}/defective/scf.in (ch. 4) ; perturbation initiale voisins {neigh[0]}-{neigh[1]}, C1 = {neigh}"]
    k = 0
    while k < len(lines) and lines[k].startswith("!"):
        k += 1
    out += title
    block = None
    for l in lines[k:]:
        s = l.strip()
        if s.startswith("&"):
            block = s.upper()
        if re.match(r"\s*calculation\s*=", l):
            l = "  calculation = 'relax'     ! ions relaxés, cellule fixe"
        elif re.match(r"\s*prefix\s*=", l):
            l = f"  prefix      = 'vac_{N}x{N}_relax_nspin{nspin}'"
        elif re.match(r"\s*outdir\s*=", l):
            l = f"  outdir      = '{SCRATCH.format(N=N, n=nspin)}'"
        elif re.match(r"\s*pseudo_dir\s*=", l):
            out += ["  forc_conv_thr = 1.0d-4", "  etot_conv_thr = 1.0d-5", "  nstep       = 200", "  tprnfor     = .true."]
        elif re.match(r"\s*ntyp\s*=", l) and nspin == 2:
            l = f"  ntyp        = 2        ! C ({2*N*N-4} atomes) + C1 (3 voisins de la lacune, même pseudo)"
        out.append(l)
        if re.match(r"\s*assume_isolated\s*=", l):
            out += ["  nosym       = .true.", "  noinv       = .true."]
            if nspin == 2:
                out += ["  nspin       = 2", "  starting_magnetization(1) = 0.0   ! C",
                        f"  starting_magnetization(2) = 0.5   ! C1 (atomes {neigh[0]}, {neigh[1]}, {neigh[2]})"]
        if s == "/" and block == "&ELECTRONS":
            out += ["", "&IONS", "  ion_dynamics = 'bfgs'", "/"]
            block = None
        if re.match(r"\s*C\s+12\.011\s+C\.upf", l) and nspin == 2:
            out.append(f"  C1 12.011 C.upf ! voisins de la lacune ({neigh[0]}, {neigh[1]}, {neigh[2]}), même pseudo")
    return "\n".join(out) + "\n"

def build(head, cell_txt, pos, Xp, kp, N, nspin, neigh):
    h = build_head(head, N, nspin, neigh)
    lines = []
    for k in range(len(pos)):
        sp = "C1" if (nspin == 2 and k + 1 in neigh) else "C"
        lines.append(f"  {sp:<2} {Xp[k,0]:.10f}  {Xp[k,1]:.10f}  {Xp[k,2]:.10f}")
    return h + "\nCELL_PARAMETERS bohr" + cell_txt + "ATOMIC_POSITIONS crystal\n" + "\n".join(lines) + "\n\nK_POINTS" + kp

SUBMIT = """#!/bin/bash
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=R2_relax_{N}x{N}_nspin{n}
#SBATCH --mail-type=NONE          # Mail events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G               # memory; default unit is megabytes
#SBATCH --time={walltime}          # time (DD-HH:MM) — BFGS de lacune, comme R1
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module restore qe

export OMP_NUM_THREADS=1

mkdir -p {outdir}
srun pw.x < relax.in > relax.out
"""

rows = []
for N in SIZES:
    head, cell_txt, pos, X, A, kp = parse(SRC.format(N=N))
    det = detect(N, X, A)
    # contrôle croisé avec la convention du dépôt
    side = json.load(open(SIDECAR.format(N=N)))
    ok_sub = side.get("vacancy_sublattice") == det["sub"]
    ok_pos = np.allclose(np.mod(np.array(side.get("vacancy_s_red")) - det["s_vac"] + 0.5, 1) - 0.5, 0, atol=1e-6)
    nbnd = int(re.search(r"nbnd\s*=\s*(\d+)", head).group(1))
    rows.append((N, det["sub"], det["s_vac"], det["neigh"], det["pair"], det["third"], nbnd, det["a_cc"], ok_sub and ok_pos))
    if N == 9:
        continue                     # R1 : détection seule (contrôle), pas de réécriture
    log = [f"R2 {N}x{N} : a_CC = {det['a_cc']:.6f} A ; lacune s_red = {det['s_vac']} ; sous-réseau {det['sub']} ; voisins {det['neigh']} (d = {det['d_vac']} A) ; paire {det['pair']}, troisième {det['third']}"]
    Xp = perturb(X, A, det["pair"], log)
    d = os.path.join(HERE, f"{N}x{N}")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "perturbation.log"), "w").write("\n".join(log) + "\n")
    wall = "03-00:00" if N >= 10 else "01-00:00"
    for n in (1, 2):
        sd = os.path.join(d, f"nspin{n}")
        os.makedirs(sd, exist_ok=True)
        open(os.path.join(sd, "relax.in"), "w").write(build(head, cell_txt, pos, Xp, kp, N, n, det["neigh"]))
        open(os.path.join(sd, "submit.relax"), "w").write(SUBMIT.format(N=N, n=n, walltime=wall, outdir=SCRATCH.format(N=N, n=n)))

print(f"{'N':>3} {'sous-res':>8} {'lacune s_red (x,y)':>22} {'voisins':>16} {'paire':>10} {'3e':>4} {'nbnd':>5} {'a_CC (A)':>9} {'sidecar':>8}")
for N, sub, s, neigh, pair, third, nbnd, acc, ok in rows:
    print(f"{N:>3} {sub:>8} {s[0]:>10.6f} {s[1]:>10.6f} {str(neigh):>16} {str(pair):>10} {third:>4} {nbnd:>5} {acc:>9.5f} {'OK' if ok else 'ECART':>8}")
