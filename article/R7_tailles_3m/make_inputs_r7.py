#!/usr/bin/env python
"""
make_inputs_r7.py -- R7 (2026-09-25) : inputs QE des super-cellules 15x15 et 18x18 (lacune non relaxée, sans spin), dérivés du 12x12.

Seuls changent : nat, CELL_PARAMETERS (cellule 12x12 x N/12), ATOMIC_POSITIONS (formule des tailles 5...12, verifiee en phase 0 :
A = ((i+1/3)/N, (j+1/3)/N), B = ((i+2/3)/N, (j+2/3)/N), i exterieur, j interieur, A puis B), nbnd = N_occ + 30, prefix, outdir,
ressources SLURM (submit.scf : ntasks-per-node, mem-per-cpu, time), prefix/outdir/filplot de pp.in. submit.pp inchange.
Lacune : sous-reseau A, atome A le plus proche du centre (0.5, 0.5) -> i = j = floor(N/2) (regle des tailles A 5, 7, 8, 9, 11).
Extension 21/24/27 (2026-09-25, GO en attente) : nbnd = N_occ + ceil(30 (N/12)^2) (meme couverture en energie au-dessus de E_D que la
12x12, ~1,96 eV ; +30 ne couvrirait que ~1,1 / 1,0 / 0,9 eV) ; noeuds entiers (--mem=0), 192 rangs = nr3 planes FFT ; pp.x 48G pour 27x27.
Regression : --check regenere 9x9 (A, i=j=4) et 12x12 (B, i=j=5) et exige l'identite byte a byte avec les scf.in existants.
Usage : python make_inputs_r7.py --check          (regression seulement, rien n'est ecrit)
        python make_inputs_r7.py --write          (ecrit super_cell/{15x15,18x18}/{defective,pristine}/{scf.in,submit.scf,pp.in,submit.pp})
        python make_inputs_r7.py --diff           (imprime les diffs contre le 12x12, positions exclues)
"""
import os, re, sys, argparse, difflib

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                     # .../graphene/qe
SC = os.path.join(GQ, "defects", "super_cell")
TEMPLATE_N = 12
SCRATCH = "/home/gregb26/links/scratch/qe_tmp"
NEXTRA = 30                                                                     # nbnd = N_occ + 30 (tailles 6...18 ; 5x5 d : +32)


def nextra(N):
    """Bandes vides ajoutees : 30 jusqu'a 18x18 (regle 5...12) ; ceil(30 (N/12)^2) au-dela (couverture ~1,96 eV au-dessus de E_D)."""
    import math
    return NEXTRA if N <= 18 else int(math.ceil(NEXTRA * (N / 12.0) ** 2))

RES = {15: dict(ntasks_per_node=8, mem_per_cpu="8G", time="03:00:00"),
       18: dict(ntasks_per_node=8, mem_per_cpu="16G", time="06:00:00"),
       21: dict(nodes=2, ntasks_per_node=96, mem="0", time="06:00:00"),          # 192 rangs, noeuds entiers (memoire totale ~0,7 To)
       24: dict(nodes=3, ntasks_per_node=64, mem="0", time="08:00:00"),          # 192 rangs (~1,2 To)
       27: dict(nodes=6, ntasks_per_node=32, mem="0", time="12:00:00")}          # 192 rangs (~1,9 To) ; 6 noeuds : exces du rang 0 mesure 0,78 x wfc a 18x18
PP_MEM = {27: "48G"}                                                            # pp.x : 32G suffit jusqu'a 24x24 (~21 G), 27x27 ~27 G


def positions(N, vac=None):
    """Liste des positions (crystal) de la parfaite N x N ; vac = (i, j, 'A'|'B') retire un atome."""
    P = []
    for i in range(N):
        for j in range(N):
            for s, f in (("A", 1.0 / 3), ("B", 2.0 / 3)):
                if vac is not None and (i, j, s) == tuple(vac):
                    continue
                P.append(((i + f) / N, (j + f) / N, 0.0))
    return P


def vacancy_A(N):
    """Atome A le plus proche du centre : i = j = argmin |(i+1/3)/N - 1/2| (= floor(N/2) pour toutes les tailles A existantes)."""
    i = min(range(N), key=lambda i: abs((i + 1.0 / 3) / N - 0.5))
    assert i == N // 2, (N, i)
    return (i, i, "A")


def read(path):
    with open(path, "r") as f:
        return f.read()


def scf_in(N, kind, vac, template):
    """kind = 'defective' | 'pristine'. Rebuilt from the 12x12 template of the same kind by targeted substitutions."""
    t = template
    nat_p = 2 * N * N; nat = nat_p - (1 if kind == "defective" else 0); nel = 4 * nat; nocc = nel // 2; nbnd = nocc + nextra(N)
    tag = f"{N}x{N}"; prefix = f"defect_{tag}_{'d' if kind == 'defective' else 'p'}"
    t = t.replace(f"{TEMPLATE_N}x{TEMPLATE_N}", tag)
    t = re.sub(r"prefix      = 'defect_\d+x\d+_[dp]'", f"prefix      = '{prefix}'", t)
    t = re.sub(r"outdir      = '[^']*'", f"outdir      = '{SCRATCH}/{prefix}'", t)
    if kind == "defective":
        t = re.sub(r"nat         = \d+\s+! number of atoms \(\d+ = 2 x \d+ x \d+ - 1 vacancy\)",
                   f"nat         = {nat}       ! number of atoms ({nat} = 2 x {N} x {N} - 1 vacancy)", t)
    else:
        t = re.sub(r"nat         = \d+\s+! number of atoms \(\d+ = 2 x \d+ x \d+\)",
                   f"nat         = {nat}       ! number of atoms ({nat} = 2 x {N} x {N})", t)
    t = re.sub(r"nbnd        = \d+\s+! \d+ valence electrons -> \d+ occupied; adjust as needed",
               f"nbnd        = {nbnd}      ! {nel} valence electrons -> {nocc} occupied; adjust as needed", t)
    # cellule : les deux premieres lignes du 12x12 x N/12 (la 9x9 existante = 12x12 x 9/12 exactement)
    m = re.search(r"CELL_PARAMETERS bohr\n(.*)\n(.*)\n(.*)\n", t)
    l1, l2, l3 = m.group(1), m.group(2), m.group(3)
    nums = re.findall(r"-?\d+\.\d+", l1)
    ax, ay = float(nums[0]) * N / TEMPLATE_N, float(nums[1]) * N / TEMPLATE_N
    new1 = l1.replace(nums[0], f"{ax:.10f}").replace(nums[1], f"{ay:.10f}")
    nums2 = re.findall(r"-?\d+\.\d+", l2)
    new2 = l2.replace(nums2[0], f"{float(nums2[0]) * N / TEMPLATE_N:.10f}").replace(nums2[1], f"{float(nums2[1]) * N / TEMPLATE_N:.10f}")
    t = t.replace(f"CELL_PARAMETERS bohr\n{l1}\n{l2}\n{l3}\n", f"CELL_PARAMETERS bohr\n{new1}\n{new2}\n{l3}\n")
    # positions
    pos = positions(N, vac if kind == "defective" else None)
    assert len(pos) == nat
    block = "".join(f"  C  {x:.10f}  {y:.10f}  {z:.10f}\n" for x, y, z in pos)
    t = re.sub(r"(ATOMIC_POSITIONS crystal\n)(?:  C  .*\n)+", lambda mm: mm.group(1) + block, t)
    return t


def submit_scf(N, template):
    r = RES[N]; t = template
    t = re.sub(r"#SBATCH --time=\S+", f"#SBATCH --time={r['time']}", t)
    if "mem" in r:                                                              # noeuds entiers : --nodes, --ntasks, --ntasks-per-node, --mem=0
        t = re.sub(r"#SBATCH --mem-per-cpu=\S+", f"#SBATCH --mem={r['mem']}", t)
        t = t.replace("#SBATCH --ntasks=64\n", f"#SBATCH --nodes={r['nodes']}\n#SBATCH --ntasks={r['nodes'] * r['ntasks_per_node']}\n#SBATCH --ntasks-per-node={r['ntasks_per_node']}\n")
    else:
        t = re.sub(r"#SBATCH --mem-per-cpu=\S+", f"#SBATCH --mem-per-cpu={r['mem_per_cpu']}", t)
        t = t.replace("#SBATCH --ntasks=64\n", f"#SBATCH --ntasks=64\n#SBATCH --ntasks-per-node={r['ntasks_per_node']}\n")
    return t


def submit_pp(N, template):
    return re.sub(r"#SBATCH --mem-per-cpu=\S+", f"#SBATCH --mem-per-cpu={PP_MEM[N]}", template) if N in PP_MEM else template


def pp_in(N, kind, template):
    tag = f"{N}x{N}"; s = "d" if kind == "defective" else "p"; prefix = f"defect_{tag}_{s}"
    t = re.sub(r"prefix   = 'defect_\d+x\d+_[dp]'", f"prefix   = '{prefix}'", template)
    t = re.sub(r"outdir      = '[^']*'", f"outdir      = '{SCRATCH}/{prefix}'", t)
    t = re.sub(r"filplot  = 'Vks_\d+x\d+_[dp]'", f"filplot  = 'Vks_{tag}_{s}'", t)
    return t


def templates(kind):
    d = os.path.join(SC, f"{TEMPLATE_N}x{TEMPLATE_N}", kind)
    return {f: read(os.path.join(d, f)) for f in ("scf.in", "submit.scf", "pp.in", "submit.pp")}


def check():
    ok = True
    for N, vac in ((9, (4, 4, "A")), (12, (5, 5, "B"))):
        for kind in ("defective", "pristine"):
            gen = scf_in(N, kind, vac, templates(kind)["scf.in"]); ref = read(os.path.join(SC, f"{N}x{N}", kind, "scf.in"))
            same = gen == ref; ok &= same
            print(f"[check] {N}x{N} {kind}: scf.in regenere == existant : {same}")
            if not same:
                sys.stdout.writelines(difflib.unified_diff(ref.splitlines(True), gen.splitlines(True), "existant", "regenere"))
        genpp = pp_in(N, "defective", templates("defective")["pp.in"]); refpp = read(os.path.join(SC, f"{N}x{N}", "defective", "pp.in"))
        print(f"[check] {N}x{N} defective: pp.in regenere == existant : {genpp == refpp}")
        ok &= genpp == refpp
    for N in (15, 18):                                                          # tailles deja ecrites : regeneration identique
        vac, out = build(N)
        for kind in ("defective", "pristine"):
            for f in ("scf.in", "submit.scf", "pp.in", "submit.pp"):
                p = os.path.join(SC, f"{N}x{N}", kind, f)
                if os.path.exists(p):
                    same = read(p) == out[kind][f]; ok &= same
                    if not same:
                        print(f"[check] {N}x{N} {kind} {f}: regenere != ecrit")
        print(f"[check] {N}x{N}: 8 fichiers ecrits reproduits a l'identique : {ok}")
    print("[check] REGRESSION", "PASS" if ok else "FAIL")
    return ok


def build(N):
    vac = vacancy_A(N); out = {}
    for kind in ("defective", "pristine"):
        T = templates(kind)
        out[kind] = {"scf.in": scf_in(N, kind, vac, T["scf.in"]), "submit.scf": submit_scf(N, T["submit.scf"]),
                     "pp.in": pp_in(N, kind, T["pp.in"]), "submit.pp": submit_pp(N, T["submit.pp"])}
    return vac, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true"); ap.add_argument("--write", action="store_true"); ap.add_argument("--diff", action="store_true")
    ap.add_argument("--sizes", default="15,18")
    a = ap.parse_args()
    if a.check and not check():
        sys.exit(1)
    for N in [int(s) for s in a.sizes.split(",")]:
        vac, out = build(N)
        i = vac[0]; print(f"[{N}x{N}] lacune A (i,j)=({i},{i}) -> ({(i+1/3)/N:.10f}, {(i+1/3)/N:.10f}), indice 1-based dans la parfaite {2*(i*N+i)+1} ; "
                          f"nat p/d {2*N*N}/{2*N*N-1} ; nbnd p/d {2*N*N*2+nextra(N)}/{(2*N*N-1)*2+nextra(N)} (+{nextra(N)}) ; ressources {RES[N]}")
        for kind in ("defective", "pristine"):
            d = os.path.join(SC, f"{N}x{N}", kind)
            if a.diff:
                T = templates(kind)
                for f in ("scf.in", "submit.scf", "pp.in", "submit.pp"):
                    ref = [l for l in T[f].splitlines(True) if not l.startswith("  C  ")]; gen = [l for l in out[kind][f].splitlines(True) if not l.startswith("  C  ")]
                    dd = list(difflib.unified_diff(ref, gen, f"12x12/{kind}/{f}", f"{N}x{N}/{kind}/{f}", n=0))
                    print("".join(dd) if dd else f"--- {N}x{N}/{kind}/{f}: identique au 12x12\n", end="")
                print(f"    ({N}x{N}/{kind}/scf.in : {sum(1 for l in out[kind]['scf.in'].splitlines() if l.startswith('  C  '))} lignes de positions, non montrees)")
            if a.write:
                os.makedirs(d, exist_ok=True)
                for f, txt in out[kind].items():
                    p = os.path.join(d, f)
                    if os.path.exists(p) and read(p) != txt:
                        print(f"[write] {p} existe et differe : NON ecrase (supprimer a la main pour regenerer)"); continue
                    with open(p, "w") as fh:
                        fh.write(txt)
                    print(f"[write] {p}")


if __name__ == "__main__":
    main()
