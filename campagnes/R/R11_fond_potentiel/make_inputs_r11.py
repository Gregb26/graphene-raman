#!/usr/bin/env python
"""
make_inputs_r11.py -- R11 (2026-10-08) : entrees QE de la paire 9x9 (parfaite, lacune) sous des variantes d'echantillonnage de la
zone de Brillouin, d'elargissement et de boite de vide, derivees des entrees de production de la 9x9 par substitutions ciblees
(meme principe que make_inputs_r7.py : rien n'est reconstruit, le fichier de production est le gabarit).

Gabarits (sur rorqual) : $GQ/defects/super_cell/9x9/{pristine,defective}/{scf.in,pp.in}  (GQ = graphene/qe, deduit de ce fichier :
WORK = $GQ/defects/R11_fond_potentiel ; --templates <dir> pour pointer ailleurs, ex. la copie locale du depot du memoire).
Ne changent que : prefix, outdir, (smearing, degauss), (assume_isolated), (c de CELL_PARAMETERS), (K_POINTS) ; positions, cutoff,
nbnd, conv_thr, mixing, pseudo : identiques a la production. Les lanceurs submit.scf / submit.pp sont ecrits par ce script.

Variantes (dict VARIANTS ; cle = nom du repertoire sous variants/) :
  phase 1 (Gamma, sigma)      g_mv0.005  g_mv0.02  g_mv0.03  g_fd0.005         (production = mv 0.01, non refaite)
  phase 2 (grille k)          k3_mv0.01  k4_mv0.01  k3_mv0.005  k3_fd0.005     (grilles Gamma-centrees, K reste inclus pour N = 3m)
  phase 2b (vide)             g_c40  g_notrunc  g_c40_notrunc                  (c = 40 bohr ; sans assume_isolated='2D')
  pp11 (phase 1)              pp.x plot_num = 11 (V_bare + V_H, sans V_xc) sur les .save de production d'une taille donnee.

Usage (depuis WORK) :
  python make_inputs_r11.py --check [--templates DIR]       porte : la variante identite reproduit les gabarits byte a byte, chaque
                                                            motif de substitution est trouve exactement une fois ; rien n'est ecrit
  python make_inputs_r11.py --diff  [--variants a,b]        diffs contre les gabarits (positions exclues), rien n'est ecrit
  python make_inputs_r11.py --write [--variants a,b]        ecrit variants/<nom>/{pristine,defective}/{scf.in,submit.scf,pp.in,submit.pp}
  python make_inputs_r11.py --pp11 9x9[,6x6,...]            ecrit pp11/<taille>/{pristine,defective}/{pp11.in,submit.pp11}
Un fichier existant et different n'est jamais ecrase (message ; supprimer a la main pour regenerer).
"""
import argparse
import difflib
import os
import re
import sys

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                     # .../graphene/qe  (repertoire de travail rorqual)
SCRATCH = "/home/gregb26/links/scratch/qe_tmp"                                  # outdir des .in de production (meme racine)
SIZE = "9x9"
ACCOUNT = "rrg-cotemich-ac"

# nom -> (K_POINTS, smearing, degauss Ry, assume_isolated 2D ?, c en bohr, limite de temps scf, description)
VARIANTS = {
    # phase 1 : diagnostic sigma au seul point Gamma (attendu : presque insensible, TODO 2026-10-03)
    "g_mv0.005": dict(k=None, smear="mv", degauss=0.005, iso2d=True, c=30.0, time="01:00:00", desc="Gamma, m-v 0.005 Ry"),
    "g_mv0.02": dict(k=None, smear="mv", degauss=0.02, iso2d=True, c=30.0, time="01:00:00", desc="Gamma, m-v 0.02 Ry"),
    "g_mv0.03": dict(k=None, smear="mv", degauss=0.03, iso2d=True, c=30.0, time="01:00:00", desc="Gamma, m-v 0.03 Ry"),
    "g_fd0.005": dict(k=None, smear="fd", degauss=0.005, iso2d=True, c=30.0, time="01:00:00", desc="Gamma, Fermi-Dirac 0.005 Ry"),
    # phase 2 : echantillonnage k de la super-cellule (candidat principal)
    "k3_mv0.01": dict(k=(3, 3, 1), smear="mv", degauss=0.01, iso2d=True, c=30.0, time="03:00:00", desc="3x3x1 Gamma-centree, m-v 0.01 Ry"),
    "k4_mv0.01": dict(k=(4, 4, 1), smear="mv", degauss=0.01, iso2d=True, c=30.0, time="04:00:00", desc="4x4x1 Gamma-centree, m-v 0.01 Ry"),
    "k3_mv0.005": dict(k=(3, 3, 1), smear="mv", degauss=0.005, iso2d=True, c=30.0, time="03:00:00", desc="3x3x1, m-v 0.005 Ry"),
    "k3_fd0.005": dict(k=(3, 3, 1), smear="fd", degauss=0.005, iso2d=True, c=30.0, time="03:00:00", desc="3x3x1, Fermi-Dirac 0.005 Ry"),
    # phase 2b : boite de vide et troncature de Coulomb
    "g_c40": dict(k=None, smear="mv", degauss=0.01, iso2d=True, c=40.0, time="02:00:00", desc="Gamma, c = 40 bohr, 2D"),
    "g_notrunc": dict(k=None, smear="mv", degauss=0.01, iso2d=False, c=30.0, time="02:00:00", desc="Gamma, c = 30 bohr, sans troncature 2D"),
    "g_c40_notrunc": dict(k=None, smear="mv", degauss=0.01, iso2d=False, c=40.0, time="02:00:00", desc="Gamma, c = 40 bohr, sans troncature 2D"),
}
IDENTITY = dict(k=None, smear="mv", degauss=0.01, iso2d=True, c=30.0, time="01:00:00", desc="production (porte --check)")
PHASES = {"1": ["g_mv0.005", "g_mv0.02", "g_mv0.03", "g_fd0.005"], "2": ["k3_mv0.01", "k4_mv0.01", "k3_mv0.005", "k3_fd0.005"],
          "2b": ["g_c40", "g_notrunc", "g_c40_notrunc"]}
SMEAR_COMMENT = {"mv": "Marzari-Vanderbilt cold smearing", "fd": "Fermi-Dirac smearing"}

SUBMIT_SCF = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --job-name=r11_{name}_{s}
#SBATCH --mail-type=NONE
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={mem}
#SBATCH --time={time}
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

module restore qe

export OMP_NUM_THREADS=1

mkdir -p {outdir}
srun pw.x < scf.in > scf.out
"""
SUBMIT_PP = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --job-name=pp_{tag}_{s}
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=00-01:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

module restore qe

export OMP_NUM_THREADS=1

srun pp.x < {ppin} > {ppout}
"""
PP11_IN = """&inputpp
  prefix   = '{prefix}'
  outdir      = '{outdir}'   ! directory of tempory files
  plot_num = 11               ! V_bare + V_Hartree (no V_xc) -- R11, reference du vide sur l'electrostatique seule
  filplot  = 'Vh_{tag}_{s}'
/
&plot
  iflag         = 3             ! 3D
  output_format = 6             ! cube
/
"""


def read(path):
    with open(path, "r") as f:
        return f.read()


def sub1(text, pattern, repl, what):
    """Substitution d'exactement une occurrence (sinon le gabarit n'est pas celui attendu : on refuse)."""
    new, n = re.subn(pattern, repl, text)
    if n != 1:
        raise SystemExit(f"[gabarit] motif « {what} » trouve {n} fois (attendu 1) : {pattern}")
    return new


def scf_in(name, v, kind, template):
    s = "d" if kind == "defective" else "p"
    prefix = f"r11_{name}_{s}" if name != "identity" else f"defect_{SIZE}_{s}"
    outdir = f"{SCRATCH}/R11/{prefix}" if name != "identity" else f"{SCRATCH}/{prefix}"
    t = template
    t = sub1(t, r"prefix      = 'defect_\d+x\d+_[dp]'", f"prefix      = '{prefix}'", "prefix")
    t = sub1(t, r"outdir      = '[^']*'", f"outdir      = '{outdir}'", "outdir")
    t = sub1(t, r"smearing    = '\w+'\s+! [^\n]*", f"smearing    = '{v['smear']}'     ! {SMEAR_COMMENT[v['smear']]}", "smearing")
    t = sub1(t, r"degauss     = [0-9.]+\s+! value of gaussian smearing in Ry", f"degauss     = {v['degauss']:g}     ! value of gaussian smearing in Ry", "degauss")
    if v["iso2d"]:
        sub1(t, r"\n  assume_isolated = '2D' ! helps for 2d systems", "", "assume_isolated (presence)")
    else:
        t = sub1(t, r"\n  assume_isolated = '2D' ! helps for 2d systems", "", "assume_isolated")
    t = sub1(t, r"(CELL_PARAMETERS bohr\n[^\n]*\n[^\n]*\n   0\.0000000000             0\.0000000000            )30\.0000000000\n",
             lambda m: m.group(1) + f"{v['c']:.10f}\n", "c de CELL_PARAMETERS")
    if v["k"] is not None:
        k1, k2, k3 = v["k"]
        t = sub1(t, r"K_POINTS gamma\n", f"K_POINTS automatic\n  {k1} {k2} {k3} 0 0 0\n", "K_POINTS")
    else:
        sub1(t, r"K_POINTS gamma\n", "K_POINTS gamma\n", "K_POINTS (presence)")
    if name != "identity":
        t = sub1(t, r"\A! ground state SCF of graphene 9x9 supercell\n", f"! R11 {name} ({v['desc']}) -- ground state SCF of graphene 9x9 supercell\n", "commentaire")
    return t, prefix, outdir


def pp_in(name, kind, template, prefix, outdir):
    s = "d" if kind == "defective" else "p"
    t = template
    t = sub1(t, r"prefix   = 'defect_\d+x\d+_[dp]'", f"prefix   = '{prefix}'", "pp prefix")
    t = sub1(t, r"outdir      = '[^']*'", f"outdir      = '{outdir}'", "pp outdir")
    if name != "identity":
        t = sub1(t, r"filplot  = 'Vks_\d+x\d+_[dp]'", f"filplot  = 'Vks_r11_{name}_{s}'", "pp filplot")
    return t


def templates(kind, tdir):
    d = os.path.join(tdir, SIZE, kind)
    return {f: read(os.path.join(d, f)) for f in ("scf.in", "pp.in")}


def build(name, v, tdir):
    out = {}
    for kind in ("defective", "pristine"):
        s = "d" if kind == "defective" else "p"
        T = templates(kind, tdir)
        scf, prefix, outdir = scf_in(name, v, kind, T["scf.in"])
        mem = "6G" if v["c"] > 30 else "4G"
        out[kind] = {"scf.in": scf, "pp.in": pp_in(name, kind, T["pp.in"], prefix, outdir),
                     "submit.scf": SUBMIT_SCF.format(account=ACCOUNT, name=name, s=s, mem=mem, time=v["time"], outdir=outdir),
                     "submit.pp": SUBMIT_PP.format(account=ACCOUNT, tag=f"r11_{name}", s=s, ppin="pp.in", ppout="pp.out")}
    return out


def check(tdir):
    ok = True
    for kind in ("defective", "pristine"):
        T = templates(kind, tdir)
        scf, prefix, outdir = scf_in("identity", IDENTITY, kind, T["scf.in"])
        same = scf == T["scf.in"]; ok &= same
        print(f"[check] {SIZE} {kind}: scf.in (variante identite) == gabarit : {same}")
        if not same:
            sys.stdout.writelines(difflib.unified_diff(T["scf.in"].splitlines(True), scf.splitlines(True), "gabarit", "regenere"))
        pp = pp_in("identity", kind, T["pp.in"], prefix, outdir)
        same = pp == T["pp.in"]; ok &= same
        print(f"[check] {SIZE} {kind}: pp.in (variante identite) == gabarit : {same}")
        n_pos = sum(1 for l in T["scf.in"].splitlines() if l.startswith("  C  "))
        print(f"[check] {SIZE} {kind}: {n_pos} lignes de positions dans le gabarit (attendu {161 if kind == 'defective' else 162})")
        ok &= n_pos == (161 if kind == "defective" else 162)
    for name, v in VARIANTS.items():                                            # chaque variante se construit (tous les motifs trouves)
        build(name, v, tdir)
    print(f"[check] {len(VARIANTS)} variantes construites sans erreur de motif")
    print("[check] PORTE", "PASS" if ok else "FAIL")
    return ok


def write_file(path, txt):
    if os.path.exists(path) and read(path) != txt:
        print(f"[write] {path} existe et differe : NON ecrase (supprimer a la main pour regenerer)"); return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(txt)
    print(f"[write] {path}")


def pp11(sizes, tdir):
    for S in sizes:
        for kind in ("defective", "pristine"):
            s = "d" if kind == "defective" else "p"
            T = read(os.path.join(tdir, S, kind, "pp.in"))
            prefix = re.search(r"prefix   = '([^']*)'", T).group(1); outdir = re.search(r"outdir      = '([^']*)'", T).group(1)
            d = os.path.join(WORK, "pp11", S, kind)
            write_file(os.path.join(d, "pp11.in"), PP11_IN.format(prefix=prefix, outdir=outdir, tag=S, s=s))
            write_file(os.path.join(d, "submit.pp11"), SUBMIT_PP.format(account=ACCOUNT, tag=f"pp11_{S}", s=s, ppin="pp11.in", ppout="pp11.out"))
            print(f"[pp11] {S} {kind}: prefix {prefix}, outdir {outdir} (le .save de production doit exister : charge-density + XML)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true"); ap.add_argument("--write", action="store_true"); ap.add_argument("--diff", action="store_true")
    ap.add_argument("--variants", default=None, help="liste de noms ou de phases (1, 2, 2b) ; defaut : toutes")
    ap.add_argument("--pp11", default=None, help="tailles (ex. 9x9,6x6,12x12) : ecrit pp11/<taille>/...")
    ap.add_argument("--templates", default=os.path.join(GQ, "defects", "super_cell"), help="repertoire des gabarits <taille>/{pristine,defective}/")
    a = ap.parse_args()
    if not os.path.isdir(os.path.join(a.templates, SIZE)):
        raise SystemExit(f"gabarits introuvables : {a.templates}/{SIZE} (option --templates)")
    if a.check and not check(a.templates):
        sys.exit(1)
    if a.pp11:
        pp11(a.pp11.split(","), a.templates)
    names = list(VARIANTS)
    if a.variants:
        names = []
        for x in a.variants.split(","):
            names += PHASES.get(x, [x])
    for name in names:
        v = VARIANTS[name]; out = build(name, v, a.templates)
        print(f"[{name}] {v['desc']} ; K_POINTS {'gamma' if v['k'] is None else v['k']} ; {v['smear']} {v['degauss']:g} Ry ; "
              f"2D {v['iso2d']} ; c {v['c']:g} bohr ; scf {v['time']}")
        for kind in ("defective", "pristine"):
            d = os.path.join(WORK, "variants", name, kind)
            if a.diff:
                T = templates(kind, a.templates)
                for f in ("scf.in", "pp.in"):
                    ref = [l for l in T[f].splitlines(True) if not l.startswith("  C  ")]; gen = [l for l in out[kind][f].splitlines(True) if not l.startswith("  C  ")]
                    dd = list(difflib.unified_diff(ref, gen, f"{SIZE}/{kind}/{f}", f"{name}/{kind}/{f}", n=0))
                    print("".join(dd) if dd else f"--- {name}/{kind}/{f}: identique au gabarit\n", end="")
            if a.write:
                for f, txt in out[kind].items():
                    write_file(os.path.join(d, f), txt)


if __name__ == "__main__":
    main()
