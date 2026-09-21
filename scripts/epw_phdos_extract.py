#!/usr/bin/env python
"""Extrait la DOS de phonons matdyn (dos=.true.) d'un répertoire phonons/ vers results/epw/phdos_<tag>.npz, avec les
paramètres réels lus dans les fichiers d'entrée (grille q de la DOS, deltaE et degauss de matdyn en cm⁻¹, degauss MV
électronique de scf.in en Ry, grilles k/q de scf.in et ph.in). Chiffres bruts, aucune renormalisation.
Usage : python scripts/epw_phdos_extract.py --dir <phonons/> --tag 24k24q_mv0.02 [--dos-in matdyn.dos.in] [--dos graphene.dos]"""
import argparse, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--dos-in", default="matdyn.dos.in"); ap.add_argument("--dos", default="graphene.dos")
ap.add_argument("--scf-in", default="scf.in"); ap.add_argument("--ph-in", default="ph.in"); ap.add_argument("--out", default=None)
a = ap.parse_args()
def nml(path):
    """Valeurs de namelist Fortran (clé -> chaîne), commentaires « ! » retirés ; plusieurs affectations par ligne acceptées."""
    d = {}
    for line in open(path, errors="replace"):
        line = line.split("!")[0]
        for k, v in re.findall(r"([A-Za-z_]\w*(?:\(\d+\))?)\s*=\s*([^,/\n]+)", line): d[k.strip().lower()] = v.strip().strip("'\"")
    return d
D = np.loadtxt(os.path.join(a.dir, a.dos))
m = nml(os.path.join(a.dir, a.dos_in)); s = nml(os.path.join(a.dir, a.scf_in)); p = nml(os.path.join(a.dir, a.ph_in))
kline = None
for i, line in enumerate(L := open(os.path.join(a.dir, a.scf_in)).read().splitlines()):
    if line.strip().upper().startswith("K_POINTS"): kline = L[i + 1].split()[:3]
out = dict(freq_cm=D[:, 0], dos=D[:, 1], pdos=D[:, 2:] if D.shape[1] > 2 else np.zeros((len(D), 0)),
           nq_dos=np.array([int(m.get("nk1", 0)), int(m.get("nk2", 0)), int(m.get("nk3", 0))]),
           deltaE_cm=float(m.get("deltae", np.nan)), degauss_dos_cm=float(m.get("degauss", np.nan)), asr=str(m.get("asr", "")),
           flfrc=str(m.get("flfrc", "")), sigma_mv_Ry=float(s.get("degauss", np.nan)), smearing=str(s.get("smearing", "")),
           nk_scf=np.array([int(x) for x in kline]) if kline else np.zeros(3, int),
           nq_ph=np.array([int(p.get("nq1", 0)), int(p.get("nq2", 0)), int(p.get("nq3", 0))]),
           units="freq cm^-1 ; DOS telle qu'ecrite par matdyn (etats/cm^-1/cellule) ; deltaE, degauss_dos en cm^-1 ; sigma_mv en Ry",
           source=os.path.abspath(os.path.join(a.dir, a.dos)))
outp = a.out or f"results/epw/phdos_{a.tag}.npz"; np.savez(outp, **out)
dE = float(np.median(np.diff(D[:, 0])))
print(f"écrit {outp} : {len(D)} points, {D[0,0]:.1f}–{D[-1,0]:.1f} cm^-1, pas {dE:.3f} cm^-1 ; intégrale DOS = {D[:,1].sum()*dE:.4f} (3·nat = 6 attendu si non tronquée) ; "
      f"grille DOS {out['nq_dos'].tolist()}, deltaE {out['deltaE_cm']} cm^-1, degauss DOS {out['degauss_dos_cm']} cm^-1 ; scf : {out['smearing']} {out['sigma_mv_Ry']} Ry, k {out['nk_scf'].tolist()} ; ph : q {out['nq_ph'].tolist()}")
