# R1b — Contrôles 3×3×1 sur les géométries relaxées du 9×9 (suite de R1, hors mémoire)

Date : 2026-09-22. Répertoire `k3x3/` : `make_k3x3.py` (générateur, positions finales extraites du bloc « Begin/End final coordinates » de
`nspin{1,2}/relax.out`, assertions nat = 161 et d(64–80) = 2.19194 / 1.98798 Å passées ; positions idéales identiques au `scf.in` du ch. 4),
`analyze_k3x3.py`, `analyse_k3x3_2026-09-22.txt`, un sous-répertoire par calcul (`scf.in`, `scf.out`, `submit.scf`, `JOBID`, `slurm-*`).
Paramètres : ceux de R1 (ecutwfc 100 Ry, mv 0.01 Ry, conv_thr 1e-10, nosym/noinv, nbnd 352, assume_isolated 2D, C1 sur 64/80/81 avec
starting_magnetization 0.5 pour nspin2) sauf calculation = 'scf', tprnfor = .true., K_POINTS automatic 3 3 1 0 0 0 (9 points k, aucune réduction).
Aucun miroir des `.save` de `k3x3/` (scratch `qe_tmp/vacancy_relaxed/k3x3/<nom>/`), rien supprimé, rien commité.

## Jobs SLURM (rrg-cotemich-ac, 64 tâches, 4 G/cœur, 12 h ; projwfc 2 h)

```
ideal_nspin1 21600427
ideal_nspin2 21600428
relax1_nspin1 21600429
relax2_nspin2 21600430
projwfc gamma_nspin2 21600431
projwfc k3x3_relax2_nspin2 21600432 (afterok:21600430)
```

## 1. Énergie de relaxation à Γ (lecture seule)

Diff `super_cell/9x9/defective/scf.in` vs `nspin1/relax.in` hors positions : seules différences = calculation, prefix/outdir, forc_conv_thr,
etot_conv_thr, nstep, tprnfor, &IONS (sans effet sur l'énergie scf) et nosym/noinv (le scf du ch. 4 utilisait 4 opérations de symétrie).
Comparaison faite avec cette réserve. scf ch. 4 idéal : F = −1939.03777059 Ry, −TS = −0.01410526 Ry, E = −1939.02366533 Ry (24 itérations).

## 2–3. Tableau par calcul

| calcul | JOB DONE | nk | F = ! total energy (Ry) | E = F+TS (Ry) | −TS (Ry) | force max comp. (Ry/bohr) | |F_atome| max (atome) | force totale | m_tot / m_abs (µB) | E_F (eV) | itér. SCF | WALL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 nspin1 (Γ, relax) | True | 1 | -1939.06466920 | -1939.05224269 | -0.01242651 | 7.33e-05 | 7.38e-05 (80) | 0.00033900 | — | -4.2239 | 12 | 3788.51s |
| R1 nspin2 (Γ, relax) | True | 1 | -1939.09828972 | -1939.08982060 | -0.00846913 | 9.04e-05 | 1.12e-04 (81) | 0.00049400 | +1.35 / 2.45 | -4.2244 | 10 | 11886.87s |
| ideal_nspin1 (3x3) | True | 9 | -1939.10377543 | -1939.10210397 | -0.00167147 | 8.04e-02 | 8.25e-02 (79) | 0.24558400 | — | -4.4035 | 23 | 3795.79s |
| ideal_nspin2 (3x3) | non (scancel) | 9 | −1939.10781398 (dernière itération, sans « ! ») | — (−TS non imprimé avant convergence) | — | — | — | — | +1.10 / 1.60 | — | 46 (annulé à la 47e) | 4 h 40 min (CANCELLED) |
| relax1_nspin1 (3x3) | True | 9 | -1939.12574149 | -1939.12258180 | -0.00315969 | 1.51e-02 | 1.73e-02 (82) | 0.04907200 | — | -4.4251 | 21 | 2845.52s |
| relax2_nspin2 (3x3) | True | 9 | -1939.15438279 | -1939.15366611 | -0.00071668 | 6.46e-03 | 6.46e-03 (61) | 0.04326500 | +1.54 / 2.13 | -4.4585 | 24 | 13617.60s |

`ideal_nspin2 (3x3)` (job 21600428) : **annulé par scancel (décision de Greg, 2026-09-22) avant convergence** — précision SCF 6.0e-10 Ry à la
dernière itération complète (46e ; la précision a stagné autour de 6e-9 Ry entre les itérations 27 et 41, puis est redescendue : 2.0e-9, 1.2e-9, 6.0e-10),
seuil 1e-10 non atteint ; énergie « total energy » de la 46e itération = −1939.10781398 Ry (valeur non finale, sans « ! »), magnétisation +1.10 / 1.60 µB.
−TS et E = F + TS ne sont pas imprimés avant convergence : les différences le concernant ne sont données qu'avec F. Pas de relance.

Forces : « force max comp. » = plus grande composante cartésienne ; « |F_atome| max » = plus grande norme par atome (numéro d'atome) ; « force totale » = « Total force » de pw.x.
Sur les positions relaxées à Γ, les forces à 3×3 ne sont pas nulles (voir tableau : relax1_nspin1 0.049, relax2_nspin2 0.043 Ry/bohr) : la géométrie bougerait sous 3×3.
Les positions idéales à 3×3 donnent une force totale de 0.246 Ry/bohr (ideal_nspin1).

## 4. Différences (Ry et meV), avec F = « ! total energy » et avec E = F + TS

| différence | avec F | avec E |
|---|---|---|
| ΔE_spin(Γ) = R1 nspin2 − R1 nspin1 | -0.03362052 Ry = -457.43 meV | -0.03757791 Ry = -511.27 meV |
| ΔE_spin(3×3) = relax2_nspin2 − relax1_nspin1 | -0.02864130 Ry = -389.68 meV | -0.03108431 Ry = -422.92 meV |
| ΔE_spin,idéal(3×3) = ideal_nspin2 − ideal_nspin1 | -0.00403855 Ry = -54.95 meV (ideal_nspin2 non convergé, 6e-10 Ry) | — (E indisponible pour ideal_nspin2) |
| ΔE_relax,1(3×3) = relax1_nspin1 − ideal_nspin1 | -0.02196606 Ry = -298.86 meV | -0.02047783 Ry = -278.62 meV |
| ΔE_relax,2(3×3) = relax2_nspin2 − ideal_nspin2 | -0.04656881 Ry = -633.60 meV (ideal_nspin2 non convergé, 6e-10 Ry) | — (E indisponible pour ideal_nspin2) |
| ΔE_relax(Γ) = R1 nspin1 − scf ch. 4 idéal | -0.02689861 Ry = -365.97 meV | -0.02857736 Ry = -388.81 meV |

Avec le smearing Marzari-Vanderbilt, F (« ! total energy ») est la colonne de référence.

## 3. Moments de Löwdin par atome (projwfc.x, lwrite_overlaps = .false., filproj ; Γ = `.save` de R1 nspin2, 3×3 = relax2_nspin2)

| moment de Löwdin (µB = up − down) | Γ (R1 nspin2) | 3×3 (relax2_nspin2) |
|---|---|---|
| atome 81 | +0.9111 | +0.9408 |
| atome 64 | +0.0434 | +0.0929 |
| atome 80 | +0.0434 | +0.0929 |
| somme sous-réseau A (80 atomes, celui de la lacune) | -0.7294 | -0.4158 |
| somme sous-réseau B (81 atomes) | +2.0610 | +1.9361 |
| total (161 atomes) | +1.3316 | +1.5203 |
| charge de Löwdin totale, Γ (R1 nspin2) : 636.983 e (161 atomes, 161 entrées) ; spilling 0.0109 |
| charge de Löwdin totale, 3×3 (relax2_nspin2) : 636.986 e (161 atomes, 161 entrées) ; spilling 0.0109 |


Composantes de l'atome 81 (polarisation, projwfc) : Γ : s 0.1474, p 0.7637 (pz 0.0924, px 0.6555, py 0.0157) ; 3×3 : s 0.1487, p 0.7921 (pz 0.1206, px 0.6553, py 0.0162).
Magnétisation totale pw.x vs somme de Löwdin : Γ 1.35 vs 1.332 µB ; 3×3 1.54 vs 1.520 µB (spilling 0.0109 dans les deux cas).
