# R1 — Relaxation de la supercellule 9x9 avec lacune (préparation article, hors mémoire)

Date de lancement : 2026-09-22. Source confirmée : `graphene/qe/defects/super_cell/9x9/defective/scf.in`
(161 atomes, lacune sous-réseau A à (13/27, 13/27), voisins sous-coordonnés = atomes 64, 80, 81).

## Jobs SLURM

| calcul | job ID | soumis | ressources |
|---|---|---|---|
| nspin1 | 21587934 | 2026-09-22 | rrg-cotemich-ac, 64 tâches, 4 G/cœur, 3 j |
| nspin2 | 21587935 | 2026-09-22 | idem |

## Géométrie de départ (reproductibilité)

- Générateur : `make_inputs.py`, md5 = 53fa3f87439ad3eb42e3a1b045c1c1f0
- Perturbation initiale (copie de `perturbation.log`) :

```
d(64-80) idéal = 4.659785 bohr = 2.465852 A
u (64->80) = [ 0. -1.  0.]
atome 64: cart. [33.62909922  2.32989234  0.        ] -> [33.62909922  2.27320055  0.        ] bohr ; delta cart = [-2.01921813e-15 -5.66917837e-02  0.00000000e+00] bohr = [-1.06852422e-15 -3.00000000e-02  0.00000000e+00] A (|d|=0.03000 A) ; delta cryst = [ 0.0013518 -0.0013518  0.       ]
atome 80: cart. [33.62909922 -2.32989234  0.        ] -> [33.62909922 -2.27320055  0.        ] bohr ; delta cart = [-2.01921813e-15  5.66917837e-02  0.00000000e+00] bohr = [-1.06852422e-15  3.00000000e-02  0.00000000e+00] A (|d|=0.03000 A) ; delta cryst = [-0.0013518  0.0013518  0.       ]
d(64-80) perturbé = 2.405852 A
```

## Résultats (2026-09-22, les deux jobs COMPLETED, ExitCode 0:0)

Sorties : `nspin1/relax.out`, `nspin2/relax.out` ; analyse brute complète : `analyse_2026-09-22.txt` (`analyze_relax.py`).
Convention : distances par image minimale ; déplacements mesurés par rapport aux positions idéales du réseau (scf.in du ch. 4) ; spglib 2.5.0 via ASE 3.26.0, symprec 1e-3 Å.

| grandeur | nspin1 | nspin2 |
|---|---|---|
| job ID / durée | 21587934 / 1 h 03 min 37 s | 21587935 / 3 h 18 min 36 s |
| pas BFGS / cycles SCF | 27 / 28 | 23 / 24 |
| itérations SCF max par pas | 26 | 54 (1er pas), ≤ 42 ensuite |
| E finale (Ry) | −1939.06466920 | −1939.09828972 |
| force totale finale (Ry/bohr) | 3.39e-4 | 4.94e-4 |
| composante de force max finale (Ry/bohr) | 7.33e-5 | 9.04e-5 |
| d(64-80) idéal → perturbé → final (Å) | 2.46585 → 2.40585 → 2.19194 | 2.46585 → 2.40585 → 1.98798 |
| d(64-81) final / d(80-81) final (Å) | 2.57319 / 2.57319 | 2.57260 / 2.57259 |
| déplacement max vs idéal (Å) | 0.14608 (atomes 64 et 80) | 0.24811 (atomes 64 et 80) |
| déplacement de 81 vs idéal, dans le plan (Å) | 0.14181 | 0.17049 |
| déplacement hors plan max, tous atomes (Å) | 5.2e-7 | 6.0e-8 |
| dz(81) (Å) | +5.2e-7 | +6.0e-8 |
| magnétisation totale / absolue finale (µB) | — | +1.35 / 2.45 (précision d'impression pw.x : 0.01) |
| spglib, positions initiales perturbées | Amm2 (n° 38, 4 op.) ; site lacune mm2 | idem |
| spglib, positions finales | Amm2 (n° 38, 4 op.) ; site lacune mm2 | Amm2 (n° 38, 4 op.) ; site lacune mm2 |

E(nspin2) − E(nspin1) = −0.03362052 Ry = **−457.43 meV** (nspin2 plus bas).
Deux des trois distances entre atomes sous-coordonnés se sont raccourcies l'une par rapport aux autres : d(64-80) < d(64-81) = d(80-81) dans les deux cas (écart 0.381 Å en nspin1, 0.585 Å en nspin2).

### Longueur 64–80 à chaque pas BFGS (avec E, force totale ; nspin2 : magnétisation totale et absolue en µB)

| pas | nspin1 : E (Ry) | Ftot | d(64-80) (Å) | nspin2 : E (Ry) | Ftot | d(64-80) (Å) | m_tot | m_abs |
|---|---|---|---|---|---|---|---|---|
| 0 | -1939.03085605 | 0.270970 | 2.40585 | -1939.04549516 | 0.264303 | 2.40585 | +1.4400 | 3.0000 |
| 1 | -1939.04909863 | 0.188660 | 2.47737 | -1939.06160392 | 0.212129 | 2.43402 | +1.2400 | 2.4200 |
| 2 | -1939.05841262 | 0.050791 | 2.48412 | -1939.07370893 | 0.065842 | 2.43027 | +1.3200 | 2.6100 |
| 3 | -1939.05921236 | 0.031723 | 2.48511 | -1939.07691819 | 0.056595 | 2.40828 | +1.4000 | 2.7600 |
| 4 | -1939.05960934 | 0.012811 | 2.49036 | -1939.08183457 | 0.062523 | 2.36300 | +1.4800 | 2.8700 |
| 5 | -1939.05976129 | 0.011168 | 2.49419 | -1939.08632218 | 0.072861 | 2.28529 | +1.3800 | 2.6600 |
| 6 | -1939.05988408 | 0.010902 | 2.49365 | -1939.08960524 | 0.052687 | 2.22845 | +1.3600 | 2.5900 |
| 7 | -1939.06009044 | 0.013398 | 2.48568 | -1939.09244996 | 0.045644 | 2.15422 | +1.3600 | 2.5600 |
| 8 | -1939.06034636 | 0.017692 | 2.47209 | -1939.09443580 | 0.037652 | 2.08781 | +1.3600 | 2.5200 |
| 9 | -1939.06069827 | 0.023219 | 2.45026 | -1939.09581547 | 0.029905 | 2.03062 | +1.3500 | 2.4900 |
| 10 | -1939.06111195 | 0.028578 | 2.41563 | -1939.09659623 | 0.024358 | 2.00332 | +1.3500 | 2.4700 |
| 11 | -1939.06170574 | 0.028375 | 2.38308 | -1939.09717708 | 0.019737 | 1.99091 | +1.3500 | 2.4600 |
| 12 | -1939.06243111 | 0.025490 | 2.34727 | -1939.09761733 | 0.017245 | 1.98459 | +1.3500 | 2.4500 |
| 13 | -1939.06337761 | 0.020328 | 2.29356 | -1939.09789877 | 0.013960 | 1.98286 | +1.3500 | 2.4500 |
| 14 | -1939.06411322 | 0.017630 | 2.21539 | -1939.09810774 | 0.010625 | 1.98157 | +1.3500 | 2.4500 |
| 15 | -1939.06430778 | 0.013264 | 2.18787 | -1939.09819740 | 0.007901 | 1.98067 | +1.3500 | 2.4500 |
| 16 | -1939.06443237 | 0.008379 | 2.17980 | -1939.09823356 | 0.004795 | 1.98108 | +1.3500 | 2.4500 |
| 17 | -1939.06451733 | 0.007285 | 2.18129 | -1939.09825374 | 0.003813 | 1.98211 | +1.3500 | 2.4500 |
| 18 | -1939.06455547 | 0.006117 | 2.18661 | -1939.09826885 | 0.003216 | 1.98378 | +1.3500 | 2.4500 |
| 19 | -1939.06460117 | 0.005370 | 2.19193 | -1939.09827883 | 0.002602 | 1.98556 | +1.3500 | 2.4500 |
| 20 | -1939.06464150 | 0.004659 | 2.19364 | -1939.09828564 | 0.001928 | 1.98736 | +1.3500 | 2.4500 |
| 21 | -1939.06465618 | 0.003473 | 2.19143 | -1939.09828850 | 0.001229 | 1.98806 | +1.3500 | 2.4500 |
| 22 | -1939.06466347 | 0.001929 | 2.19007 | -1939.09828941 | 0.000798 | 1.98810 | +1.3500 | 2.4500 |
| 23 | -1939.06466631 | 0.001271 | 2.19081 | -1939.09828972 | 0.000494 | 1.98798 | +1.3500 | 2.4500 |
| 24 | -1939.06466772 | 0.000930 | 2.19208 | — | — | — | — | — |
| 25 | -1939.06466857 | 0.000629 | 2.19292 | — | — | — | — | — |
| 26 | -1939.06466898 | 0.000497 | 2.19264 | — | — | — | — | — |
| 27 | -1939.06466920 | 0.000339 | 2.19194 | — | — | — | — | — |

## .save finaux

Scratch : `qe_tmp/vacancy_relaxed/nspin{1,2}/vac_9x9_relax_nspin{1,2}.save` : nspin1 2.1 G (wfc1.hdf5), nspin2 4.2 G (wfcup1.hdf5 + wfcdw1.hdf5), plus charge-density.hdf5, data-file-schema.xml, C.upf.
Miroir /project (règle du 2026-09-17) : `graphene/qe/qe_tmp_backup/vacancy_relaxed/nspin{1,2}/` avec `MD5SUMS_nspin{1,2}_2026-09-22.txt` vérifiés.

## Fichiers du répertoire

`make_inputs.py`, `perturbation.log`, `analyze_relax.py`, `analyse_2026-09-22.txt`, `nspin{1,2}/{relax.in,relax.out,submit.relax,JOBID,slurm-*.out,slurm-*.err}`, `R1_rapport.md`.
Aucun post-traitement effectué (pas de ΔV, pw2wannier90, Wannier90 ni pipeline e-d). Rien n'est commité dans le dépôt.
