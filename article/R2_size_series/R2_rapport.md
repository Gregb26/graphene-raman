# R2 — Série en taille 5×5 → 12×12 de la lacune relaxée (préparation article, hors mémoire)

Statut : PRODUCTION. Protocole : celui de R1 (`../9x9/`, `article/R1_vacancy_relaxed/`), inchangé : pw.x v7.5, calculation = 'relax',
cellule fixe, bfgs, nstep 200, forc_conv_thr 1e-4, etot_conv_thr 1e-5, tprnfor, nosym/noinv, ecutwfc/mv/conv_thr/assume_isolated/nbnd du
`scf.in` du ch. 4 de la même taille (`super_cell/NxN/defective/scf.in`, phase 0 : seuls nat/ntyp, nbnd, cellule et positions diffèrent
du 9×9) ; nspin1 et nspin2 (C1 sur les 3 voisins sous-coordonnés, starting_magnetization(2) = 0.5) ; perturbation initiale 0.03 Å par atome
sur la paire des deux voisins d'indice le plus bas, l'un vers l'autre dans le plan (`NxN/perturbation.log`).
Générateur : `make_inputs_series.py` (md5 dc001ce0d92dbd75f6eb376d0ad54343), dérivé de `make_inputs.py` de R1 ; régression : le 9×9
régénéré reproduit les inputs R1 hors titre/prefix/outdir. Table de la phase 0 : `R2_phase0_table.txt`.
SLURM : rrg-cotemich-ac, 64 tâches, 4 G/cœur, 3 j (10, 11, 12) ou 1 j (5–8) ; jobs et état final dans `R2_jobs.md`.

## Emplacements (règle 5 de CLAUDE.md)

- Répertoire de travail : `graphene/qe/defects/super_cell_relaxed/series/NxN/{nspin1,nspin2}/` (`relax.in`, `relax.out`, `submit.relax`, `JOBID`, `slurm-*`), `NxN/nspin2/projwfc/` (projwfc.x).
- outdir scratch : `qe_tmp/vacancy_relaxed/series/NxN/nspinN/vac_NxN_relax_nspinN.save` (+ `.xml`).
- Miroir /project : `graphene/qe/qe_tmp_backup/vacancy_relaxed/series/NxN/nspinN/` (§ « .save finaux »).
- Copie versionnée : `ab-initio-defects/article/R2_size_series/` (inputs, submit, scripts, rapports, `relax.out` < 5 Mo ; pas de `.save`, `slurm-*`, `pdos_*`, `proj_*`, `projwfc.out`).

## Analyse (2026-09-24)

Script : `analyze_relax_series.py` (md5 ec518e72535ff87a5df59944439fdbdb), généralisation de `9x9/analyze_relax.py` (R1) ; sortie brute
complète par (N, spin), avec suivi par pas BFGS (E, ΔE, force totale, composante max, itérations SCF, d(paire), m_tot/m_abs) :
`analyse_R2_2026-09-24.txt`. Le 9×9 est relu depuis `../9x9/` avec le même script : ses chiffres coïncident avec `R1_rapport.md` et `R1b_rapport.md`.
Conventions : F = « ! total energy » (Ry) comme référence ; E scf idéal = F du `scf.out` du ch. 4 (nspin1, 4 opérations de symétrie,
réserve déjà notée en R1b §1) ; distances par image minimale ; déplacements vs positions idéales du `scf.in` du ch. 4 ;
spglib 2.5.0 via ASE 3.26.0, symprec 1e-3 Å, « site » = symétrie du site de la lacune (atome factice à la position idéale) ;
sous-réseaux par la règle s_uc = mod(N·s, 1) (A ≈ 1/3, B ≈ 2/3), même règle que `tag_vacancy_sublattice.py` ;
« K replié sur Γ » = oui si N multiple de 3. « 3e » = troisième voisin (celui hors paire), « paire » = (p1, p2).
Magnétisations pw.x imprimées à 0.01 µB près. projwfc.x v7.5 : `NxN/nspin2/projwfc/projwfc.in` (identique à R1b hors prefix/outdir/noms
de fichiers : lwrite_overlaps = .false., Emin −30, Emax 10, DeltaE 1.0, degauss 0.01, ngauss −1), jobs 21735568–21735574 (5 s à 2 min).

## Tableau récapitulatif (une ligne par (N, spin) ; F = « ! total energy » ; 9x9 = R1 relu avec le même script)

| N | nat | K replié sur Γ | spin | job ID | état | durée (sacct) | pas BFGS / cycles SCF | SCF non conv. | F finale (Ry) | F scf idéal ch. 4 (Ry) | ΔE_relax (Ry) | ΔE_relax (meV) | force totale finale (Ry/bohr) | comp. max (Ry/bohr) | d(paire) idéal → perturbé → final (Å) | d(p1–3e) / d(p2–3e) final (Å) | dépl. max vs idéal (Å) (atome) | dépl. 3e dans le plan (Å) | hors plan max (Å) | m_tot / m_abs (µB) | spglib initial | spglib final |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 49 | non | nspin1 | 21666594 | COMPLETED | 00:02:47 | 16 / 17 | 0 | -589.72287142 | -589.69692552 | -0.02594590 | -353.01 | 4.09e-04 | 9.43e-05 | 2.46585 → 2.40585 → 2.53032 | 2.53079 / 2.53079 | 0.03761 (25) | 0.03761 | 6.0e-07 | — | Amm2 (n° 38, 4 op.), site mm2 | P-6m2 (n° 187, 12 op.), site -6m2 |
| 5 | 49 | non | nspin2 | 21666595 | COMPLETED | 00:07:14 | 19 / 20 | 0 | -589.73483173 | -589.69692552 | -0.03790621 | -515.74 | 2.25e-04 | 4.78e-05 | 2.46585 → 2.40585 → 2.20393 | 2.59215 / 2.59212 | 0.14660 (25) | 0.14660 | 6.0e-08 | +2.00 / 2.89 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 6 | 71 | oui | nspin1 | 21666592 | COMPLETED | 00:06:05 | 21 / 22 | 0 | -854.72638978 | -854.69845523 | -0.02793455 | -380.07 | 4.69e-04 | 9.66e-05 | 2.46585 → 2.40585 → 2.24641 | 2.58226 / 2.58226 | 0.13439 (40) | 0.13439 | 1.9e-07 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 6 | 71 | oui | nspin2 | 21666593 | COMPLETED | 00:13:50 | 17 / 18 | 0 | -854.76177739 | -854.69845523 | -0.06332216 | -861.54 | 3.71e-04 | 8.59e-05 | 2.46585 → 2.40585 → 2.06819 | 2.57963 / 2.57962 | 0.21016 (30) | 0.15972 | 7.5e-08 | +1.97 / 5.27 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 7 | 97 | non | nspin1 | 21666590 | COMPLETED | 00:09:17 | 17 / 18 | 0 | -1168.09808835 | -1168.07147324 | -0.02661511 | -362.12 | 4.82e-04 | 8.46e-05 | 2.46585 → 2.40585 → 2.50243 | 2.50270 / 2.50270 | 0.04301 (60) | 0.02136 | 1.2e-06 | — | Amm2 (n° 38, 4 op.), site mm2 | P-6m2 (n° 187, 12 op.), site -6m2 |
| 7 | 97 | non | nspin2 | 21666591 | COMPLETED | 00:27:37 | 22 / 23 | 0 | -1168.11414258 | -1168.07147324 | -0.04266934 | -580.55 | 3.53e-04 | 6.24e-05 | 2.46585 → 2.40585 → 2.10980 | 2.58904 / 2.58903 | 0.19041 (48) | 0.16133 | 1.0e-07 | +2.00 / 2.84 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 8 | 127 | non | nspin1 | 21666588 | COMPLETED | 00:15:01 | 21 / 22 | 0 | -1529.49014442 | -1529.46517508 | -0.02496934 | -339.73 | 5.47e-04 | 8.08e-05 | 2.46585 → 2.40585 → 2.45234 | 2.53385 / 2.53385 | 0.06048 (73) | 0.06048 | 1.2e-05 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 8 | 127 | non | nspin2 | 21666589 | COMPLETED | 00:46:42 | 23 / 24 | 0 | -1529.51310712 | -1529.46517508 | -0.04793204 | -652.15 | 3.49e-04 | 8.00e-05 | 2.46585 → 2.40585 → 2.07452 | 2.58478 / 2.58477 | 0.20690 (72) | 0.16479 | 1.1e-07 | +1.94 / 2.85 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 9 | 161 | oui | nspin1 | 21587934 | COMPLETED | 01:03:37 | 27 / 28 | 0 | -1939.06466920 | -1939.03777059 | -0.02689861 | -365.97 | 3.39e-04 | 7.33e-05 | 2.46585 → 2.40585 → 2.19194 | 2.57319 / 2.57319 | 0.14608 (64) | 0.14181 | 5.2e-07 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 9 | 161 | oui | nspin2 | 21587935 | COMPLETED | 03:18:36 | 23 / 24 | 0 | -1939.09828972 | -1939.03777059 | -0.06051913 | -823.40 | 4.94e-04 | 9.04e-05 | 2.46585 → 2.40585 → 1.98798 | 2.57260 / 2.57259 | 0.24811 (80) | 0.17049 | 6.0e-08 | +1.35 / 2.45 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 10 | 199 | non | nspin1 | 21666586 | COMPLETED | 00:45:26 | 15 / 16 | 0 | -2396.97400552 | -2396.94841705 | -0.02558847 | -348.15 | 5.15e-04 | 7.65e-05 | 2.46585 → 2.40585 → 2.48474 | 2.51089 / 2.51089 | 0.04843 (91) | 0.03301 | 9.0e-08 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 10 | 199 | non | nspin2 | 21666587 | COMPLETED | 02:41:26 | 26 / 27 | 0 | -2396.99673446 | -2396.94841705 | -0.04831741 | -657.39 | 3.62e-04 | 6.08e-05 | 2.46585 → 2.40585 → 2.04976 | 2.58525 / 2.58525 | 0.21911 (90) | 0.16917 | 1.3e-07 | +1.94 / 2.80 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 11 | 241 | non | nspin1 | 21666584 | COMPLETED | 02:43:46 | 26 / 27 | 0 | -2902.96097163 | -2902.93599922 | -0.02497241 | -339.77 | 5.91e-04 | 9.12e-05 | 2.46585 → 2.40585 → 2.40959 | 2.54232 / 2.54232 | 0.07790 (121) | 0.07790 | 7.1e-05 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 11 | 241 | non | nspin2 | 21666585 | COMPLETED | 07:34:46 | 27 / 28 | 0 | -2902.98659119 | -2902.93599922 | -0.05059197 | -688.34 | 5.05e-04 | 8.63e-05 | 2.46585 → 2.40585 → 2.02325 | 2.58144 / 2.58145 | 0.23162 (100) | 0.17113 | 1.2e-07 | +1.77 / 2.57 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 12 | 287 | oui | nspin1 | 21666582 | COMPLETED | 03:33:17 | 31 / 32 | 0 | -3457.11840254 | -3457.09139568 | -0.02700686 | -367.45 | 5.38e-04 | 9.79e-05 | 2.46585 → 2.40585 → 2.15980 | 2.56934 / 2.56934 | 0.16084 (131) | 0.14636 | 3.1e-07 | — | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |
| 12 | 287 | oui | nspin2 | 21666583 | COMPLETED | 07:54:45 | 29 / 30 | 0 | -3457.15209844 | -3457.09139568 | -0.06070276 | -825.90 | 4.78e-04 | 9.77e-05 | 2.46585 → 2.40585 → 1.94876 | 2.56827 / 2.56828 | 0.26691 (132) | 0.17447 | 2.4e-07 | +1.13 / 1.73 | Amm2 (n° 38, 4 op.), site mm2 | Amm2 (n° 38, 4 op.), site mm2 |

## ΔE_spin = F(nspin2) − F(nspin1) et moments de Löwdin (projwfc.x sur le .save nspin2, lwrite_overlaps = .false.)

| N | K replié sur Γ | lacune / paire / 3e | ΔE_spin (Ry) | ΔE_spin (meV) | m_tot / m_abs pw.x (µB) | Löwdin atome 3e (µB) | Löwdin paire (µB) | Löwdin somme sous-réseau de la lacune (µB, n) | Löwdin somme autre sous-réseau (µB, n) | Löwdin total (µB) | comp. 3e : s / p (pz, px, py) | charge Löwdin (e) | spilling |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | non | A / (16, 24) / 25 | -0.01196031 | -162.73 | +2.00 / 2.89 | +0.9954 | +0.1929 / +0.1928 | -0.6514 (A, 24) | +2.6271 (B, 25) | +1.9757 | 0.1524 / 0.8429 (0.1801, 0.6455, 0.0173) | 193.854 | 0.011 |
| 6 | oui | B / (29, 30) / 40 | -0.03538761 | -481.47 | +1.97 / 5.27 | +0.9738 | +0.1217 / +0.1217 | -2.1947 (B, 35) | +4.1454 (A, 36) | +1.9507 | 0.1513 / 0.8225 (0.1569, 0.1754, 0.4903) | 280.893 | 0.0109 |
| 7 | non | A / (36, 48) / 49 | -0.01605423 | -218.43 | +2.00 / 2.84 | +1.0220 | +0.1978 / +0.1978 | -0.6075 (A, 48) | +2.5805 (B, 49) | +1.9730 | 0.1520 / 0.8700 (0.2016, 0.6509, 0.0176) | 383.764 | 0.0109 |
| 8 | non | A / (58, 72) / 73 | -0.02296270 | -312.42 | +1.94 / 2.85 | +0.9842 | +0.1566 / +0.1566 | -0.6459 (A, 63) | +2.5662 (B, 64) | +1.9203 | 0.1498 / 0.8343 (0.1677, 0.6496, 0.0169) | 502.465 | 0.0109 |
| 9 | oui | A / (64, 80) / 81 | -0.03362052 | -457.43 | +1.35 / 2.45 | +0.9111 | +0.0434 / +0.0434 | -0.7294 (A, 80) | +2.0610 (B, 81) | +1.3316 | 0.1474 / 0.7637 (0.0924, 0.6555, 0.0157) | 636.983 | 0.0109 |
| 10 | non | B / (89, 90) / 108 | -0.02272894 | -309.24 | +1.94 / 2.80 | +0.9957 | +0.1607 / +0.1607 | -0.6142 (B, 99) | +2.5276 (A, 100) | +1.9134 | 0.1503 / 0.8454 (0.1757, 0.1759, 0.4937) | 787.333 | 0.0109 |
| 11 | non | A / (100, 120) / 121 | -0.02561956 | -348.57 | +1.77 / 2.57 | +0.9654 | +0.1215 / +0.1215 | -0.5583 (A, 120) | +2.3086 (B, 121) | +1.7503 | 0.1491 / 0.8162 (0.1455, 0.6543, 0.0165) | 953.511 | 0.0109 |
| 12 | oui | B / (131, 132) / 154 | -0.03369590 | -458.46 | +1.13 / 1.73 | +0.9043 | +0.0334 / +0.0334 | -0.3756 (B, 143) | +1.4907 (A, 144) | +1.1151 | 0.1470 / 0.7572 (0.0837, 0.1762, 0.4975) | 1135.510 | 0.0109 |

### Faits bruts lisibles dans les tableaux (sans interprétation)

- 14/14 jobs COMPLETED, « bfgs converged » partout, 0 « convergence NOT achieved », aucune relance. Forces totales finales 2.3e-4 à 5.9e-4 Ry/bohr, composante max ≤ 9.8e-5 Ry/bohr.
- Hors plan : |dz| max ≤ 1.2e-6 Å partout sauf 8×8 nspin1 (1.2e-5 Å) et 11×11 nspin1 (7.1e-5 Å).
- nspin1 : d(paire) finale = 2.53 (5), 2.25 (6), 2.50 (7), 2.45 (8), 2.19 (9), 2.48 (10), 2.41 (11), 2.16 (12) Å ; spglib final P-6m2 (12 op.) pour 5×5 et 7×7, Amm2 (4 op.) pour les autres.
- nspin2 : d(paire) finale = 2.20 (5), 2.07 (6), 2.11 (7), 2.07 (8), 1.99 (9), 2.05 (10), 2.02 (11), 1.95 (12) Å ; Amm2 (4 op.) partout ; d(paire–3e) 2.568–2.592 Å.
- ΔE_relax nspin1 : −340 à −380 meV (8 tailles) ; ΔE_relax nspin2 : −516 à −862 meV ; ΔE_spin : −163 (5), −481 (6), −218 (7), −312 (8), −457 (9), −309 (10), −349 (11), −458 (12) meV.
- m_tot pw.x : 2.00 / 1.97 / 2.00 / 1.94 / 1.35 / 1.94 / 1.77 / 1.13 µB pour N = 5…12 ; m_abs 6×6 = 5.27 µB (2.45–2.89 ailleurs, 1.73 pour 12×12) ; sommes de Löwdin par sous-réseau du 6×6 : −2.19 / +4.15 µB (−0.38 à −0.73 / +1.49 à +2.63 ailleurs).
- Löwdin, atome 3e : +0.90 à +1.02 µB ; paire : +0.03 (12), +0.04 (9), +0.12 (6, 11), +0.16 (8, 10), +0.19 (5), +0.20 (7) µB. Composantes p du 3e : (px, py) ≈ (0.65, 0.02) pour lacune sur A et ≈ (0.18, 0.49) pour lacune sur B (orientation de la paire différente, cf. `perturbation.log`).
- Durées nspin2 / nspin1 : 2.6 (5), 2.3 (6), 3.0 (7), 3.1 (8), 3.1 (9), 3.6 (10), 2.8 (11), 2.2 (12).

## .save finaux (2026-09-24)

Scratch : `qe_tmp/vacancy_relaxed/series/NxN/nspinN/vac_NxN_relax_nspinN.save` (+ `.xml`) : nspin1 253 M (5) → 6.3 G (12), nspin2 503 M (5) → 13 G (12) ;
les `.save` nspin2 contiennent en plus `atomic_proj.xml` écrit par projwfc.x (2.6–70 Mo). Total 51 G, 84 fichiers.
Miroir /project (règle du 2026-09-17, rsync `-a -r --no-o --no-g --open-noatime`, atimes du scratch intacts) :
`graphene/qe/qe_tmp_backup/vacancy_relaxed/series/` = 84 fichiers, 54.45 Go, `MD5SUMS_series_2026-09-24.txt` (md5 calculés avec O_NOATIME des deux côtés,
84/84 identiques). Note : la première passe rsync a recréé les dossiers sans bit setgid (groupe `gregb26`) ; corrigé par `chgrp -R rrg-cotemich-ac` +
`chmod g+s` sur les dossiers ; seconde passe avec `--chmod=Dg+s`. Rien supprimé, ni sur scratch ni ailleurs.

## Fichiers du répertoire

`make_inputs_series.py`, `R2_phase0_table.txt`, `R2_jobs.md`, `analyze_relax_series.py`, `analyse_R2_2026-09-24.txt`, `R2_rapport.md`,
`NxN/perturbation.log`, `NxN/nspin{1,2}/{relax.in,relax.out,submit.relax,JOBID,slurm-*.out,slurm-*.err}`,
`NxN/nspin2/projwfc/{projwfc.in,submit.projwfc,JOBID,projwfc.out,slurm-*,pdos_*,proj_*}`.
Copie versionnée filtrée : `ab-initio-defects/article/R2_size_series/` (README.md). Aucun autre post-traitement (pas de ΔV, pw2wannier90, Wannier90 ni pipeline e-d). Rien n'est commité.
