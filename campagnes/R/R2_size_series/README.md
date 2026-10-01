# R2 — série en taille 5×5 → 12×12 de la lacune relaxée (préparation article, hors mémoire)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les scripts `make_inputs_series.py`, `analyze_relax_series.py` (aucune figure). Gardés : inputs, `relax.out`, rapports, analyses `.txt`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R2_size_series/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R2 — Série en taille 5×5 → 12×12 de la lacune relaxée
- **Question** : comment évoluent avec la taille N = 5…12 (9×9 = R1) la distorsion de la paire de voisins, l'énergie de relaxation ΔE_relax, l'écart ΔE_spin = F(nspin2) − F(nspin1) et le moment magnétique de la lacune relaxée ?
- **Pourquoi** : préparation de l'article, hors mémoire ; on étend le protocole R1 aux tailles du ch. 4 (`super_cell/NxN/defective/scf.in`). Motivation plus précise : non consignée.
- **Méthode** :
  - Protocole R1 inchangé : pw.x v7.5 `relax`, cellule fixe, bfgs, nstep 200, forc_conv_thr 1e-4, etot_conv_thr 1e-5, nosym/noinv. Le reste des paramètres vient du scf.in du ch. 4 de chaque taille (nbnd 130 à 604).
  - Spin : nspin1 et nspin2 (C1 sur les 3 voisins, starting_magnetization 0.5). Paire rapprochée de 0.03 Å par atome.
  - Jobs : 14 relax, 21666582–21666595, soumis le 2026-09-23 (Rorqual, 64 tâches, 1 j pour 5–8 et 3 j pour 10–12). Le 9×9 est repris de R1. projwfc.x : 21735568–21735574.
- **Résultat** :
  - 14/14 COMPLETED, « bfgs converged » partout, aucune relance. Forces finales 2.3e-4 à 5.9e-4 Ry/bohr ; durées de 00:02:47 à 07:54:45.
  - d(paire) finale, nspin2 : 2.20 (5), 2.07 (6), 2.11 (7), 2.07 (8), 1.99 (9), 2.05 (10), 2.02 (11), 1.95 (12) Å, symétrie Amm2 partout.
  - d(paire) finale, nspin1 : 2.53, 2.25, 2.50, 2.45, 2.19, 2.48, 2.41, 2.16 Å. Symétrie P-6m2 (12 op.) pour 5×5 et 7×7.
  - ΔE_spin : −163 (5), −481 (6), −218 (7), −312 (8), −457 (9), −309 (10), −349 (11), −458 (12) meV. ΔE_relax : −340 à −380 meV en nspin1, −516 à −862 meV en nspin2.
  - m_tot pw.x : 2.00 / 1.97 / 2.00 / 1.94 / 1.35 / 1.94 / 1.77 / 1.13 µB (N = 5…12). Moment de Löwdin du 3e voisin : +0.90 à +1.02 µB. Le 6×6 sort du lot : m_abs 5.27 µB.
- **Interprétation / décision** :
  - Interprétation : non consignée. R2_rapport ne donne que des « faits bruts … (sans interprétation) ».
  - Décisions consignées : campagne classée PRODUCTION le 2026-09-24 ; `.save` miroités, md5 84/84 identiques.
  - Suite : la 12×12 nspin1 sert de géométrie source au transplant R7c (`campagnes/R/R7_tailles_3m/README.md`, 2026-09-26).
- **Où sont les données** :
  - Tables : `R2_rapport.md` (tableau récapitulatif, ΔE_spin et Löwdin), `R2_jobs.md`, `R2_phase0_table.txt`, `analyse_R2_2026-09-24.txt`.
  - Entrées/sorties, pour NxN ∈ {5x5, 6x6, 7x7, 8x8, 10x10, 11x11, 12x12} : `NxN/perturbation.log`, `NxN/nspin1/{relax.in,relax.out,submit.relax}`, `NxN/nspin2/{relax.in,relax.out,submit.relax}`, `NxN/nspin2/projwfc/{projwfc.in,submit.projwfc}`.
  - Pas de `9x9/` (voir `../R1_vacancy_relaxed/`). Aucun json/npz ni `projwfc.out` : les valeurs de Löwdin ne figurent que dans le rapport et l'analyse.
  - `.save` hors dépôt : scratch et miroir `qe_tmp_backup/vacancy_relaxed/series/` (`MD5SUMS_series_2026-09-24.txt`).

---

Copie versionnée (règle 5 de CLAUDE.md, § « Campagnes de calcul ») du répertoire de travail
`graphene/qe/defects/super_cell_relaxed/series/` (hors dépôt, à côté de celui-ci). Statut : PRODUCTION.
Versionné ici : `make_inputs_series.py`, `analyze_relax_series.py`, `R2_phase0_table.txt`, `R2_jobs.md`, `R2_rapport.md`,
`analyse_R2_2026-09-24.txt`, `NxN/perturbation.log`, `NxN/nspin{1,2}/{relax.in,relax.out,submit.relax}` (relax.out de pw.x ≤ 2.3 Mo),
`NxN/nspin2/projwfc/{projwfc.in,submit.projwfc}`.
Non versionné : `JOBID`, `slurm-*`, `projwfc.out`, `pdos_*`, `proj_*`, et les `.save` (scratch `qe_tmp/vacancy_relaxed/series/NxN/nspinN/`,
miroir `graphene/qe/qe_tmp_backup/vacancy_relaxed/series/` avec `MD5SUMS_series_2026-09-24.txt`).
Le 9×9 de la série est celui de R1 (`campagnes/R/R1_vacancy_relaxed/`, répertoire de travail `super_cell_relaxed/9x9/`).
Les scripts portent des chemins absolus vers `graphene/qe/` et s'exécutent depuis le répertoire de travail.
