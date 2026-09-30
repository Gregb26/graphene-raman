# R1 / R1b — relaxation de la supercellule 9x9 avec lacune (préparation article, hors mémoire)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les scripts `make_inputs.py`, `analyze_relax.py`, `k3x3/make_k3x3.py`, `k3x3/analyze_k3x3.py` (générateurs d'inputs et analyses ponctuelles ; aucune figure). Gardés : inputs, sorties pw.x, rapports, analyses `.txt`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R1_vacancy_relaxed/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R1 — Relaxation de la lacune 9×9 (Γ, nspin 1/2) et contrôles 3×3×1 (R1b)
- **Question** : quelle géométrie relaxée (cellule fixe) prend la lacune simple de la supercellule 9×9 (161 atomes, lacune sur le sous-réseau A, voisins 64/80/81), avec et sans spin ? Quel état est le plus bas ? Ces résultats obtenus à Γ tiennent-ils avec une grille 3×3×1 (R1b) ?
- **Pourquoi** : préparation de l'article, hors mémoire. Le ch. 4 utilisait les positions idéales (`super_cell/9x9/defective/scf.in`). Motivation plus précise : non consignée.
- **Méthode** :
  - Calcul : pw.x `relax` (BFGS, cellule fixe, Γ) avec forc_conv_thr 1e-4, etot_conv_thr 1e-5, ecutwfc 100 Ry, mv 0.01 Ry, conv_thr 1e-10, nosym/noinv, nbnd 352, assume_isolated 2D. En nspin2, espèce C1 sur 64/80/81 avec starting_magnetization 0.5. Géométrie de départ : 64 et 80 rapprochés de 0.03 Å chacun.
  - Jobs : 21587934 (nspin1, 1 h 03 min 37 s) et 21587935 (nspin2, 3 h 18 min 36 s).
  - R1b : scf 3×3×1 (9 points k) sur les géométries idéale et relaxées, jobs 21600427–21600430 ; projwfc.x (Löwdin), jobs 21600431–21600432.
- **Résultat** :
  - E(nspin2) − E(nspin1) = −0.03362052 Ry = **−457.43 meV** à Γ ; −389.68 meV à 3×3×1 (relax2_nspin2 − relax1_nspin1, avec F).
  - d(64-80) : 2.46585 → 2.19194 Å (nspin1) et → 1.98798 Å (nspin2). d(64-81) = d(80-81) = 2.57319 Å (nspin1) et 2.57260 Å (nspin2). Hors plan max 5.2e-7 Å. Symétrie Amm2 (n° 38), site mm2.
  - nspin2 : m_tot / m_abs = +1.35 / 2.45 µB. Moment de Löwdin de l'atome 81 : +0.9111 µB (Γ), +0.9408 µB (3×3).
  - Énergie de relaxation ΔE_relax, nspin1, par rapport au scf idéal du ch. 4 : −365.97 meV à Γ ; −298.86 meV à 3×3.
  - Forces à 3×3 sur les géométries relaxées à Γ : 0.049 Ry/bohr (nspin1) et 0.043 Ry/bohr (nspin2). ideal_nspin2 (3×3) a été annulé avant convergence (6.0e-10 Ry).
- **Interprétation / décision** :
  - Interprétation : non consignée dans R1/R1b. Seul constat écrit : « la géométrie bougerait sous 3×3 » (R1b).
  - Décisions consignées : ideal_nspin2 annulé par Greg, pas de relance. R1b classé « TEST consigné » ; ses `.save` ont été supprimés le 2026-09-23 (`admin/CLEANUP.md`, étage 7).
  - Suites : les géométries R1 sont figées et réutilisées comme 9×9 de R2 et pour R4 J0 (pp.x).
- **Où sont les données** :
  - Tables : `R1_rapport.md`, `R1b_rapport.md`, `analyse_2026-09-22.txt`, `k3x3/analyse_k3x3_2026-09-22.txt`.
  - Entrées/sorties : `perturbation.log`, `nspin1/{relax.in,relax.out,submit.relax}`, `nspin2/{relax.in,relax.out,submit.relax}`, `k3x3/{ideal_nspin1,ideal_nspin2,relax1_nspin1,relax2_nspin2}/{scf.in,scf.out,submit.scf}`, `k3x3/JOBIDS`, `projwfc/{gamma_nspin2,k3x3_relax2_nspin2}/{projwfc.in,submit.projwfc}`.
  - Aucun json/npz. Les `pp.in`, `pp_up.in`, `pp_dw.in` et `submit.pp` cités plus bas (R4 J0) sont absents de la copie.
  - `.save` hors dépôt : scratch et miroir `qe_tmp_backup/vacancy_relaxed/nspin{1,2}/`.

---

Copie versionnée des inputs, scripts et rapports de `graphene/qe/defects/super_cell_relaxed/9x9/` (déplacé le 2026-09-23 depuis `graphene/qe/vacancy_relaxed/`) (hors dépôt, à côté de celui-ci, même
convention que `graphene/qe/epw/` pour le chapitre 5). Les `relax.out`/`scf.out` de pw.x (< 5 Mo) sont versionnés depuis le 2026-09-23 ; les autres sorties (`projwfc.out`, `slurm-*`, `pdos_*`,
`proj_*`) et les `.save` ne sont pas versionnés : `.save` de R1 sur scratch `qe_tmp/vacancy_relaxed/nspin{1,2}` et miroir
`graphene/qe/qe_tmp_backup/vacancy_relaxed/` (md5) ; `.save` de R1b (`k3x3/`) sur scratch seulement (tests reproductibles).
Les scripts `make_inputs.py`, `analyze_relax.py`, `k3x3/make_k3x3.py`, `k3x3/analyze_k3x3.py` portent des chemins absolus
vers `graphene/qe/` et s'exécutent depuis le répertoire d'origine.
- 2026-09-25 (R4, J0) : potentiels locaux pp.x (plot_num 1) des géométries relaxées, calculés pour la campagne R4 (`graphene/qe/defects/R4_quasi_lie/`) : `nspin1/Vks_R1_nspin1` (spin_component 0) et `nspin2/Vks_R1_nspin2_{up,dw}` (spin_component 1, 2), entrées `pp.in`, `pp_up.in`, `pp_dw.in`, `submit.pp`, jobs 21796633 / 21796634 ; fichiers Vks (241 Mo) et `pp*.out` (184 Mo) non versionnés.
