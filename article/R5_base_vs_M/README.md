# R5 — Base ou M ? Diagnostic de la marche (a1) de R4 ; vérité DFT selon la taille (préparation article, hors mémoire)

- But : savoir si l'écart entre H_p + M en ondes planes (16 bandes × 81 k : état π à −1,350 eV, pas de doublet σ) et QE
  (−0,737 et +0,101 eV) vient de la troncature de la base repliée (bandes) ou de la matrice M elle-même, et où la DFT
  place l'état quasi-lié quand la taille de super-cellule varie (N = 5…12).
- Prompt d'origine : R5 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) A → C → (GO nscf) B → rapport → STOP.
- Statut : **TEST** (post-traitement ; un seul calcul QE autorisé, le nscf 128 bandes de B, sur GO séparé, dans un outdir neuf ;
  données de production en lecture seule ; routines de production intouchées ; git en lecture seulement).
- Dates : phase 0 le 2026-09-25 ; GO reçu le 2026-09-25 ; jobs J1–J5 exécutés et rapport terminé le 2026-09-25 (STOP).
- Répertoire de travail (celui-ci, hors dépôt, règle 5 de CLAUDE.md) : `graphene/qe/defects/R5_base_vs_M/`, à côté de
  `R4_quasi_lie/` qu'il prolonge. Copie versionnée prévue : `graphene-raman/article/R5_base_vs_M/` (rapport, pilote, tables,
  figures ; jamais les `.npy` de M, les `slurm-*`, `JOBID`, ni les fichiers > 5 Mo).
- Contenu : `R5_rapport.md` (phase 0, méthodes, A, C, B, fichiers, manifestes), modules de campagne `r5_sc_projection.py`,
  `r5_deltav_pw.py`, `r5_basis_diagnostics.py`, `r5_alignment_ext.py` (Greg a refusé l'écriture dans `src/`), pilote `r5_driver.py`,
  `submit_r5.sh`, `submit_r5_m128.sh`, `b/nscf_nb128.in` + `.diff` + `submit_nscf_nb128.sh` + `nscf_nb128.out`, résultats `a/ c/ b/`
  (json, npz, tables ; `b/M_*_nb128.npy` 3 × 1,72 Go non versionnés), `fig/`, `r5_log.txt`, `JOBID`, `slurm-*`.
- Jobs : J3 nscf 21808943 ; J4 M 128 bandes 21810702 ; J1 A 21811069 (tables par `atables`) ; J2 C 21811223 ; J5 B 21810855.
- Constat principal (A.2) : M^L de production en norme super-cellule, M^NL en norme maille (facteur 81) ; production non touchée.
- Aucune suppression ; tout nettoyage passe par un manifeste et un GO séparé.
