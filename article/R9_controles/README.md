# R9 — Contrôles avant le chapitre 4

- But : (A) décalage uniforme de ΔV^L dans M2 (alignement de potentiel), (B) position de la résonance quand N_k^int converge,
  (C) chaîne repliée N×N contre la série QE de R7 et limite diluée, (D) grandeur exacte comparée aux ~70 eV Å² de Kaasbjerg (Fig. 3).
- Prompt : R9 (Greg, 2026-09-27). Ordre : phase 0 → STOP → étape G (Greg) → (GO) A → B → C → D → rapport → STOP.
- Statut : **TEST** (post-traitement seul ; aucun calcul QE ; résultats consignés dans `R9_rapport.md`).
- Date : 2026-09-27 (phase 0).
- Lecture seule : `results/M2/`, `results/M/`, caches Wannier, `.save`, `Vks_*`, `config/production.json`. Rien n'est écrit dans `results/`.
- Pilote unique prévu : `r9_driver.py` (ce répertoire). Copie versionnée : `article/R9_controles/` (rapport, pilote, tables, figures ;
  jamais les npz > 5 Mo ni les slurm).
- État : TERMINÉ le 2026-09-28 (phase 0, étape G, A, B, C, D ; rapport complet `R9_rapport.md`), STOP.
