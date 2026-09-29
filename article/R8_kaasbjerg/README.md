# R8 — Reproduction des figures 13 et 14 de Kaasbjerg, PRB 101, 045433 (2020)

- But : DOS (c_i = 0,1 % et 1 %, Fig. 13) et fonction spectrale Γ–K–M (c_i = 1 %, Fig. 14) du graphène avec lacunes sur un seul
  sous-réseau, avec M2 de production et la matrice t locale (bloc π), protocole de Kaasbjerg : Σ_k = c_i T̄_kk, c_i par MAILLE.
- Prompt : R8 (Greg, 2026-09-27) + ajout du 27 sept. + amendement du 28 sept. après R9 (N_k^int 900, variantes tel quel / aligné plateau, 2 bis, 7a avant les DOS).
- Statut : **TEST** (post-traitement seul ; résultats consignés dans `R8_rapport.md` et `results/M2/R8/`) — à confirmer par Greg au GO.
- Date : 2026-09-27 (phase 0) ; mise à jour du 2026-09-28 (amendement après R9).
- Lecture seule : `results/M2/`, `wannier/27x27`, `wannier/24x24`, json de R9 (jamais `R9_controles/cache/`) ; g₀ recalculé. Aucun calcul QE. Git en lecture.
- Copie versionnée : `article/R8_kaasbjerg/` (rapport, pilote, tables, figures ; jamais npz > 5 Mo ni slurm).
- État : phase 0 + mise à jour du 28 sept. + 7a (extraction Kaasbjerg) + fonctions Q1–Q6 écrites (non commitées), STOP ; attente : audit de l’image minimale, GO des calculs (J1…).
