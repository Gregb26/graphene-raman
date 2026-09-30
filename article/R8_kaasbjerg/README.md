# R8 — Reproduction des figures 13 et 14 de Kaasbjerg, PRB 101, 045433 (2020)

- But : DOS (c_i = 0,1 % et 1 %, Fig. 13) et fonction spectrale Γ–K–M (c_i = 1 %, Fig. 14) du graphène avec lacunes sur un seul
  sous-réseau, avec M2 de production et la matrice t locale (bloc π), protocole de Kaasbjerg : Σ_k = c_i T̄_kk, c_i par MAILLE.
- Prompt : R8 (Greg, 2026-09-27) + ajout du 27 sept. + amendement du 28 sept. après R9 (N_k^int 900, variantes tel quel / aligné plateau, 2 bis, 7a avant les DOS).
- Statut : **TEST** (post-traitement seul ; résultats consignés dans `R8_rapport.md`, `out/` et `fig/` ; rien dans `results/`).
- Date : 2026-09-27 (phase 0) ; mise à jour du 2026-09-28 (amendement après R9).
- Lecture seule : `results/M2/`, `wannier/27x27`, `wannier/24x24`, json de R9 (jamais `R9_controles/cache/`) ; g₀ recalculé. Aucun calcul QE. Git en lecture.
- Copie versionnée : `article/R8_kaasbjerg/` (rapport, pilote, tables, figures ; jamais npz > 5 Mo ni slurm).
- État : 1er GO (étapes 1, 2, 2 bis, 3, 5, 7b) et GO article (étapes 4, 6) exécutés le 2026-09-29 ; résultats dans `R8_rapport.md`, `out/`, `fig/` ; STOP.
- Portée du mémoire (2026-09-30) : étapes 1–5 et 7 ; figures `dos_c`, `spectral_GKM`, `sigma_K`, `sensibilites`, `superposition` ; l'étape 6 (Dirac) reste pour l'article.
