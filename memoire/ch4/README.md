# memoire/ch4 — chiffres du chapitre 4 : table v1 → final, régions d'alignement, NOTES_TGAMMA au final

- But : source unique pour la réécriture du ch. 4 (`memoire/défauts.tex`) avec la production finale (v2 + alignement de Kumagai–Oba,
  `results/M2_plateau/`, R10). Prompt de Greg du 2026-09-30 (« Chiffres du ch. 4 — table v1 → final, NOTES_TGAMMA au final, régions d'alignement »),
  GO du même jour avec les décisions D1–D5 par défaut et trois précisions (1.b, 1.g, 3).
- Statut : PRODUCTION (calcul local de quelques secondes, aucun job, aucun calcul QE ni t/Γ ; pas de répertoire de travail hors dépôt).
- Date : 2026-09-30 ; HEAD au moment de l'écriture : 66602ed (rien n'est commité par Code ; Greg relit et commit).
- Nomenclature (donnée une fois dans `table_v1_final.md`) : non aligné (= tel quel, `results/M2`) ; alignement à site unique (= Lu) ;
  alignement de Kumagai–Oba (= plateau (i), ≥ 0,75 r_max) ; ΔV_PA^(N) = C_N = `alignment.C_N_eV` de `config/production.json` ; final = v2 + Kumagai–Oba.

| fichier | rôle | sources |
|---|---|---|
| `ch4_chiffres.py` | pilote unique (sous-commandes `check`, `table`, `regions`, `notes`, `tex`, `all` ; docstring) ; racine du dépôt déduite de `__file__` | — |
| `table_v1_final.md` | partie 1 : table principale v1 → non aligné → final (207 lignes, statut mécanique), écarts entre les deux tables sources, chiffres nouveaux a–l, grandeurs retirées | `article/R6_production_corrigee/etape3/table_v1_v2.md`, `article/R10_plateau/c/table_v2_plateau.md`, json/npz/csv de R5, R6, R8, R9, R10, `results/M2`, `results/M2_plateau`, `config/production.json` |
| `alignement_regions.md`, `.npz`, `.pdf`, `.png` | partie 2 : régions (A) ≥ 0,75 r_max, (B) Kumagai–Oba 2D ≥ N·a/2, (C) Kumagai–Oba 3D min(N·a/2, c/2), (D) site unique, 13 tailles ; figure au style du mémoire, **non installée dans `figures/`** | `article/R10_plateau/a/profiles_<S>.npz`, `a1_results.json` ; `scf.in` des super-cellules (a, c) ; `wannier/27x27/wannier.wout` (contrôle) |
| `NOTES_TGAMMA.diff` | partie 3 : `git diff` de `NOTES_TGAMMA.md` (racine du dépôt, modifié en place : valeurs finales en tête, [v1 : …], ligne « Alignement du potentiel », C14, §6e, §8 R10) | `results/M2_plateau/*.csv`, `*.npz`, tables R10, R8, R5 |
| `defauts_nombres.md` | partie 4 : les 553 nombres de `memoire/défauts.tex` (lecture seule, non modifié), appariés à la table principale | `memoire/défauts.tex` (md5 dans l'en-tête) |
| `CLAUDE_md.diff` | diff **proposé, non appliqué** de `CLAUDE.md` (results_dir, sampling_table, ligne alignement, ligne de classement pour `memoire/ch4/`) | — |
| `ch4_rapport.md` | rapport de session : phase 0 → 4, comptes, portes, écarts entre sources, erreurs d'archives, diffs | — |

- Relancer : `.venv/bin/python memoire/ch4/ch4_chiffres.py all` (depuis la racine, `PYTHONPATH=src` ; `regions` lit les `scf.in` hors dépôt,
  chemins tirés d'`a1_results.json`). `check` est exécuté avant chaque sous-commande (porte de la partie 2 : les 13 C_N de la config au bit).
- Règles : tout chiffre copié d'un fichier source (chemin et clé/ligne dans la table) ; seuls calculs : partie 2 et les colonnes arithmétiques
  demandées (différences, produits) ; aucun arrondi (format de la source) ; deux sources différentes → les deux ; erreur d'archive → rapportée.
- Lecture seule : `results/`, `article/`, `config/`, `src/`, `scripts/`, `figures/`, `memoire/défauts.tex`. Aucune suppression.
