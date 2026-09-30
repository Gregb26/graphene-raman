# campagnes/M/ch4 — chiffres du chapitre 4 : table v1 → final, régions d'alignement, NOTES_TGAMMA au final

- But : source unique pour la réécriture du ch. 4 (`memoire/chapitres/défauts.tex` du dépôt du mémoire) avec la production finale (v2 + alignement de Kumagai–Oba,
  `results/M2_plateau/`, R10). Prompt de Greg du 2026-09-30 (« Chiffres du ch. 4 — table v1 → final, NOTES_TGAMMA au final, régions d'alignement »),
  GO du même jour avec les décisions D1–D5 par défaut et trois précisions (1.b, 1.g, 3).
- Statut : PRODUCTION (calcul local de quelques secondes, aucun job, aucun calcul QE ni t/Γ ; pas de répertoire de travail hors dépôt).
- Date : 2026-09-30 ; HEAD au moment de l'écriture : 66602ed (rien n'est commité par Code ; Greg relit et commit).
- Nomenclature (donnée une fois dans `table_v1_final.md`) : non aligné (= tel quel, `results/M2`) ; alignement à site unique (= Lu) ;
  alignement de Kumagai–Oba (= plateau (i), ≥ 0,75 r_max) ; ΔV_PA^(N) = C_N = `alignment.C_N_eV` de `config/production.json` ; final = v2 + Kumagai–Oba.

| fichier | rôle | sources |
|---|---|---|
| `ch4_chiffres.py` | pilote unique (sous-commandes `check`, `table`, `regions`, `notes`, `tex`, `all` ; docstring) ; racine du dépôt déduite de `__file__` | — |
| `table_v1_final.md` | partie 1 : table principale v1 → non aligné → final (207 lignes, statut mécanique), écarts entre les deux tables sources, chiffres nouveaux a–l, grandeurs retirées ; partie 5 : section « Compléments » (5.1 tab:rcut_M en %, 5.2 tableau niveau 1, 5.3 tab:échantillonnage, 5.4 scalaires lus sur les courbes des npz avec porte v1/v2) ; compléments lus dans `results/M`, `results/M2`, `results/M2_plateau` (`m_rcut_convergence.csv`, `level1_summary.csv`, `sampling_table.csv`, `resonance_9x9.npz`, `resonance_criteria_9x9.npz`) | `campagnes/R/R6_production_corrigee/etape3/table_v1_v2.md`, `campagnes/R/R10_plateau/c/table_v2_plateau.md`, json/npz/csv de R5, R6, R8, R9, R10, `results/M2`, `results/M2_plateau`, `config/production.json` |
| `alignement_regions.md`, `.npz`, `.pdf`, `.png` | partie 2 : régions (A) ≥ 0,75 r_max, (B) Kumagai–Oba 2D ≥ N·a/2, (C) Kumagai–Oba 3D min(N·a/2, c/2), (D) site unique, 13 tailles ; figure au style du mémoire, **non installée dans `figures/`** | `campagnes/R/R10_plateau/a/profiles_<S>.npz`, `a1_results.json` ; `scf.in` des super-cellules (a, c) ; `wannier/27x27/wannier.wout` (contrôle) |
| `NOTES_TGAMMA_partie5.diff` | partie 5 : `git diff` des trois lignes de `NOTES_TGAMMA.md` complétées par 5.4 (Born/T min et max, Γ_T à c = 0,1 %, ħ/Γ ; §2 et C13) | `table_v1_final.md`, compléments 5.4 |
| `NOTES_TGAMMA.diff` | partie 3 : `git diff` de `NOTES_TGAMMA.md` (alors à la racine du dépôt, dans `notes/` depuis le 2026-09-30 ; modifié en place : valeurs finales en tête, [v1 : …], ligne « Alignement du potentiel », C14, §6e, §8 R10) | `results/M2_plateau/*.csv`, `*.npz`, tables R10, R8, R5 |
| `defauts_nombres.md` | parties 4 et 5.5 : les 606 nombres de `défauts.tex` (lu dans le dépôt du mémoire, lecture seule), appariés à la table principale et aux compléments (pourcentages × 100, notation scientifique lue en entier, entiers < 100 comparés dans les `tabular`) ; liste des non appariés restants | `msc-graphene-raman-defects/memoire/chapitres/défauts.tex` (HEAD du mémoire et md5 dans l'en-tête) |
| `CLAUDE_md.diff` | diff **proposé, non appliqué** de `CLAUDE.md` (results_dir, sampling_table, ligne alignement, ligne de classement pour `campagnes/M/ch4/`) | — |
| `ch4_rapport.md` | rapport de session : phase 0 → 4, comptes, portes, écarts entre sources, erreurs d'archives, diffs | — |

- Relancer : `.venv/bin/python campagnes/M/ch4/ch4_chiffres.py all` (depuis la racine, `PYTHONPATH=src` ; `regions` lit les `scf.in` hors dépôt,
  chemins tirés d'`a1_results.json` ; `all` réécrit aussi la date et le HEAD des en-têtes ;
  `tex` lit le chapitre dans le dépôt du mémoire, voisin de celui-ci, sinon `$MSC_THESIS` ou `--tex CHEMIN` : la copie `memoire/défauts.tex` a été retirée le 2026-09-30, `6340c2f`). `check` est exécuté avant chaque sous-commande (porte de la partie 2 : les 13 C_N de la config au bit).
- Règles : tout chiffre copié d'un fichier source (chemin et clé/ligne dans la table) ; seuls calculs : partie 2 et les colonnes arithmétiques
  demandées (différences, produits) ; aucun arrondi (format de la source) ; deux sources différentes → les deux ; erreur d'archive → rapportée.
- Lecture seule : `results/`, `article/`, `config/`, `src/`, `scripts/`, `figures/`, et le dépôt du mémoire (`défauts.tex`). Aucune suppression.
