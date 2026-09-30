# Cartographie du dépôt `graphene-raman` — 2026-09-30 (ménage, étape 1, lecture seule)

Réécriture complète de ce fichier. L'ancien contenu (inventaire P12 du stockage Rorqual, 2026-09-16 : quotas, purge du
scratch, classement irremplaçable / reproductible / jetable, chemins codés en dur) est dépassé par CLEANUP.md (tous les
étages exécutés le 2026-09-23/24) et reste lisible dans git : `git show 193e25d:INVENTAIRE_2026-09-16.md`.

Session strictement en lecture seule : rien déplacé, renommé ni supprimé ; seule écriture, ce fichier. Aucun commit.
État cartographié : première passe sur HEAD `193e25d` (« updated plan », 2026-09-30 08:46) ; **complété le même jour sur
HEAD `656c05c`** (quatre commits de plus : `241739e`, `66602ed`, `90d2022`, `656c05c`), arbre propre. Les ajouts de la
seconde passe sont au §0 bis et reportés dans les sections concernées (§5 R8, §6.2, §8, §9, §10, §12, §13). Les numéros
de ligne de CLAUDE.md cités au §9 sont ceux de `193e25d` (décalés de +4 après la l. 243 depuis `90d2022`). **Troisième
passe (HEAD `7ac2c69`, 2026-09-30 11:45)** : seul l'emplacement et l'état du dépôt du mémoire ont été mis à jour (§0 bis,
§8, §10, §12, §13 étape F) ; rien d'autre n'a changé dans ce dépôt. Les calculs rorqual du mémoire sont terminés (R10
clos le 2026-09-30). **Quatrième passe (HEAD `e916500`, 2026-09-30 11:53, arbre propre, `main` = `origin/main`)** : le
dossier local de ce dépôt a été renommé `raman-graphene` → **`graphene-raman`** vers 11 h 46 (même nom que le remote
GitHub et que le clone rorqual) ; ce renommage **casse l'installation éditable du venv** (§0 bis, §11, §13 étape 0).
Aucun fichier suivi n'a changé depuis la troisième passe ; tout le reste a été revérifié tel quel (comptes de fichiers
suivis, 199 tests collectés avec `PYTHONPATH=src`, 20 figures du mémoire dont 8 diffèrent, md5 de `défauts.tex`).
**Cinquième passe, courte (HEAD `7446981`, 2026-09-30 12:31, arbre propre, `main` = `origin/main`)** : trois commits de
plus (`214ed02`, `6340c2f`, `7446981`). L'**étape 0 est faite** (installation éditable réparée par Greg, 199 tests
collectés sans `PYTHONPATH`) ; la campagne `campagnes/M/ch4/` a reçu sa **partie 5** (compléments, trois lignes de
`NOTES_TGAMMA.md`) ; **`memoire/défauts.tex` a été retiré du dépôt**, ce qui casse `ch4_chiffres.py tex` (§6.2, §11) ;
côté mémoire, `PROVENANCE.md` est commité et le PDF compilé n'est plus suivi (HEAD `4da4d62`). Détail au §0 bis ;
sections touchées : §6.2, §7, §8, §9, §11, §12, §13. Rien d'autre n'a changé (comptes de fichiers suivis, `results/`,
`scripts/`, `src/`, `article/`, `figures/` identiques ; les 8 mêmes figures diffèrent).

**Sixième passe = exécution (2026-09-30 après-midi, HEAD `32efd9e` « cleanup: untangled broken, obsolete and incoherent
code », qui contient la cinquième passe, la réparation de `ch4_chiffres.py tex` et l'arbitrage du §11).** L'**étape B est
faite** (non commitée) : voir le tableau « Étape B — exécutée » au §13. Conséquence pour la lecture de ce fichier : partout
ailleurs, `NOTES_TGAMMA.md`, `NOTES_EPW.md`, `NOTES_EPW_REPERES.md` désignent maintenant `notes/<nom>` ; `CLEANUP.md` et
`REECRITURE_HISTORIQUE_2026-09-29.md` sont dans `admin/` ; `reports/` et `research_notes/` sont sous
`article/C_optique_lacunes/biblio/` ; `requirements.txt` n'existe plus ; les numéros de ligne de CLAUDE.md cités au §9 sont
ceux d'avant sa correction.

**Déplacement du 2026-09-30 (après l'étape « campagnes ») : `article/R*` → `campagnes/R/R*`, `memoire/EM` → `campagnes/EM`,
`memoire/ch4` → `campagnes/M/ch4` ; `article/` ne garde que `C_optique_lacunes/`. Les chemins de ce fichier ont été réécrits
mécaniquement (y compris dans les descriptions de commits du §0 bis, qui parlaient d'`article/…` et `memoire/…` à l'époque).
Index des campagnes : `campagnes/README.md`.**

Méthode : lecture de chaque fichier de `scripts/`, `src/`, `tests/`, `article/`, `memoire/`, `results/`, `config/`,
`.md` racine ; résolution **par import réel** de tous les `from graphene_raman… import …` (635 imports,
script AST + `importlib`) et par AST des imports inter-scripts (`sys.path`) ; grep des figures et tables dans le dépôt
du mémoire (lecture seule ; aujourd'hui `~/projects/msc-graphene-raman-defects`, voir §0 bis) ; collecte pytest (`--collect-only`, 199 tests, 0 erreur).

Vocabulaire des statuts : **actif** = fait partie de la base vivante R10 ou d'un usage courant ; **figé** = campagne
terminée, gardé pour la provenance, encore exécutable ; **obsolète** = remplacé, non rejouable ou désactivé
(`SystemExit`) ; **cassé** = plante ou lit le mauvais répertoire à l'exécution. Campagnes : **chaîne M** (calcul de M),
**chaîne T** (t local, Γ), **P** (série mémoire), **R1–R10** (article), **EM** (électron-photon, §2.5), **EPW** (ch. 5).

---

## 0. Constats principaux

1. **Aucun import cassé.** `alignment_C` et `matrices_dir` existent dans `src/graphene_raman/config.py`
   (l. 30 et 35), ajoutés au commit `c2bc733` « checkpoint R10 » (2026-09-29 20:55), poussé sur `origin/main`. Les 635
   imports du paquet se résolvent, ainsi que tous les imports inter-scripts (`_palette`, `_bands`, `_paths`,
   `r4_driver`/`r5_driver`, `em2_B_compare`, `test_ks_reconstruction`). La ligne « ImportError en local, 2026-09-29 »
   de `article/C_optique_lacunes/PLAN.md` (dernière ligne) date du commit précédent `7414c2c` et est périmée. Seul nom
   inexistant : `qe_io.get_fermi` dans `validate_wannier_bands.py:67`, protégé par `hasattr`.
2. **La casse réelle est un problème de chemins, introduit par R10.** `config/production.json` sépare depuis R10
   `results_dir = results/M2_plateau` (produits npz/csv seulement) de `matrices_dir = results/M2` (les 43 matrices M2,
   présentes uniquement sur `/project` de rorqual ; **aucun `.npy` dans `results/M2` en local**). Les scripts qui
   lisent ou écrivent des `M_*.npy` via `results_dir()` visent donc un répertoire sans matrices. Liste exacte au §11.
3. **Base vivante** : `results/M2_plateau` (R10, `MD5SUMS_2026-09-30.txt`) ; `results/M2` = matrices v2 (R6) + produits
   gelés ; `results/M` = v1 gelé. CLAUDE.md (l. 243), NOTES_TGAMMA.md, PLAN C et `submit_post.sh` disent encore que
   `results/M2` est vivant.
4. **Documents périmés** : CLAUDE.md (scripts `run.py`/`compute_M_cluster.py` inexistants depuis `cf05e92`, carte des
   modules incomplète, `single_defect` dit orphelin, tableau des campagnes sans R2–R10), README.md (v1, `results/M`,
   `notebooks/`, EM absente), NOTES_TGAMMA.md (état R6, pas R10), requirements.txt (freeze de l'ère ABINIT), CLEANUP.md
   (en-tête « dry-run » alors que tout est exécuté). Détail au §9.
5. **Figures du mémoire** : le dépôt LaTeX inclut 20 figures ; 8 d'entre elles diffèrent (md5) de la version R10
   régénérée le 2026-09-30 dans `figures/` (le mémoire a encore la version R6). `fig_em_coupling` n'est pas encore
   incluse (§2.5 non rédigé). Détail au §8.
6. **Doublons** : `results/M/ved_analysis.npz` = `results/M2/ved_analysis.npz` (md5 identique, 13,5 Mio ×2) ;
   `ks_reconstruction.npz` identique entre M2 et M2_plateau ; 33 figures de `figures/` ont une copie octet pour octet
   dans `campagnes/R/R10_plateau/fig/` ou `campagnes/R/R6_production_corrigee/etape3/figures_v2/` ; scripts en doublon au §12.
7. **Deux systèmes de tests** : `tests/` (pytest, 199 tests, tous sans données externes, EM + R8/R9/R10 + matrix_io)
   et `scripts/test_*.py` (8 scripts PASS/FAIL, chaîne M/T, 6 exécutables en local). Détail au §7.

## 0 bis. Nouveautés depuis la première passe (commits `241739e` → `656c05c`, 2026-09-30 09:11 → 10:34)

| Commit | Contenu | Effet sur cette carte |
|---|---|---|
| `241739e` « added sigma_k figure to thesis » | `campagnes/R/R8_kaasbjerg/{R8_rapport.md,README.md}` : l'étape 4 (`fig/sigma_K`) entre dans la portée du mémoire (étapes 1–5 et 7) ; l'étape 6 (Dirac) reste pour l'article | §5 (R8), §10 ; aucune figure refaite ; rien n'est encore installé dans le mémoire |
| `66602ed` « added chap4 tex file » | `memoire/défauts.tex` (750 lignes) : copie du chapitre 4 du mémoire, état **v1** | §6.2 ; copie identique à la version du dépôt du mémoire à une ligne vide près |
| `90d2022` « conversion v1 -> final chap4 thesis done » | nouvelle campagne `campagnes/M/ch4/` (11 fichiers, pilote `ch4_chiffres.py`) ; `NOTES_TGAMMA.md` porté à l'état final R10 (en-tête, §1, §2, §3, §6e, §8) ; `CLAUDE.md` corrigé sur 4 points (`results_dir` = M2_plateau, `sampling_table.csv`, ligne « Alignement », ligne `campagnes/M/ch4/` du tableau des campagnes) | §6.2, §9 (NOTES_TGAMMA n'est plus « périmé depuis R10 » ; CLAUDE.md l. 243 et 255 corrigées) |
| `656c05c` « re-did inventory for eventual clean up » | ce fichier (première passe) | — |
| `7ac2c69` « updated inventory for cleanup » | ce fichier (seconde passe) | — |
| `e916500` « initialisation de PROVENANCE.md » | ce fichier (troisième passe) ; malgré son titre, ne touche que cet inventaire (le `PROVENANCE.md` est dans le dépôt du mémoire) | — |
| `214ed02` « finished all calculations on the cluster » (12:19) | `campagnes/M/ch4/` **partie 5** : `ch4_chiffres.py` (+187 l., fonctions `compl_5_1` … `compl_5_4`, appariement du tex refait), `table_v1_final.md` (section « Compléments »), `defauts_nombres.md` (refait), `ch4_rapport.md` (partie 5, annexe C), `README.md`, nouveau `NOTES_TGAMMA_partie5.diff` ; `NOTES_TGAMMA.md` : trois lignes (§2 ×2, C13). Aucun fichier de `results/`, `article/`, `scripts/`, `src/` ni `figures/` : le titre annonce la fin des calculs, le commit ne contient que des chiffres relus en local | §6.2, §9, §12 |
| `6340c2f` « removed thesis chapter » (12:21) | suppression de `memoire/défauts.tex` (la copie du ch. 4 ; l'original reste dans le dépôt du mémoire) | §6.2, §11 (`ch4_chiffres.py tex` cassé), §12 (doublon résolu), §13 étape D |
| `7446981` « inventory update » (12:22) | ce fichier (quatrième passe) | — |

Hors dépôt, le même matin : **le dépôt du mémoire a changé deux fois de place, et de forme.** Emplacement actuel :
**`~/projects/msc-graphene-raman-defects/`** (remote GitHub `Gregb26/msc-graphene-raman-defects`, HEAD **`4da4d62`**
à la cinquième passe, arbre propre ; `e07e0d5` aux troisième et quatrième). Des deux anciens emplacements, `~/LaTeX/master_thesis`
(première passe) existe encore, **vide** depuis 11 h 21 ; `~/LaTeX/msc-graphene-raman-defects` (seconde passe) **n'existe
plus**. Le dossier parent s'appelle `~/projects` (minuscule), et `~/projects/qe_pp` est le voisin de ce dépôt.

**Renommage local de ce dépôt (quatrième passe, vers 11 h 46)** : `~/projects/raman-graphene` → **`~/projects/graphene-raman`**.
Le nom local rejoint celui du remote (`Gregb26/graphene-raman`), du clone rorqual, de CLAUDE.md l. 46-49 et de
`PROVENANCE.md` ; `~/projects/raman-graphene` n'existait plus à 11 h 53 ; à la cinquième passe c'est un **dossier vide**,
recréé à 11 h 55 (ni dépôt ni fichier ; origine inconnue, probablement l'ancienne fenêtre de l'éditeur ouverte sur ce
chemin) : à supprimer par `rmdir` s'il n'est pas voulu. Effets vérifiés (état de la quatrième passe, puis de la cinquième) :

| Élément | État après renommage | Remarque |
|---|---|---|
| fichiers suivis | aucun chemin `raman-graphene` ni `/Users/gregou` (`git grep`, hors ce fichier) ; `config.ROOT` se déduit de `__file__` → `/Users/gregou/projects/graphene-raman` | rien à corriger dans le dépôt |
| **installation éditable** | **cassée** : `.venv/lib/python3.13/site-packages/__editable__.graphene_raman-0.1.0.pth` pointe vers `/Users/gregou/Projects/raman-graphene/src` (inexistant) ; `import graphene_raman` → `ModuleNotFoundError` ; `.venv/bin/python -m pytest tests` échoue au chargement de `conftest.py` | contournement : `PYTHONPATH=src` (199 tests collectés) ; **réparée à la cinquième passe** (étape 0 faite par Greg) : `.pth` et `direct_url.json` → `/Users/gregou/projects/graphene-raman`, import et collecte pytest (199 tests) sans `PYTHONPATH` |
| exécutables du venv | `.venv/bin/pytest` et `py.test` ont un shebang vers l'ancien chemin (morts, **encore à la cinquième passe** ; `pip` a été réécrit) ; `activate` et `pyvenv.cfg` citent un chemin encore plus ancien (`~/graphene_raman/.venv`) | `.venv/bin/python -m …` fonctionne ; venv jamais recréé depuis août 2025 |
| liens symboliques, `data/` | aucun lien dans l'arbre (hors `.venv`) ; rien de cassé | — |
| variables `PROJECTS`, `GRAPHENE_RAMAN` | non définies en local (elles le sont dans `~/.bashrc` de rorqual) ; le défaut `$PROJECTS/graphene-raman` des pilotes et d'EM2 correspond maintenant aussi au nom local | EM3 (`make_em3_data.py`) les pose lui-même depuis `__file__` |
| mémoire persistante de Claude Code | rangée par chemin : l'ancienne est sous `~/.claude/projects/-Users-gregou-Projects-raman-graphene/memory/` (7 notes), la nouvelle (`-Users-gregou-projects-graphene-raman`) était vide | notes recopiées et mises à jour à la quatrième passe ; les transcriptions des anciennes sessions restent sous l'ancien chemin |

État du dépôt du mémoire à `e07e0d5`, mis à jour à `4da4d62` (deux commits de 11:59-12:00, « fin du suivi git du pdf
compilé pour linstant » : `c844c22` retire `master_thesis.pdf` du suivi, `4da4d62` l'ajoute au `.gitignore` et commite
`PROVENANCE.md`) :

| Élément | Contenu | Remarque |
|---|---|---|
| `memoire/` | source LaTeX (`master_thesis.tex`, `chapitres/`, `annexes/`, `pages_liminaires/`, `tikz/`, `dms.cls`, `.latexmkrc`, `build/`) ; `figures/` réduit aux **20 PDF inclus** (les 19 figures non incluses ont été retirées) | `latexmk` depuis `memoire/`, PDF écrit à la racine (`master_thesis.pdf`, **ignoré par git depuis `4da4d62`**) |
| `analyse/`, `calculs/`, `donnees/`, `scripts_figures/` | dossiers d'accueil, **encore vides** (un README d'une ligne chacun) | à remplir : §13 étape F |
| `README.md` | structure, compilation, installation du venv ; « Scripts : à venir », « Données lourdes : à venir » | — |
| `PROVENANCE.md` (172 l. ; vide dans `dfc3180`, **commité dans `4da4d62`**) | registre de provenance : les 20 figures incluses (label, `.tex:ligne`, producteur et données dans ce dépôt, état R10 ou R6), les figures à venir, les 16 tables, les chiffres cités par chapitre, annexe ligne par ligne | **s'appuie sur ce fichier** (`INVENTAIRE_2026-09-16.md` §8, §10, §13 étape F, HEAD `656c05c`) ; sa colonne « Dans msc- » est vide et se remplira à la copie, avec le commit source |
| `requirements.txt`, `.venv/` (Python 3.12.1) | gel de 154 paquets d'un environnement généraliste (abipy, netCDF4, torch, jupyterlab, numpy 1.26.4, matplotlib 3.8.2) | **ne correspond pas** au venv de ce dépôt (Python 3.13, numpy 2.3.5) ; sans `mpi4py`, `pytest` ni le paquet `graphene_raman` : à refaire quand `analyse/` et `scripts_figures/` seront remplis |

Conséquences pour cette carte : §8 (chemin et contenu de `figures/`), §10 et §13 étape F (où copier quoi, et
articulation avec `PROVENANCE.md`). Le skill `thesis-section-pass` citait encore `~/LaTeX/master_thesis` ; **corrigé à la
cinquième passe** (description et SKILL.md l. 9 : `~/projects/msc-graphene-raman-defects`).

Ce que la seconde passe change dans les constats ci-dessus : le point 3 est corrigé dans CLAUDE.md et NOTES_TGAMMA.md
(reste README.md, PLAN C, `submit_post.sh`) ; le point 4 ne vaut plus pour NOTES_TGAMMA.md ; le point 5 est inchangé
(les 8 mêmes figures diffèrent, et le texte du ch. 4 cite encore les valeurs v1 : voir `campagnes/M/ch4/defauts_nombres.md`).

---

## 1. Configuration et répertoires de résultats

### 1.1 `config/production.json` (commit `c2bc733`, R10) et `src/…/config.py`

| Clé / fonction | Valeur / rôle | Lu par |
|---|---|---|
| `results_dir` = `results/M2_plateau` | produits vivants | `config.results_dir(cfg)` : 39 fichiers |
| `matrices_dir` = `results/M2` | matrices M2 brutes, lecture seule | `config.matrices_dir(cfg)` : `analyze_M`, `compute_spectral_wannier`, `mwr_locality_coarse_vs_dense` ; `dense_paths()["mfile"]` |
| `results_dir_frozen` = `[results/M, results/M2]` | liste depuis R10 | `r6_compare_v1_v2.py:14` seulement (qui attend une chaîne → TypeError) |
| `alignment` (plateau (i), 13 C_N en eV, source `R10_plateau/a/a1_results.json`, md5 vérifié) | décalage C_N soustrait sur M_W(R,R) | `config.alignment_C(cfg, size)` : 10 fichiers |
| `R_cut 3`, `grid 240`, `eta_eV 0.02`, `nk_int 300`, `e_window_eV 3`, `ne_per_eta 8`, `N_min 7`, `reference_size 9x9` | paramètres gelés le 2026-09-05 | chaîne T |
| `dense` : 5x5→25, 6x6→24, 7x7→28, 8x8→32, 9x9→27, 12x12→24 (p, D) | grilles denses | `config.dense_paths(cfg, size)` : 17 fichiers ; `KeyError` pour 10x10, 11x11, ≥15x15 |
| `M_normalization` = v2 (R6) ; `units` ; `K_red` ; `families` ; `defect_concentration_for_dos 0.01` | conventions | `matrix_io` (porte v2), `resonance_*` |
| non lues par le code | `labels_offgrid`, `production_sizes`, `nbnd_dense`, `M_normalization_note`, `K_note`, `source_md5`, `units` | documentaires |

`config.py` : `HA2EV` (dupliqué dans `matrix_io.py:26` et redéfini dans 8 scripts), `ROOT` (dérivé de `__file__`),
`load_production` (exige 12 clés), `results_dir`, `matrices_dir`, `alignment_C`, `dense_paths` (scratch rorqual codé
en dur `/home/gregb26/links/scratch/qe_tmp` : tout script qui l'appelle ne tourne que sur le cluster). Statut : actif.

### 1.2 Les trois `results/`

| Répertoire | Contenu | Écrit par | Statut |
|---|---|---|---|
| `results/M` | v1 (M^L en norme super-cellule, facteur N_cells manquant, R5-A.2) : 23 produits suivis ; `M_ed_{5..12}.npy` + `M_ed.npy` (juin, orphelin) ignorés par git ; **pas de README local** (le texte « GELÉ » n'est que dans `campagnes/R/R6_production_corrigee/phase0/README_results_M_gele.md`, la liste blanche `.gitignore` de `results/M/` exclut README.md) | chaîne M/T jusqu'au 2026-09-24 | gelé v1 ; encore lu par R4 (`RES`), R5, `r6_kernel_check`, `r6_m_rcut_resigma`, `r6_compare_v1_v2`, `make_inputs_series` (R2, sidecars), `submit_post.sh ksrec` (l. 31 copie `ved_analysis.npz` v1) |
| `results/M2` | v2 (R6) : README, `MD5SUMS_2026-09-25.txt` (43 M npy, **tous sur rorqual seulement**), 28 produits suivis | R6 (`runbook_3.sh`, `submit_post.sh`, `assemble_M2.py`, `r6_*`) | matrices vivantes (`matrices_dir`) ; produits gelés (« v2 tel quel ») |
| `results/M2_plateau` | R10 : README, `MD5SUMS_2026-09-30.txt` (27/33 vérifiables en local, 6 `resigma_*` absents et non versionnés), 28 produits suivis ; `resonance_*.npz` sans `ML_diag_mean`/`G_S`/`Gamma_T_noshift`, avec `C_N_eV` | R10 (`submit_r10.sh` → lanceurs de `scripts/`, `r10_driver c1post/c5`) | **vivant** (`results_dir`) |
| `results/epw` | 25 npz suivis (validation, selfen, phself, dfpt_path_freq, phdos, ring_check, d2_extract ; chaîne `_mv0.02` = production, sans suffixe = ancienne chaîne degauss 0.002) ; `ed_vs_ep_*.npz` = v1 obsolètes (l'écrivain écrit dans `results_dir` depuis R6) | `epw_*.py` | figé (ch. 5) |
| `results/test_recon` | `ks_recon_*.npz`, `summary.txt` (juin), ignoré | aucun script (l'écrivain `run_test_A_batch.py` vise `results/test_A`) | orphelin |
| fichiers isolés (ignorés) | `KB_projectors_C.pdf` (sans producteur ; inclus par le mémoire ch. 2), `compare_bands_{qe,w90_qe}.png`, `wannier_bands_path.png`, `wannier_vs_dft_coarse.png` | `compare_bands_*`, `validate_wannier_bands` (écrit dans le cwd) | figé |

Écrivains / lecteurs par produit (scripts de `scripts/` sauf mention) :

| Produit | Écrit par | Lu par |
|---|---|---|
| `M_analysis.npz`, `M_tests_summary.csv` | `analyze_M.py` (+ lignes de porte : `r6_tests_gate_row`, `r10 c1post`) | `make_figures*`, `level2_families`, R9, R10, `r6_compare` |
| `ved_analysis.npz` | `analyze_Ved.py` | `make_figures*`, `sampling_table`, R9, R10 |
| `mwr_locality.npz` | `mwr_locality_coarse_vs_dense.py` | `make_figures*`, `level2_families`, R9, R10 |
| `ks_reconstruction.npz` | `ks_reconstruction_all.py` (R10 : copié de M2) | `make_figures`, `sampling_table` |
| `level1_summary.csv`, `level2_summary.csv` | `make_figures.py` (effet de bord) | R10, `r6_compare` |
| `level2_families.csv` | `level2_families.py` | R10, `r6_compare` |
| `lnl_frobenius.csv` | `lnl_frobenius_all.py` | R10 (tab:L_NL) |
| `m_rcut_convergence.csv` | `m_rcut_convergence.py` (**en ajout** : une relance duplique les lignes) | R9, R10 (tab:rcut_M) |
| `m_rcut_resigma.csv` | `campagnes/R/R6…/etape3/r6_m_rcut_resigma.py` ou `r10 c1post` (**aucun script de `scripts/`**) | R10 |
| `nkint_check_9x9.csv` | `nkint_check_post.py` (lit `resigma_9x9_rc3_nk*.npz`, non versionnés) | R10 |
| `resonance_<S>.npz`, `resonance_9x9_shiftL.npz` | `resonance_metrics.py` (shiftL : R10 ±9,05 meV ; `submit_post.sh c14` a encore ±25) | `make_figures*`, `resonance_criteria`, `epw_ed_vs_ep`, R8, R9, R10 |
| `resonance_criteria_<S>.npz` | `resonance_criteria.py` | `make_figures*` |
| `specwd_<S>_prod.npz` | `compute_spectral_wannier.py` (via `submit_spectral_wannier_dense.sh`) | `level2_families`, `make_figures*`, `r6_level1_gate`, R9 |
| `sampling_table.csv` | `sampling_table.py` | CLAUDE.md, mémoire (tab:échantillonnage) |
| `ed_vs_ep_*.npz` | `epw_ed_vs_ep.py` | R10, `r6_compare` |
| `resigma_9x9_*.npz` (non versionnés) | `rcut_resigma.py` | `nkint_check_post`, `r6_level1_gate`, `r6_m_rcut_resigma`, `r10 c1post` |

### 1.3 `wannier/` (72 fichiers suivis, propre)

12 grilles × 6 fichiers (`wannier.eig`, `wannier.wout`, `wannier_manifest.json`, `wannier_tb.dat`, `wannier_u.mat`,
`wannier_u_dis.mat`). Denses de production (`dense_paths`) : 25x25 (5×5), 24x24 (6×6 et 12×12), 28x28 (7×7), 27x27
(9×9, toute la série EM, `tests/conftest.py` avec sha256), 32x32 (8×8). Grossières 5x5, 7x7, 8x8 :
`mwr_locality_coarse_vs_dense.py`, `test_local_tmatrix_real.py`. `*_nb16` (25, 27, 28, 32) : héritage nbnd 16,
**aucun lecteur** (candidats au ménage, 4 × ~4 Mo). Les `24x24` et `*_nb16` « non suivis » du 2026-09-16 le sont depuis.

---

## 2. `src/graphene_raman/` (42 fichiers suivis)

Paquet d'espace de noms (pas de `__init__.py` à la racine ni dans `io/`, `wannier/`, `wavefunctions/`, `plotting/`,
`many_body/`), installé en éditable (`.pth` → `src/` ; refait après le renommage du dossier, §0 bis). `mpi4py` est importé au chargement de `local_G`, `local_R`,
`non_local` (donc `deltav_pw`) et absent de `requirements.txt`.

| Module | Rôle | Importé par | Campagne | Statut |
|---|---|---|---|---|
| `config.py` | chargeur unique de `production.json`, chemins | 52 fichiers | R10 | actif |
| `io/qe_io.py` | lecteur QE (`.save` : XML Ha, wfc hdf5, filplot pp.x Ry→Ha) | 46 fichiers | chaîne M | actif (figé 2026-09-02) ; `get_typat` morte |
| `io/matrix_io.py` | `save_M`/`load_M_checked`/`check_manifest`/`read_manifest`, sidecar JSON (norme de Bloch, unité, `M_normalization`), conversion Ha→eV unique | 35 fichiers | chaîne M, R6 | actif ; docstring et défaut `require_bloch_norm=SUPERCELL` de `load_M_checked` antérieurs à R6 (piège) |
| `io/pseudo_io.py` | `read_upf` (Ry→Ha), `fq_from_fr` ; `read_psp8` ABINIT | non_local, deltav_pw, R4, R5, 5 scripts | chaîne M | actif ; `read_psp8` obsolète |
| `io/wannier_io.py` | `read_w90_mat`, `read_w90_tb` (H(R) eV brut, r(R) Å), `read_w90_hr`, `read_w90_HR` (wrapper), `check_hermicity_HR` | 25 fichiers | chaîne M, EM | actif |
| `io/wannier_provenance.py` | porte de jauge sha256 (`write_wannier_manifest`, `load_wannier_checked`) | 9 fichiers | chaîne T | actif |
| `io/qe_gamma_io.py` | wfc « gamma trick » Γ, `mirror_parity_z`, `density_2d`, `inplane_disc_mask` | R4, R5, R6, R7, R9, R10, `gate_M_normalization` | R4 | figé |
| `io/projwfc_io.py` | lecteurs projwfc.x | `r4_driver` seulement | R4 | figé ; `read_pdos_m` morte |
| `utils/fft_utils.py` | G → indice FFT (`np.rint`) | local_G, local_R, fold_wfk_to_sc, wfk | chaîne M | actif (interne) |
| `utils/lattice.py` | `red_to_cart`, `build_k_path` (R8) ; `monkhorst_pack_grid`, `generate_mp_grid`, `write_kpoints` mortes | non_local, R5, R8, 3 scripts, test_r8 | chaîne M, R8 | actif en partie |
| `utils/planewaves.py` | `mask_invalid_G` ; `make_Cdicts_for_k` morte | non_local, fft_utils, 2 scripts | chaîne M | actif |
| `utils/interpolation.py` | tri-linéaire / spline périodiques | **personne** | — | **orphelin** |
| `wavefunctions/wfk.py` | ψ réel depuis C_nk | `r6_kernel_check`, `ks_reconstruction_all`, `test_ks_reconstruction` | chaîne M | figé |
| `wavefunctions/fold_wfk_to_sc.py` | repli ψ sur la super-cellule, norme maille (v2 depuis `a223687`) | `local_R.compute_ML_R` | chaîne M | figé (interne) |
| `wavefunctions/sc_projection.py` | états Γ de super-cellule sur grille FFT, projection Bloch (promu de R5) | deltav_pw, `gate_M_normalization` | R5→R6 | figé |
| `wannier/wannier_hamiltonian.py` | `Hwr_to_Hwk(Hwr, Rw, k, ndegen)` → (H(k), ε, U) | 16 fichiers (+ `lt.Hwr_to_Hwk`) | chaîne M | actif |
| `wannier/wannier_interpolation.py` | `Mbk_to_Mwk`, `Mwk_to_Mwr`, `Mwr_to_Mwk`, `Mwr_to_Mwk_pairs` (R9), `Mwk_to_Mbk`, `_match_kpoint_order`, `_infer_mp_grid`, `wannier_interpolate`, `ws_images`, `ws_phase` (R10) | 23 fichiers | chaîne M, R9, R10 | actif ; `wannier_interpolate` n'est appelée que par `test_wannier.py` |
| `wannier/supercell_fold.py` | H de super-cellule replié à Γ, LDOS (R9), `sc_planewave_index` | R4, R5, R6, R9, `gate_M_normalization`, test_r9 | R4–R9 | figé ; `kvec_to_rvec` morte |
| `defects/alignment.py` | alignement de potentiel : `vacancy_site`, `far_atom_alignment` (R4), `true_min_image_dist` (R9), `atom_sphere_shifts` (R9, source des C_N) | 10 fichiers | R4, R9, R10 | figé (C_N gelés dans le JSON) |
| `defects/deltav_pw.py` | porte A.2 : ΔV appliqué aux Bloch purs, `RemovedAtomProjector`, `check_pure_bloch` | `gate_M_normalization`, `r6_states_identity` | R5→R6 | figé |
| `defects/local_R.py` | M^L espace réel : `compute_ML_R` (série ; défaut `subtract_mean=True` ≠ production), `prep_realspace_inputs`, `compute_ML_R_mpi` (grossier), `compute_ML_R_mpi_shared` (dense, zero-padding, production), `fourier_resample` | `compute_M`, `compute_M_dense_stages`, `gate_M_normalization`, R5, R6, R9, 2 tests | chaîne M (R6 v2) | actif |
| `defects/local_G.py` | M^L espace réciproque (`compute_ML_G*`, `zero_pad_potential`) | `compute_M_dense_stages --kernel G`, `test_pad_vs_full_supercell`, `test_zero_pad_dense` | chaîne M (historique) | **obsolète** (remplacé par `compute_ML_R_mpi_shared`, `9b58a02`) ; non mis à l'échelle N_cells en R6 mais étiqueté v2 par `--kernel G` (**à vérifier**, chemin non pris par défaut) ; `compute_ML_G_mpi` jamais appelée |
| `defects/non_local.py` | M^NL Kleinman-Bylander (préfacteur 4π/√Ω_uc, norme maille), `build_K_vectors`, `compute_phase`, `compute_angular_part`, `compute_M_NL_mpi` | deltav_pw, local_G, R4, R5, 5 scripts | chaîne M | actif ; `split_counts`, `local_slice` mortes ; `_mpi` n'est appelée que par `_diag_mnl_mpi.py` |
| `defects/many_body/local_tmatrix.py` | cœur de la chaîne T (§3) | 19 fichiers + 3 tests | chaîne T, R9, R10 | actif ; `scattering_rate_from_wannier` (pré-R10, sans C_N) **orpheline** ; argument `perdef` mort |
| `defects/many_body/disorder_average.py` | moyenne sur le désordre Kaasbjerg : `tbar_reduce`, `tbar_k`, `green_k`, `dos_average`, `spectral_path`, `spectral_maxima`, `sigma_eff`, `dirac_*` | `r8_driver`, `test_r8` | R8 | actif (aucun script de `scripts/`) |
| `defects/many_body/pole_criterion.py` | critère de pôle (det, λ), `local_t_cache`, `tbar_pair` | R4, R6, R8, R9, `test_r8` | R4 | figé |
| `defects/many_body/tb_models.py` | banc synthétique graphène π à 5 WF (D6) | R4, `r6_d4`, `test_r8` | R4 | figé (« test material », docstring) |
| `defects/many_body/single_defect.py` | T-matrice de Bloch dense historique (`compute_G0`, `compute_T`, `compute_G`, M en norme super-cellule) | `test_local_tmatrix{,_real}`, `test_local_rcut` (référence du test d'or) ; `compute_spectral`, `compute_tmatrix`, `_eta_scan` (désactivés), `_normtest` | chaîne T pré-2026-09-05 | figé (référence de test) ; **pas orphelin**, contrairement à CLAUDE.md l. 131 |
| `electron_phonon/{phself,selfen}.py` | post-traitement EPW (`linewidth.phself`, elecselfen ; Γ^ep = 2 Im Σ) | `epw_d2_extract`, `epw_phself_post`, `epw_selfen_post`, `make_figures_epw` | EPW | figé ; `RY2EV` morte |
| `electron_photon/{tb_model,kgrid,velocity_operator,ring,kubo,diagnostics}.py` | série EM (k cartésien, eV/Å ; `__all__` de 31 noms complet) | `campagnes/EM/*`, `tests/`, `dev.ipynb` | EM | actif |
| `plotting/plot_psp_radial_proj.py` (+ `these.mplstyle`) | facteurs de forme radiaux | **personne** | — | **orphelin** (probable producteur manuel de `KB_projectors_C.pdf`) |

Fonctions mortes ou orphelines (récapitulatif) : `utils/interpolation.py`, `plotting/`, `local_G.compute_ML_G_mpi`,
`local_tmatrix.scattering_rate_from_wannier`, `non_local.{split_counts,local_slice}`, `pseudo_io.read_psp8`,
`lattice.{monkhorst_pack_grid,generate_mp_grid,write_kpoints}`, `planewaves.make_Cdicts_for_k`, `qe_io.get_typat`,
`projwfc_io.read_pdos_m`, `supercell_fold.kvec_to_rvec`, `selfen.RY2EV`.

---

## 3. Chaîne matrice T (M → M_wk → M_wr → V_loc → g₀ → t → Γ/Σ → DOS/A_k) — pour la réappropriation

Référence canonique = ce qu'exécutent `compute_spectral_wannier.py`, `rcut_resigma.py`, `resonance_metrics.py`,
`resonance_criteria.py` et `r10_driver.py`. Unités : fichiers M en Hartree (sidecar obligatoire), tout le reste en eV
après `load_M_checked(units=EV)` ; H(R) Wannier en eV ; Γ, Σ, η, C_N en eV (×10³ → meV). Index : `M[bra, k', ket, k]`.

Recette minimale (reconstituée de `resonance_metrics.py:25-58` et `compute_spectral_wannier.py:66-109`) :

```python
cfg = load_production(); dp = dense_paths(cfg, S); paths = wannier_provenance.load_wannier_checked(dp["manifest"])
M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=UNIT_CELL, units=EV, require_normalization=M_NORM_V2)  # eV
k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k)
U, kU = read_w90_mat(paths["u"]);       U  = U[_match_kpoint_order(kU, k)]        # (nk, 5, 5)
Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_kpoint_order(kUd, k)]     # (nk, 20, 5)
Hwr, Rw, nd, _, _ = read_w90_tb(paths["tb"])
d = lt.defect_mwr(M, U, Ud, k, MP, n_box=N, C_N=alignment_C(cfg, S)); Mwr, Rn = d["Mwr"], d["Rn"]
lt.mwr_locality(Mwr, Rn)                                                            # garde-fou a
Rloc = Rn[np.linalg.norm(Rn, axis=1) <= cfg["R_cut"] + 1e-9]; V, _ = lt.extract_V_loc(Mwr, Rn, Rloc)
Hk_int, _, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(cfg["nk_int"]), ndegen=nd)
g0 = lt.local_green_batch(Hk_int, k_int, Rloc, egrid, eta)                         # (nE, dim, dim)
t  = lt.local_t(V, g0[j])                                                           # V (1 - g0 V)^-1
Gamma_nk = -2 Im phi_nk^† t(eps_nk) phi_nk                                          # lt.scattering_rate_fast
```

| Étape | Fonction (fichier:ligne) | Formule / convention | Appelée par |
|---|---|---|---|
| 0. Portes | `config.load_production`, `dense_paths`, `alignment_C` (`config.py:10/44/35`) ; `wannier_provenance.load_wannier_checked` (`:44`, sha256 tb/u/u_dis) | toute incohérence lève une erreur | tous les scripts de production |
| 1. M | construction en amont (Ha) : `local_R.compute_ML_R_mpi` (grossier, `compute_M.py`) ou `compute_ML_R_mpi_shared` (dense zero-padding D = pN, `compute_M_dense_stages.py`) + `non_local.compute_M_NL` ; réassemblage v2 `assemble_M2.py` (M2 = N_cells·M^L(v1) + M^NL). Chargement `matrix_io.load_M_checked` (`matrix_io.py:70-106`) | M[n,k',m,k] = ⟨ψ_nk'|ΔV|ψ_mk⟩ sur la super-cellule, ψ normée sur la maille ; refuse sidecar absent, `bloch_norm ≠ unit_cell`, `units ≠ hartree`, `M_normalization` ≠ v2 ; seule conversion Ha→eV | 10 scripts + 8 pilotes |
| 1b. Grille et U | `qe_io.get_k_red` ; `_infer_mp_grid` (`wannier_interpolation.py:208-227`, grille MP complète Γ-centrée) ; `read_w90_mat` + `_match_kpoint_order` (`:186-206`, tol 1e-5) | ligne r de U_dis = bande absolue r (vrai car `dis_win_min = −25 eV` sous la bande la plus basse) | idem |
| 2. M_wk = V†MV | `Mbk_to_Mwk(Mbk, U, U_dis)` (`wannier_interpolation.py:11-54`) | V(k) = U_dis(k)·U(k) ; double boucle (k',k) ; vérifie tr(VV†) = nw ; sortie (5, nk, 5, nk) eV | `defect_mwr`, `wannier_interpolate`, 4 scripts, 7 pilotes |
| 3. M_wr | `Mwk_to_Mwr(Mwk, k_red, MP)` (`:56-89`) ; inverse `Mwr_to_Mwk` (`:91-123`, sans 1/N), rectangulaire `Mwr_to_Mwk_pairs` (`:125-156`, R9) | M_W(R,R') = (1/N_k²) Σ e^{+2πik'·R} M_W(k',k) e^{−2πik·R'} ; boîte R centrée, meshgrid 'ij', nR = D² ; **pas de ndegen, pas de Wigner-Seitz** (étiquettes modulo D) ; commentaire l. 87 annonce une mauvaise forme | `defect_mwr`, 3 scripts, 6 pilotes |
| 3b. Alignement + recentrage (R10) | `local_tmatrix.defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N)` (`local_tmatrix.py:309-344`) → `recenter_mwr` (`:82-97`) ; garde-fou a `mwr_locality` (`:39-54`) | `in_box = all((R mod MP)[:2] < n_box)` sur les étiquettes brutes ; `Mwr[w,in_box,w,in_box] -= C_N` (diagonale R = R' seulement : approximation (i)) ; R_d = argmax ‖M_W[:,R,:,R]‖_F ; Rn = ((R − R_d + D//2) mod D) − D//2 ; `AssertionError` si le max n'est pas en R₀ ; distance = norme des **étiquettes réduites** | `defect_mwr` : 8 scripts + r10 ; `recenter_mwr` seul (chaîne brute pré-R10, sans C_N) : `check_onsite_and_NL`, `test_local_tmatrix_real`, r4–r9 |
| 4. V_loc | troncature **dans les scripts** : `Rloc = Rn[‖Rn‖ ≤ R_cut]` ; `extract_V_loc(Mwr, R, R_local, herm_atol=1e-10)` (`:57-79`) | norme euclidienne des indices réduits (R_cut 3 → 29 mailles, dim 145) ; index plat L·nw + w ; garde-fou b : symétrise ½(V+V†) si résidu > 1e-10 | 6 scripts, 6 pilotes |
| 5. g₀ | `read_w90_tb` ; `Hwr_to_Hwk` (`wannier_hamiltonian.py:7-42`, H(k) = Σ_R e^{+2πik·R} H(R)/ndegen, `eigh`) ; `mp_grid` (`:27`, grille [0,1) non décalée, nk_int = 300 → 90 000 k, **découplée** de la sortie 240²) ; `local_green` (`:100-113`, une énergie par `inv`) ; `local_green_batch(Hwk, k_int, R_local, egrid, eta, deriv=False)` (`:125-159`, restructuré par différences D = R_L − R_L', `_diff_table` `:116`) | g₀[(L,w),(L',w')](ε) = (1/N) Σ_k e^{2πik·R_L} [(ε+iη) − H(k)]⁻¹_ww' e^{−2πik·R_L'} ; `deriv=True` → dg₀/dε ; pas d'énergie η/ne_per_eta = 2,5 meV | `rcut_resigma`, `resonance_*`, `test_local_green_batch`, r4, r8, r9, r10 |
| 6. t | `local_t(V, g0)` (`:162-165`) ; lots : `pole_criterion.local_t_cache` (`:116-121`, threads `R4_EIG_WORKERS`/`SLURM_CPUS_PER_TASK`) ; Born t_B = V + V g₀ V en ligne (`resonance_metrics.py:60`) | t = V (1 − g₀V)⁻¹ par `solve` ; exact si R_local couvre le support de V | `rcut_resigma` ; r4, r6, r8 |
| 7. Σ, Γ sur couche | `scattering_rate` (`:199-230`, exact, lent) ; **`scattering_rate_fast`** (`:264-306`, production : g₀ une fois sur la grille d'énergie, état → point le plus proche, hors fenêtre → NaN, défaut `ne_per_eta=4` ≠ config 8) ; Σ complexe en ligne dans `rcut_resigma.py` | φ_nk[(L,w)] = U_out[k,w,n] e^{2πik·R_L} (`_phase` `:34`) ; Γ_nk = −2 Im φ†t(ε_nk)φ ; `AssertionError` si min Γ < −1e-8 ; Γ « par défaut » intensif = Γ·N_cells ; taux physique = c·Γ | `compute_spectral_wannier`, r9, r10 (`fast`) ; tests (`scattering_rate`) |
| 8. DOS / A_k (R8, Kaasbjerg PRB 101, 045433) | `disorder_average.tbar_reduce` (`:29-67`, τ(D;ε)), `tbar_k` (`:70-104`, T̄_k = Σ_D e^{−2πik·D} τ ; base de bandes U†T̄U), `green_k` (`:107`), `dos_average` (`:127-176`, Σ_k = c·T̄_k, ρ = −(1/πN_k) Σ Im Tr G_k ; `linear=True` = Lloyd), `spectral_path` (`:179`, A_k = −2 Im Tr G_k), `spectral_maxima`, `sigma_eff` (`:233`, Schur 2×2, éq. 44-45), `dirac_*` (`:259-306`, éq. 47-48) | états/eV/maille/spin ; c = `defect_concentration_for_dos` (0,01) ; `resonance_metrics.py` n'utilise que Lloyd (ρ_dis = ρ₀ + c·δρ, δρ = (1/π) Im Tr[t dg₀/dε]) | `r8_driver`, `test_r8` seulement |
| 9. Résonance, LDOS | E_D = milieu π/π* au minimum du gap sur `mp_grid(90)` ; E_res = argmax |Γ| dans ±1,5 eV ; `pole_criterion.det_eig_criterion` (`:43-74`, A = I − V g₀, slogdet, λ_min), `local_minima`, `sign_changes`, `tbar_pair` ; `resonance_criteria.py:~50-58` refait det/λ en ligne ; `local_tmatrix.cluster_ldos` (`:168-196`, G = g₀ + g₀Tg₀, R9) ; `supercell_fold.ldos_from_eigenpairs` (`:156`) | Friedel : ∫δρ (Tr[t g₀′]) contrôlé par Lloyd | `resonance_criteria` ; r4, r6, r8, r9 |

Chemins annexes : `wannier_interpolate` (Mbk→Mwk→Mwr→Mwk fin→`Mwk_to_Mbk`, sans recentrage ; `test_wannier.py`
seulement) ; `ws_images`/`ws_phase` (`wannier_interpolation.py:290-372`, R10 : images de Wigner-Seitz pour les k hors
grille, poids 1/n_tie ; utilisées par `m_rcut_convergence.py` et r10 pour tab:rcut_M, **pas par la matrice t locale**,
d'où la clé `labels_offgrid` documentaire) ; référence dense `single_defect.compute_T` (test d'or
`test_local_tmatrix_real.py`, local = dense × N_cells à 1e-13) ; H replié à Γ `supercell_fold.bloch_folded_hamiltonian`
(diag(ε) + M/N_cells, r4–r9) ; banc synthétique `tb_models` (r4, test_r8).

Scripts de la chaîne, dans l'ordre de production (tous sous `submit_r10.sh` → lanceurs de `scripts/`) :

| Maillon | Script | Sortie |
|---|---|---|
| M (amont) | `compute_M.py` (grossier, MPI), `compute_M_dense_stages.py` (dense), `assemble_M2.py` (v2), porte `gate_M_normalization.py` (A.2) | `results/M2/M_*.npy` + sidecars |
| M → M_wr, localité | `mwr_locality_coarse_vs_dense.py` ; `m_rcut_convergence.py` (M_wr → M fin, R_cut) ; `analyze_M.py` (cartes, échelle, tests, D6) ; `lnl_frobenius_all.py` (L/NL) | `mwr_locality.npz`, `m_rcut_convergence.csv`, `M_analysis.npz`, `M_tests_summary.csv`, `lnl_frobenius.csv` |
| V_loc → g₀ → t → Γ, niveau 1 | `compute_spectral_wannier.py` (`scattering_rate_fast`, carte R_cut × grille × η) | `specwd_<S>_prod.npz` |
| Σ par R_cut, N_k^int | `rcut_resigma.py` (`--rcut`, `--npe`, `--nk-int`) ; `nkint_check_post.py` | `resigma_*.npz`, `nkint_check_9x9.csv`, `m_rcut_resigma.csv` (via r6/r10) |
| Γ(ε), Born, δρ, ρ_dis, T̄(K), C14 | `resonance_metrics.py` | `resonance_<S>.npz`, `_shiftL` |
| det/λ, Friedel/Lloyd, Γ à c = 0,1 % | `resonance_criteria.py` | `resonance_criteria_<S>.npz` |
| niveau 2, familles | `level2_families.py` ; `make_figures.py` (csv) | `level2_families.csv`, `level{1,2}_summary.csv` |
| pont ch. 5 | `epw_ed_vs_ep.py` | `ed_vs_ep_*.npz` |
| figures | `make_figures.py`, `make_figures_memoire.py`, `make_figures_epw.py` (`fig_epw_vs_ed`) | `figures/` |
| validation | `test_local_tmatrix.py`, `test_local_rcut.py` (synthétiques), `test_local_green_batch.py` (H 27×27), `test_local_tmatrix_real.py` (test d'or réel, **cassé depuis R10**) | PASS/FAIL |
| DOS / A_k | `campagnes/R/R8_kaasbjerg/r8_driver.py` seulement | `campagnes/R/R8_kaasbjerg/out/` |
| ancienne chaîne de Bloch (obsolète) | `compute_spectral.py`, `compute_tmatrix.py`, `compute_convergence.py`, `_eta_scan.py`, `_normtest.py` | — |

Paramètres réels et contrôles C1–C18 : NOTES_TGAMMA.md §1–§3, **à l'état final R10 depuis `90d2022`** (valeur finale en
tête, « [v1 : …] » entre crochets ; valeur non alignée seulement au §8) ; table complète v1 → non aligné → final :
`campagnes/M/ch4/table_v1_final.md`. Nomenclature fixée le 2026-09-30 (à utiliser dans le mémoire ; les archives gardent
leurs étiquettes) : **non aligné** (= « tel quel », `results/M2`), **alignement à site unique** (= « Lu »), **alignement
de Kumagai–Oba** (= « plateau (i) », ≥ 0,75 r_max), **final** = v2 + Kumagai–Oba (`results/M2_plateau`), symbole
ΔV_PA^(N) = C_N. Points à retenir en se réappropriant : R_cut est une norme sur les indices réduits, pas une distance
cartésienne ; deux « médianes de Γ » coexistent (états vs courbe) ; E_res est discret et **retiré du ch. 4** (R10 B.1) ;
N_k^int 300 → 600 change Γ_T(E_D) de +4,8 % ; la matrice t ne voit pas les images de Wigner-Seitz ; C14 est redéfini
(±rms du plateau, ±9,05 meV à 9×9) ; la production n'a jamais soustrait la moyenne de ΔV (`subtract_mean=False`).

---

## 4. `scripts/` (81 fichiers suivis : 67 `.py`, 14 `.sh`)

`scripts/run.py` et `scripts/compute_M_cluster.py`, cités par CLAUDE.md (l. 59, 63-65, 217) et `link_data.sh:2`,
**n'existent pas** (supprimés par `cf05e92`, 2026-06-19 ; remplacés par `compute_M.py` et `compute_M_dense_stages.py`).
Colonnes : rôle · campagne · qui l'utilise · statut · mémoire.

### 4.1 Production de M (chaîne M)

| Script | Rôle | Campagne | Utilisé par | Statut | Mémoire |
|---|---|---|---|---|---|
| `compute_M.py` | pilote M grossier en 3 étapes `ml` (MPI, `compute_ML_R_mpi`) / `nl` / `combine` (v2) | chaîne M, R6 | `submit_M.sh`, `submit_r6_kernel.sh`, `_test_mnl_mpi.sh`, r9 | figé (outil valide) ; l. 1-2 avant le shebang, `RES` inutilisé, commentaires l. 124/141 vs 150-152 contradictoires | ch. 4 |
| `compute_M_dense_stages.py` | M dense (pN)² par zero-padding, `ml`/`nl`/`combine`, noyau R (défaut) ou G, `--coarse` | chaîne M dense, R4, R6 | `submit_M_dense.sh`, `submit_r6_kernel.sh`, r4, r9 (`paths`) | figé ; `PFAC` et scratch en dur (l. 23-24) dupliquent la config ; `--kernel G` étiquette v2 un M probablement v1 | ch. 4 |
| `assemble_M2.py` | M2 = N_cells·M^L(v1) + M^NL par blocs memmap, sidecars v2, refus d'écrasement | R6 | `etape{1,2}/submit_r6.sh` ; `results/M2/README.md` | figé (outil officiel de reconstruction des 43 M2, non miroités) | ch. 4 |
| `gate_M_normalization.py` | porte A.2 : ΔV appliqué aux Bloch purs vs M/N_cells (1e-6 eV, code 3 = refus) | R6 | `submit_r6.sh`, `submit_r6_kernel.sh`, R10 c1post | actif ; `--m2-dir results/M2` en dur cohérent avec `matrices_dir` | ch. 4 (tab:tests_M) |
| `tag_vacancy_sublattice.py` | écrit `vacancy_sublattice` A/B dans les sidecars | chaîne M | manuel | **cassé silencieusement** (glob l. 12 sur `results_dir` : « 0 sidecars tagged ») | métadonnée |
| `finalize_wannier.py` | rapport d'étalements, copie tb/u/u_dis dans `wannier/{N}`, écrit le manifeste | chaîne M/T, EM1 | manuel (`submit_dense_nb20.sh` rorqual) | figé (défauts l. 85-90 périmés : tailles 5/7/8, sans `--dense`) | ch. 4 (entrées), §2.5 |
| `link_data.sh` | construit `data/graphene/**` en liens vers le scratch et `Vks_*` | infra | manuel (README) | actif sur rorqual ; cite `run.py` (l. 2) ; nommage `defect_<N>.save` ≠ local `defect_unit_cell_<N>.save` | infra |
| `migrate_M_norm.py` | `M_ed_*` → `*_norm.npy` (÷N_cells) | pré-2026-09-05 | — | **obsolète** (`SystemExit` l. 2) | — |
| `_diag_mnl_mpi.py`, `_test_mnl_mpi.sh` | diagnostic de `compute_M_NL_mpi` (corruption de tas nk ≥ 81) | chaîne M | — | figé (clos) ; `.sh` écrit dans `$RES/_test_mnl` (M2_plateau) | — |

### 4.2 Contrôles de M

| Script | Rôle | Campagne | Utilisé par | Statut | Mémoire |
|---|---|---|---|---|---|
| `analyze_M.py` | carte |M| à K, V_ed ligne, Re M^L/M^NL à K, échelle max|M|(N), table des tests, variante alignée D6 avec portes | chaîne M, R6, R10 | `submit_post.sh analyze`, `submit_r10.sh` ; lu par `make_figures*`, `level2_families` | actif ; l. 231 chiffre codé en dur (job 20238555) | ch. 4 (§4.1.5, tab:tests_M) |
| `analyze_Ved.py` | V_ed^L : carte, profils, moyenne azimutale, anomalie 1,42 Å | chaîne M, R10 | `submit_r10.sh` (l. 132) | actif | ch. 4 (fig_Ved) |
| `ks_reconstruction_all.py` | reconstruction KS H = T + ⟨V_p⟩ + V^NL sur grossiers et denses | chaîne M, R6 | `submit_post.sh ksrec`, `submit_r6_kernel.sh` | figé (R10 a copié le npz) ; `sys.path.insert(0,"scripts")` relatif exige cwd = racine | ch. 4 (fig_ks_reconstruction) |
| `sampling_table.py` | table N, N mod 3, ΔE_F, V_ed frontière/radial | chaîne M, R10 | `submit_post.sh ksrec`, `submit_r10.sh p_ved` | actif (lancé en local, écraserait le CSV avec des NaN : seule la 5×5 est présente) | ch. 4 (tab:échantillonnage) |
| `lnl_frobenius_all.py` | ⟨‖M^NL‖⟩/⟨‖M^L‖⟩ blocs π/π*, brut et aligné D6 | chaîne M, R10 | `submit_post.sh analyze`, `submit_lnl_frobenius.sh`, r10 | actif | ch. 4 (tab:L_NL) |
| `_bz_ratio_LNL.py` | même rapport L/NL (ponctuel R6 J8) | R6 | manuel | figé (doublon de `lnl_frobenius_all`) | ch. 4 (chiffre) |
| `check_ML_coarse_kernel.py` | noyau partagé sur la grille grossière = `M_L_<S>` (1e-10) | chaîne M dense, R6 | — (cité R6_rapport) | **cassé depuis R10** (l. 5, 8-9 : `results_dir`) ; doublon d'`analyze_M` l. 197-205 | tab:tests_M (via analyze_M) |
| `check_M_dense_vs_coarse.py` | k communs dense/grossier : valeurs singulières, Frobenius | chaîne M dense, R6 | `submit_golden_dense.sh` | **cassé depuis R10** (l. 8, 15-17) ; doublon d'`analyze_M` l. 173-183 | idem |
| `check_M_dense_vs_coarse_bands.py` | idem résolu en bandes (nb_sub 4-16) | chaîne M dense | personne | obsolète (orphelin, sans porte v2, code mort l. 20-21) | — |
| `check_M_dense_nb20_vs_nb16.py` | non-régression nbnd 16 → 20 | chaîne M dense | — | **obsolète / non rejouable** (nb16 sur nearline, `SystemExit` l. 6) + chemin R10 | tab:tests_M (« non rejouable ») |
| `check_onsite_and_NL.py` | sur-site p_z–p_z de V_loc, M^L/M^NL à K et K' | chaîne M dense, R4 | — (rapports R4, R6, R10) | **cassé depuis R10** (l. 9) et obsolète (chaîne brute sans C_N) ; remplacé par `mwr_locality_coarse_vs_dense` et `analyze_M` §3 | — |
| `validate_ML_grid_7x7.py` | M^L de référence 7×7 (noyau série) vs MPI et ancienne grille 216 | chaîne M (bogue de grille) | — | **cassé / non rejouable** (argv sans défaut l. 9, garde l. 15, chemin l. 16) | — |
| `_old_vs_new_7x7.py` | M dense 7×7 ancienne vs nouvelle grille | chaîne M | — | obsolète (`obsolete_grid_7x7` supprimé, garde l. 6) | — |
| `_mcheck.py` | normes max|M|, ‖M‖_F de `M_ed_{5..8}` | chaîne M | — | **cassé depuis R10** (l. 4, 7, 11) | — |
| `_normtest.py` | médiane Γ avec M brut vs M/nk (`compute_T`) | pré-2026-09-05 | — | **cassé** (chemin l. 11 ; mélange Ha/eV l. 11-12) → à classer obsolète | — |

### 4.3 Chaîne T (production, ch. 4)

| Script | Rôle | Campagne | Utilisé par | Statut | Mémoire |
|---|---|---|---|---|---|
| `compute_spectral_wannier.py` | niveau 1 : `defect_mwr` → V_loc → `scattering_rate_fast`, carte (R_cut, grille, η), porte de jauge | chaîne T (gel 2026-09-05), R6, R10 | `submit_spectral_wannier{,_dense}.sh`, `submit_r10.sh` | actif | ch. 4 (fig_convergence, fig_rcut, fig_plateau, fig_level2) |
| `resonance_metrics.py` | Γ_T, Γ_Born, Γ/ρ₀, δρ, ρ_dis, T̄_ππ(K), `--shift-L-meV` (C14) | chaîne T, R6, R10 (D10) | `submit_post.sh resonance|c14`, `submit_r10.sh` | actif (`np.trapz` déprécié l. 62) | ch. 4 (fig_spectral*), ch. 5 (fig_epw_vs_ed) |
| `resonance_criteria.py` | det/λ par bloc (full, π, σ), Friedel/Lloyd, Γ à c = 0,1 % | chaîne T, R6, R10 | `submit_post.sh resonance|criteria`, r10 | actif ; exige `resonance_<S>.npz` avant | ch. 4 (fig_spectral*) |
| `rcut_resigma.py` | Σ_nk complexe par R_cut, `--npe`, `--nk-int` | chaîne T, P13, R6, R10 | `submit_rcut_resigma.sh`, `submit_nkint_check.sh`, r9, r10 | actif | ch. 4 (C9-C11, tab:rcut_M via csv) |
| `nkint_check_post.py` | tableaux N_k^int (médianes, E_res, Γ_T(E_D)) | P13, R6, R10 | `submit_post.sh figures`, `submit_r10.sh c1f` | actif ; l. 1-2 avant le shebang | ch. 4 (NOTES_TGAMMA §6) |
| `m_rcut_convergence.py` | convergence de R_cut au niveau de M (images WS) | chaîne T, R6, R10 (b) | `submit_post.sh locality`, r9, r10 | actif ; csv **en ajout** (l. 62-66) | ch. 4 (tab:rcut_M) |
| `mwr_locality_coarse_vs_dense.py` | sur-site p_z–p_z et décroissance ‖M_W(R,R₀)‖, grossier vs dense | chaîne T, R6, R10 (b) | `submit_post.sh locality`, r10 | actif | ch. 4 (fig_locality*) |
| `level2_families.py` | table niveau 2 par famille 3m / non-3m | niveau 2, R10 (D6) | `submit_post.sh figures`, `submit_r10.sh`, r10 | actif | ch. 4 (tableau des familles) |
| `summarize_level1_maps.py` | synthèse niveau 1 depuis les **logs** `specwd_*.out` | chaîne T | personne | obsolète (remplacé par `level1_summary.csv` ; lit les logs, contraire à la règle) | — |
| `compute_spectral.py`, `compute_tmatrix.py`, `compute_convergence.py` | ancienne T-matrice de Bloch dense (`compute_T`, `compute_G0/G`) et sa convergence ; lisaient `M_ed_*_norm.npy` (supprimés) | pré-2026-09-05 | `submit_spectral.sh`, `submit_tmatrix.sh` | **obsolètes** (`SystemExit` l. 2 / garde l. 36) | — |
| `_eta_scan.py` | sensibilité à η (`compute_T`) | pré-2026-09-05 | — | obsolète (`SystemExit` l. 7) | — |
| `epw_ed_vs_ep.py` | Γ^ed(c)/Γ^ep(T) point par point | EPW/ch. 5, R6, R10 | `submit_post.sh figures`, `submit_epw_p2_post_mv.sh`, `submit_r10.sh` | actif ; l. 1-2 avant le shebang | ch. 5 |

### 4.4 Figures et tables

| Script | Rôle | Figures produites (lues dans `results_dir` = M2_plateau) | Statut | Mémoire |
|---|---|---|---|---|
| `make_figures.py` | figures de travail ch. 4 (14) + `level{1,2}_summary.csv` | fig_rcut, fig_plateau, fig_level2 (specwd) ; fig_locality (mwr_locality) ; fig_spectral (resonance_9x9 + criteria) ; fig_M_map, fig_Ved_boundary, fig_M_scaling (M_analysis) ; fig_ks_reconstruction ; fig_Ved_map, fig_Ved_profile_mean, fig_Ved_radial, fig_Ved_zoom, fig_Ved_radial_masked (ved_analysis) | actif ; code mort l. 25 ; style par chemin relatif (cwd = racine) | ch. 4 (fig_level2, fig_convergence via memoire, fig_ks_reconstruction, fig_Ved_zoom inclus) |
| `make_figures_memoire.py` | versions finales ch. 4 (6) | fig_convergence (specwd ×6), fig_locality_final, fig_spectral_final, fig_M_map_final, fig_M_scaling_final, fig_Ved (ved_analysis + C_N) | actif ; aucune garde d'existence ; docstring 6,5×7,2 vs code 6,5×6,0 | ch. 4 (5 incluses + fig_Ved) |
| `make_figures_epw.py` | figures ch. 5 (8) | fig_epw_validation, fig_epw_gamma, fig_epw_vs_ed (+ resonance_9x9), fig_epw_g_control, fig_epw_phonselfen, fig_epw_decay, fig_epw_kohn_degauss, fig_epw_phonons (`results/epw/*`) | actif ; défauts `--prod-tag prod` et `phdos_24k24q` inexistants : sans les options de `submit_post.sh` l. 44-45, deux figures manquent et l'ancienne chaîne mv0.002 est tracée ; docstring « trois figures » | ch. 5 (7 incluses) + ch. 2 (fig_epw_phonons) |
| `make_figures_em.py` | figure §2.5 | fig_em_coupling (`campagnes/EM/{M4_sigma,EM3,EM2}`) | actif ; `ROOT` absolu (robuste) | §2.5 (pas encore incluse) |
| `_palette.py` | palette fixe (NAVY #000080, …, `CMAP_SEQ/DIV`) | — | actif ; importé par 4 `make_figures*`, r4/r5/r8/r9/r10, EM2/EM3/B | toutes |
| `_bands.py` | helpers de chemins de bandes | — | actif (`compare_bands_*`, r9) ; l. 16 : image de K (2/3, −1/3) n'est pas un point K (étiquettes seulement) | — |
| `_paths.py` | chemins locaux `data/graphene` (`EDI_DATA`, nommage `defect_unit_cell_<N>.save`) | — | actif (6 scripts de validation) | — |

### 4.5 EPW (ch. 5)

| Script | Rôle | Campagne | Utilisé par | Statut |
|---|---|---|---|---|
| `epw_pp_save.py` | équivalent pp.py : dvscf/dyn/phsave → `save/` | EPW P0/P1 | manuel | figé ; `_ph0` supprimé, non rejouable ; `save/` irremplaçable |
| `epw_extract_gkk.py` | découpe prtgkk d'`epw.out` | P1 | `submit_epw_p1_post.sh` | figé |
| `epw_validate.py` | bandes, phonons, décroissance, |g| DFPT vs EPW (fsthick) → `validation_<tag>.npz` | P1, P18 | `submit_epw_p1_post.sh` ; lu par `make_figures_epw`, `epw_d2_extract` | figé ; maille en dur l. 29 |
| `epw_selfen_post.py` | elecselfen → Γ^ep = 2 Im Σ, médianes → `selfen_<tag>_T<T>.npz` | P2/P4 | `submit_epw_p2_post_mv.sh` | figé |
| `epw_phself_post.py` | phonselfen → ω, γ, λ → `phself_<tag>.npz` | P6/P7 | idem | figé ; `E_D = −4.2389` en dur l. 34 |
| `epw_d2_extract.py` | ⟨D²_Γ⟩, ⟨D²_K⟩ par trois routes → `d2_extract_<tag>.npz` | P7 | idem | figé ; chemin cluster absolu l. 28 ; maille en dur l. 35 |
| `epw_ring_check.py` | convention gamma de EPW sur anneaux → `ring_check_<tag>.npz` | P11 | idem | figé ; chemin absolu l. 14 |
| `epw_dfpt_path_freq.py` | fréquences ph.x sur Γ–K–M–Γ → `dfpt_path_freq_<tag>.npz` | P18 | `make_figures_epw` | figé |
| `epw_phdos_extract.py` | DOS de phonons matdyn → `phdos_<tag>.npz` | P18 | `make_figures_epw` | figé |

### 4.6 Validation Wannier et tests autonomes (PASS/FAIL)

| Script | Rôle | Données | Exécutable en local | Statut |
|---|---|---|---|---|
| `test_ks_reconstruction.py` | A : H = T + V_loc + V^NL = diag(ε) ; B : défaut nul ⇒ M = 0 ; C : V_p restreint | `_paths` 5×5 (présent) | oui | actif ; fournit `sc_pot_on_uc_grid` et `reconstruct_ks_hamiltonian` à `ks_reconstruction_all`, `run_test_A_batch`, `r6_kernel_check` |
| `test_wannier.py` | 5 tests d'interpolation Wannier (parseurs, aller-retour FT, spectre, grille fine, M réel) | 11×11 locale ; `hr.dat` et `M_ed.npy` absents → sous-tests ignorés | oui | actif ; docstring « 5×5 » périmée |
| `test_zero_pad_dense.py` | exactitude de `zero_pad_potential` ; non-régression dense vs `compute_ML_G` | test 1 sans données ; test 2 5×5 | oui | actif |
| `test_pad_vs_full_supercell.py` | M^L zero-padded ×p vs super-cellule p fois plus grande | paires 6×6 et 12×12 | non | figé |
| `test_local_tmatrix.py` | test d'or synthétique (t local sur site = `compute_T` dense × N_k) | aucune | oui | actif |
| `test_local_rcut.py` | R_cut synthétique, `extract_V_loc`, `mwr_locality` ; importe `random_H` | aucune (lancer `python scripts/…`) | oui | actif |
| `test_local_green_batch.py` | `local_green_batch` = `local_green` ; `scattering_rate_fast` = `scattering_rate` (C8) | `wannier/27x27` (versionné), cwd = racine | oui | actif ; **pas de code de sortie** |
| `test_local_tmatrix_real.py` | test d'or réel bloquant (C6) sur M grossier ou dense | scratch + M | non | **cassé depuis R10** (l. 18, 27, 29 : `results_dir`) ; **aucun `SystemExit` avec code** : FAIL sort en 0 (l. 63) ; scratch en dur l. 26, table `PF` dupliquée l. 23 |
| `validate_wannier_bands.py` | Wannier vs DFT sur la grille MP grossière (avec/sans ndegen), 2 PNG dans le cwd | 11×11 locale | oui | actif ; docstring périmée ; `get_fermi` sous `hasattr` |
| `compare_bands_qe.py` | H(R) interpolé vs `bands.dat` sur un chemin → `results/compare_bands_qe.png` | 11×11 locale | oui | figé |
| `compare_bands_w90_qe.py` | `wannier_band.dat` vs `bands.dat` → `results/compare_bands_w90_qe.png` | 11×11 locale | oui | figé |
| `run_test_A_batch.py` | test A sur les 8 mailles → `results/test_A/` | nommage cluster | non | figé (répertoire `test_A` inexistant ; `results/test_recon` orphelin) |

### 4.7 Lanceurs SLURM (`submit_*.sh`, 14)

Traits communs : `PROJ=${GRAPHENE_RAMAN:-git rev-parse}`, `PYTHONPATH=$PROJ/src`, `RES=results_dir(cfg)`
(= M2_plateau depuis R10), `#SBATCH --output` sur le littéral `results/M2/logs/` (12 lanceurs). Aucun chemin
`ab-initio-defects` ni `results/M/` dans les options SBATCH.

| Lanceur | Lance | Statut |
|---|---|---|
| `submit_M.sh` | `compute_M.py ml → nl → combine` (tableau 5…12) → `$RES/M_*` | figé ; **écrirait les M dans M2_plateau** (l. 40-42) au lieu de `matrices_dir` |
| `submit_M_dense.sh` | `compute_M_dense_stages.py` → `$RES/M_*dense_*` | figé ; même défaut (l. 25-33) |
| `submit_golden_dense.sh` | `test_local_tmatrix_real.py --dense` puis `check_M_dense_vs_coarse.py $RES/M_ed… $RES/M_dense…` | **cassé depuis R10** (l. 25) ; `rc=$?` vaut 0 sur FAIL |
| `submit_spectral_wannier_dense.sh` | `compute_spectral_wannier.py --dense` → `specwd_<S>_prod.npz` | actif (R10 C1 ×6) |
| `submit_spectral_wannier.sh` | version grossière → `specw_<S>.npz` | obsolète (aucun `specw_*` conservé) |
| `submit_rcut_resigma.sh`, `submit_nkint_check.sh` | `rcut_resigma.py` par R_cut / par N_k^int | actifs (R10 C1) |
| `submit_lnl_frobenius.sh` | `lnl_frobenius_all.py` | actif |
| `submit_post.sh` | tâches `locality`, `analyze`, `ksrec`, `resonance`, `criteria`, `c14`, `figures` | actif en partie ; périmé : `ksrec` copie `results/M/ved_analysis.npz` (v1, l. 31) ; `c14` `--shift-L-meV 25,-25` (R6 ; R10 = ±9,05, l. 36) → relancer écraserait `resonance_9x9_shiftL.npz` ; `figures` sans `--outdir` |
| `submit_test_A.sh` | `run_test_A_batch.py` | figé ; pas de `PYTHONPATH` |
| `submit_spectral.sh`, `submit_tmatrix.sh` | `compute_spectral.py`, `compute_tmatrix.py` | **obsolètes** (scripts désactivés) |
| `submit_epw_p1_post.sh` | `epw_extract_gkk.py`, `epw_validate.py` | figé ; chemin EPW absolu l. 14 ; pas de `PYTHONPATH` |
| `submit_epw_p2_post_mv.sh` | `epw_selfen_post` ×3, `epw_phself_post` ×2, `epw_d2_extract`, `epw_ring_check`, `epw_ed_vs_ep` | figé ; chemin absolu l. 16 ; **pas de `PYTHONPATH`** alors que 4 scripts importent le paquet (risque `ImportError` sur rorqual sans install éditable) |
| `_test_mnl_mpi.sh` | diagnostic M^NL MPI vs série | figé (clos) |

---

## 5. `article/` (720 fichiers suivis) — campagnes R

Toutes les copies `article/` sont des copies versionnées de répertoires de travail rorqual : les pilotes calculent
`PROJ` et `R4DIR`/`R5DIR` à partir de l'arborescence `graphene/qe/defects/…` et **ne s'exécutent pas sur place en
local** (sauf `GRAPHENE_RAMAN` + arborescence rorqual). Aucun import du paquet cassé. `ase` et `spglib` (R1, R2) ne
sont pas dans le venv local.

| Campagne | Répertoire | But | Statut déclaré | Pilote(s) et imports notables | Figures / tables → pour | Statut réel |
|---|---|---|---|---|---|---|
| R1 (+R1b) | `R1_vacancy_relaxed` | relaxation 9×9 lacune nspin 1/2 ; scf 3×3×1 de contrôle | non déclaré (prépa article), miroité, 2026-09-22 | `make_inputs.py`, `analyze_relax.py` (ase, spglib), `k3x3/*` ; aucun import du paquet | tables des rapports → article | figé (géométrie relaxée ; sert à R4 J0, plan C5) |
| R2 | `R2_size_series` | relaxations 5×5…12×12 | PRODUCTION 2026-09-24 | `make_inputs_series.py` (lit sidecars `results/M/*.json` rorqual l. 15), `analyze_relax_series.py` | tables → article | figé (source du transplant R7c) |
| R4 | `R4_quasi_lie` | états π/σ de la lacune 9×9 en DFT vs chaîne M → Wannier → T | TEST 2026-09-25 | `r4_driver.py` (prep, d6, d1, d5, r, d4, d3, g0, tables, figs) : qe_io, matrix_io, wannier_provenance, wannier_io, qe_gamma_io, projwfc_io, pseudo_io, wannier_interpolation, supercell_fold, local_tmatrix, pole_criterion, tb_models, alignment, non_local ; `_palette` | `fig/d1_states, d3_alpha_*, d4_ladder, d5_boundary, d6_*` → article ; `sections/tables.md` | **obsolète pour M** (v1 ; erratum R5 ; D3/D4 refaits par R6 étape 2) ; rejeu incohérent : `RES = results/M` (l. 42) mais `dense_paths` → M2 (l. 109, 216) |
| R5 | `R5_base_vs_M` | écart H_p+M vs QE : base ou M ? ; nscf 128 bandes | TEST 2026-09-25 | `r5_driver.py` (**importe `r4_driver`**, 29 noms), `r5_deltav_pw.py` (promu en `defects/deltav_pw.py`), `r5_sc_projection.py` (promu en `wavefunctions/sc_projection.py`), `r5_basis_diagnostics.py`, `r5_alignment_ext.py` | `fig/a_delta, a_ladder_corrected, size, delta_vs_n, ladder_bands` → article ; tables A/B/C | figé ; A.2 (facteur 81) à l'origine de R6 ; modules `r5_*` restent importés par R6, R7, R9 ; rejeu l. 171 mélange v1/v2 |
| R6 | `R6_production_corrigee` | M2 = N_cells·M^L + M^NL, porte A.2, D4, régénération ch. 4 | PRODUCTION 2026-09-25→27 | `phase0/*.diff` (appliqués) ; `etape1/r6_kernel_check.py` (importe `scripts/test_ks_reconstruction`), `r6_states_identity.py` (importe `r5_driver`) ; `etape2/r6_d4.py` (importe `r5_driver`, lit `results/M2` l. 40) ; `etape3/r6_compare_v1_v2.py`, `r6_level1_gate.py`, `r6_m_rcut_resigma.py`, `r6_tests_gate_row.py`, `refactor_3_0.py` (patcheur appliqué), `runbook_3.sh`, `submit_r6_final.sh`, `install_figures_R6.sh` | `etape2/fig` (validation) ; `etape3/figures_v2` (28 fig + `_mv0.002`) et `figures_fix` → mémoire ch. 4/5 **remplacées par R10** ; `table_v1_v2.md`, `level1_gate.md`, `gate_table.md` ; `csv_v2` = `results/M2/*.csv` | figé ; **`etape3/r6_compare_v1_v2.py:14` TypeError latent** (`results_dir_frozen` est une liste) ; les 3 autres scripts etape3 visent M2_plateau (dérive de sens) |
| R7 (+R7c) | `R7_tailles_3m` | scf 15×15…27×27 famille 3m ; R7c relaxations nspin1 | PRODUCTION 2026-09-25/26 | `r7_driver.py` (importe `r4_driver`, `r5_driver`, `r5_alignment_ext`), `make_inputs_r7{,_relax}.py`, `phase0_inventory.py` | `fig/size_3m, localized_3m` (+ `_d1_8pts`, `_d1_reg`) → article ; `d1*/D1_tables.md` | figé (D1 à 8 points = référence R9 C, R10 A.2) ; **R7c inachevé dans la copie** : `relax.out` pour 15 et 18 seulement |
| R8 | `R8_kaasbjerg` | Fig. 13/14 de Kaasbjerg (DOS, A_k) avec M2 et t local | TEST 2026-09-27→29 | `r8_driver.py` (extract, prep, g0, gate, dos, spec, sens, fig, sigeff, dirac) : `disorder_average.*`, `pole_criterion.{local_t_cache, sign_changes}`, `lattice.build_k_path`, chaîne T ; lit `results/M2/{MD5SUMS, resonance_9x9}` explicitement (l. 575, 897), C_N de `R9/a/a1_results.json` ; sorties dans `out/` | `fig/dos_c` (étape 2), `spectral_GKM` (3), `sigma_K` (4), `sensibilites` (5), `superposition` (7b) → **mémoire ch. 4** (portée élargie à l'étape 4 le 2026-09-30, `241739e`) ; hors mémoire : `7a_controle`, `*_rho0_300` ; étape 6 (Dirac, tableau) → article | figé (base pré-R10, sans WS ; même C_N plateau) ; seul utilisateur de `disorder_average` ; ses 5 figures ne sont ni dans `figures/` ni dans le dépôt du mémoire |
| R9 | `R9_controles` | A : C_N ; B : E_res vs N_k^int ; C : chaîne repliée vs R7 ; D : Kaasbjerg ; clôture plateau ; audit image minimale | TEST, clos 2026-09-29 | `r9_driver.py` (17 sous-commandes) : chaîne T complète, `alignment.atom_sphere_shifts`, `local_R.*`, `supercell_fold`, `cluster_ldos` ; importe `scripts/compute_M_dense_stages.paths`, `_bands`, `_palette`, `r5_driver` ; `cloture/r9_driver_propose.py` + `submit_r9_propose.sh` (proposition remplacée, 525 lignes de diff) | `fig/offset_profiles, resonance_vs_nkint{,_R9,_plateau}, rcut_aligned, folded_vs_R7, kaasbjerg_fig3_map` → contrôles ch. 4 (non installées) ; tables A1/A3/B/C/D, `cloture/synthese.md`, `audit_image_minimale.md` | figé ; plateau promu par R10 ; **rejeu cassé** : `r9_driver.py:855, 903` lisent `M_ed_<S>.npy` dans `results_dir` (= M2_plateau) ; l. 607/724/973/995 changent de sens ; `cloture/*_propose*` obsolètes |
| R10 | `R10_plateau` | base unique du ch. 4 : A C_N plateau (13 tailles), B E_res/familles, C rejeu de la production (`defect_mwr`, C_N config, images WS) | TEST (sorties C = production ch. 4), **CLOS 2026-09-30** | `r10_driver.py` (a0…c6) : `config.{alignment_C, results_dir, dense_paths}`, chaîne T (`defect_mwr`, `scattering_rate_fast`, `ws_images`, `ws_phase`), `alignment.*` ; aucun import de r4–r9 (lit leurs JSON) ; `submit_r10.sh` (c1 = lanceurs de `scripts/`, p_ved, p_c14, c1f) ; lit `results/M2` (`M2_DIR`), écrit `results/M2_plateau` | `fig/fig_*` (28 PNG **identiques au bit à `figures/`**) → mémoire ch. 4 + ch. 5 ; contrôles `offset_profiles_13, levels_vs_invN, eres_vs_grid, kaasbjerg_plateau_ws, M_map_brut_aligne, dV_z_profiles`, `fig/avant_apres/` ; tables A1-3, B, C0/C2/C3/C6, `c/table_v2_plateau.md`, `manifeste_R10.md` | **actif** (définit la base vivante) ; README l. 6 et rapport l. 628/658 citent le hash `23ee3bb` **réécrit** (devenu `c2bc733`) ; README se termine sur « - Lecture : » vide (comme R9) |
| C | `C_optique_lacunes/PLAN.md` | σ(ω) Kubo-Greenwood avec Σ = c·T̄ (puis SCTMA, spin) ; réutilise la chaîne `r8_driver` | plan, 2026-09-29 (perspective C d'EM.md) | — | → article (futur) | plan ; désigne `results/M2/resonance_*` (périmé) ; dernière ligne « ImportError » périmée |

Dépendances entre pilotes (à garder à l'esprit avant tout déplacement) : `r5_driver` ← `r4_driver` ;
`r6_states_identity`, `r6_d4`, `r7_driver`, `r9_driver` ← `r5_driver` (et `r5_sc_projection`, `r5_alignment_ext`) ;
`r6_kernel_check`, `ks_reconstruction_all`, `run_test_A_batch` ← `scripts/test_ks_reconstruction` ; `r9_driver` ←
`scripts/compute_M_dense_stages`, `scripts/_bands` ; r4/r5/r8/r9/r10, EM2/EM3/B ← `scripts/_palette`.

---

## 6. `memoire/` (168 fichiers suivis) — série EM (§2.5) et chiffres du ch. 4

### 6.1 `campagnes/EM/` (série EM, §2.5)

`campagnes/EM/EM.md` : plan M0–M4, fonctions F1–F19, décisions verrouillées, statut. Toutes les fonctions existent dans
`electron_photon/` sauf **F16 `shift_home_cell`** (reportée, absente de src). M0–M4 faits (2026-09-28/29) ;
perspectives A (kT, dopage) et B (σ complexe) faites, C → `article/C_optique_lacunes`, D ouverte ; P28 (texte du
mémoire) à faire. `EM2_prompt.md` (prompt EM2), `notes_numpy_pytest.md` (pièges NumPy/pytest de M0).

| Campagne | Répertoire | But | Statut déclaré | Scripts (imports `electron_photon`) | Produit → pour | Statut réel |
|---|---|---|---|---|---|---|
| EM1 | `EM1_tb` | r(R) et `_tb.dat` de la wannierisation 27×27 (`restart = plot`) | PRODUCTION 2026-09-25 | `em1_check.py` (aucun import du paquet) | données §2.5 (`wannier/27x27/wannier_tb.dat`) | figé ; README l. 14 et rapport l. 5 à corriger sur rorqual puis `cp -p` (mémoire persistante) |
| EM2 | `EM2` | DFT directe (`bands.x lp`) et postw90 `kubo` vs M4 | PRODUCTION 2026-09-29 | `em2_kpoints.py`, `em2_A_compare.py`, `em2_A_figure.py`, `em2_B_ours.py` (`kgrid`, `kubo`), `em2_B_compare.py` (aucun), `em2_B_decompose.py`, `em2_B_figure.py` (importent `em2_B_compare`) ; `REPO` exige `PROJECTS` ou `GRAPHENE_RAMAN` | `em2_A_anneaux`, `em2_B_sigma` (validation §2.5) ; `em2_A.npz`, `em2_postw90_sigma.npz` | figé ; `__pycache__` local ignoré |
| EM3 | `EM3` | données de la figure et tableau du §2.5 (anneau 2,33 eV, `transl_inv`) | PRODUCTION 2026-09-29 (local, < 1 s) | `make_em3_data.py` (importe `EM2/em2_B_compare`) | `em3_ring_2p33.npz`, `em2_postw90_sigma_ti.npz`, `em_table.{md,csv}` → `scripts/fig/make_figures_em.py` → `fig_em_coupling` | actif |
| M4_sigma | `M4_sigma` (+ `pilote/`) | σ(ω) trois variantes, carte K, statistiques d'anneaux | PRODUCTION 2026-09-29 (local, ~4 min) | `m4_prod.py`, `make_table.py` → `em_table.tex` ; `pilote/{sweep,check_equiv,convergence}.py` | `em_sigma_*_N1200_eta0.04.npz`, `em_map_K.npz`, `em_table.tex` → mémoire §2.5 | actif (fermé) ; `em_table.csv` diffère de celui d'EM3 |
| B | `B_sigma_complex` | σ complexe (perspective B) | PRODUCTION 2026-09-29 (local, 3 min) | `b_prod.py`, `b_figure.py` (`_palette`) | `b_sigma_complex.{pdf,png,npz}` → hors mémoire | figé |

Tests associés : `tests/test_{tb_model,kgrid,velocity_operator,ring,kubo,diagnostics,wannier_io}.py` (§7).
`dev.ipynb` (racine, ignoré) : carnet brouillon EM (`from electron_photon import *`), jetable.

### 6.2 `campagnes/M/ch4/` et `memoire/défauts.tex` (chapitre 4, ajoutés le 2026-09-30, commits `66602ed`, `90d2022`)

Campagne « Chiffres du ch. 4 » : statut déclaré PRODUCTION, calcul local de quelques secondes, aucun job, seul
emplacement (pas de répertoire de travail hors dépôt). But : source unique pour réécrire le ch. 4 avec la production
finale (v2 + alignement de Kumagai–Oba, `results/M2_plateau`). Rien n'a été touché dans `src/`, `scripts/`, `config/`,
`results/`, `article/`, `figures/` ; seul fichier modifié hors du répertoire : `NOTES_TGAMMA.md`.

**Cinquième passe** : le tableau ci-dessous décrit l'état de `90d2022` ; ce qui a changé avec `214ed02` (partie 5) et
`6340c2f` (chapitre retiré), 12 fichiers suivis dans `campagnes/M/ch4/` :

| Fichier | Changement | Conséquence |
|---|---|---|
| `ch4_chiffres.py` (904 → 1 091 l.) | fonctions `compl_5_1` … `compl_5_4`, `complements`, `curve_scalars` ; appariement du tex refait (pourcentages × 100, notation scientifique lue en entier, entiers < 100 comparés dans les `tabular`) ; `all` réécrit la date et le HEAD des en-têtes | `check` rejoué à la cinquième passe : PASS (13 C_N, écart 0) ; `tex` (et donc `all`) cassé par `6340c2f`, **réparé** : lit le chapitre dans le dépôt du mémoire (`THESIS`/`TEX`, `$MSC_THESIS`, `--tex`), §11 |
| `table_v1_final.md` (631 → 738 l.) | section « Compléments » en fin de fichier : 5.1 tab:rcut_M en % (trois csv, 9×9 et 12×12, R_cut 0…6), 5.2 tableau niveau 1 (6 tailles × R_cut 0…4), 5.3 tab:échantillonnage (ΔE_F, V̄ à a_CC, aligné = non aligné − C_N), 5.4 scalaires lus sur les courbes des npz (Born/T min 1,050 / max 202,851 ; Γ_T à c = 0,1 % : 2,48 / 40,58 / 8,47 meV ; ħ/Γ 26 / 212 fs), porte PASS sur les 8 valeurs v1 et non alignées publiées | table principale (207 lignes) et sections a–l inchangées ; l'en-tête dit encore HEAD `90d2022` (généré avant le commit) |
| `defauts_nombres.md` (565 → 685 l.) | refait : **606 nombres** ; appariés **234**, non appariés 63, entiers courts hors `tabular` non comparés 309 ; lignes de la table jamais citées 22/207 (avant : 553 ; 87 / 67 / 399 ; 59/207) | régénéré après réparation sur l'original du dépôt du mémoire (HEAD `4da4d62`, md5 6d6c7c60…) : mêmes comptes, numéros de ligne +1 |
| `NOTES_TGAMMA_partie5.diff` (29 l., nouveau) | `git diff` des trois lignes de `NOTES_TGAMMA.md` | **redondant** avec l'annexe C de `ch4_rapport.md` et `git show 214ed02 -- NOTES_TGAMMA.md` |
| `ch4_rapport.md` (333 → 403 l.) | « Partie 5 — compléments » et annexe C ; constats : aucune ligne en double dans les trois `m_rcut_convergence.csv` ; `results/M/level1_summary.csv` n'a pas de R_cut 4 (le 2 468 du mémoire vient de `m_rcut_resigma.csv`) | se termine sur « STOP — Greg relit et commit » (fait : `214ed02`) |
| `README.md` | lignes `table_v1_final.md`, `defauts_nombres.md`, relance ; ligne `NOTES_TGAMMA_partie5.diff` | mentions de `memoire/défauts.tex` corrigées (dépôt du mémoire) ; dit toujours « `CLAUDE_md.diff` non appliqué » |
| `memoire/défauts.tex` | **supprimé** (`6340c2f`) | unique exemplaire : `~/projects/msc-graphene-raman-defects/memoire/chapitres/défauts.tex` (md5 6d6c7c60…, 750 l., toujours v1, 0 occurrence de « Kumagai ») ; les mentions de l'ancien chemin ne subsistent que dans le corps historique de `ch4_rapport.md` (addendum en fin de rapport) |

| Fichier | Rôle | Lit | Utilisé par / pour | Statut |
|---|---|---|---|---|
| `ch4/ch4_chiffres.py` (904 l.) | pilote unique, sous-commandes `check` (porte : 13 C_N de la config redonnés au bit par la règle ≥ 0,75 r_max), `table`, `regions`, `notes`, `tex`, `all` ; racine déduite de `__file__` ; **aucun import du paquet** (numpy, matplotlib ; `scripts/_palette` par `sys.path` l. 725 pour la figure) | `campagnes/R/R{4,5,6,8,9,10}` (tables, json, npz), `results/M2`, `results/M2_plateau` (csv, npz, `MD5SUMS_2026-09-30.txt`), `config/production.json`, `memoire/défauts.tex`, `wannier/27x27/wannier.wout` | relance : `.venv/bin/python campagnes/M/ch4/ch4_chiffres.py all` | **actif** ; `check` rejoué en local le 2026-09-30 : PASS (écart 0,0 eV) ; `regions` lit les `scf.in` des super-cellules **hors dépôt** (`/lustre09/project/…/graphene/qe/defects/super_cell/<N>/defective/scf.in`, chemins tirés d'`a1_results.json`) → rorqual seulement ; `notes` appelle `git diff` |
| `ch4/table_v1_final.md` (631 l.) | partie 1 : table principale v1 → non aligné → final, **207 lignes** (117 appariées, 9 v1 seulement, 81 final seulement ; statuts : remplacé 85, nouveau 87, retiré par décision 24, inchangé 6, sans équivalent final 5) ; 31 écarts entre les deux tables sources (rapportés, non corrigés) ; chiffres nouveaux a–l (porte A.2, escalier D4, convergence en bandes, niveaux QE π/σ par taille, chaîne repliée − QE, E_res vs grille, familles N mod 3, Kaasbjerg Fig. 3/13/14, C14 ±9,05 meV, Re Σ vs R_cut, ⟨ΔV⟩_3D, Γ^ed/Γ^ep du ch. 5) ; grandeurs retirées | `R6/etape3/table_v1_v2.md` (md5 f8972c4d…), `R10/c/table_v2_plateau.md` (md5 807bb352…) + sources ci-dessus | **source des chiffres du ch. 4** (et de la ligne Γ^ed/Γ^ep du ch. 5) ; cité par NOTES_TGAMMA.md et CLAUDE.md | actif (généré, ne pas éditer) |
| `ch4/alignement_regions.{md,npz,pdf,png}` | partie 2 : régions d'échantillonnage de l'alignement, 13 tailles : (A) ≥ 0,75 r_max (production), (B) Kumagai–Oba 2D ≥ N·a/2, (C) Kumagai–Oba 3D min(N·a/2, c/2), (D) site unique publié et site vraiment le plus loin ; figure au style du mémoire | `R10/a/profiles_<S>.npz`, `a1_results.json`, `scf.in` (a, c) | justification de la région dans le ch. 4 ; figure **non installée dans `figures/`** ni dans le mémoire | actif |
| `ch4/defauts_nombres.md` (565 l.) | partie 4 : les 553 nombres de `memoire/défauts.tex` appariés à la table : 87 appariés (valeurs **v1**), 67 non appariés (surtout des paramètres), 399 entiers courts non comparés ; 59 lignes de la table jamais citées | `memoire/défauts.tex` (md5 70140bc0…) | liste de travail de la réécriture du ch. 4 | actif ; à régénérer après chaque révision du chapitre |
| `ch4/NOTES_TGAMMA.diff` (180 l.) | partie 3 : `git diff` de `NOTES_TGAMMA.md` au moment de la session | — | trace | figé ; **redondant** avec `git show 90d2022 -- NOTES_TGAMMA.md` et avec l'annexe A de `ch4_rapport.md` |
| `ch4/CLAUDE_md.diff` (42 l.) | diff « proposé, non appliqué » de `CLAUDE.md` | — | trace | figé ; **périmé** : il a été appliqué tel quel dans `90d2022` (le README de `ch4/` dit encore « non appliqué ») |
| `ch4/ch4_rapport.md` (333 l.) | rapport de session : phase 0, porte, comptes, écarts entre sources, erreurs d'archives, annexe A (diff complet de NOTES_TGAMMA) | — | provenance | figé ; dit « rien n'est commité » et « fichiers non suivis » (vrai à l'écriture, faux depuis `90d2022`) |
| `ch4/README.md` | README de campagne (but, prompt, statut, date, table des fichiers, relance, règles) | — | — | à jour, sauf « `CLAUDE_md.diff` non appliqué » |
| `memoire/défauts.tex` (750 l.) | copie du chapitre 4 du mémoire (état v1 : aucune mention de Kumagai–Oba, médianes 2 524 / 2 509 / 2 597 meV, « quasi-annulation à ∼ −2,5 eV ») | — | lu par `ch4_chiffres.py tex` (lecture seule) | **copie de travail** : le rapport (erreur d'archive n° 4) note que le prompt la voulait hors dépôt ; identique à `~/projects/msc-graphene-raman-defects/memoire/chapitres/défauts.tex` à une ligne vide près (md5 70140bc0… vs 6d6c7c60…) → doublon à surveiller dès que l'un des deux est édité |

Décisions consignées par cette campagne (à reporter dans le mémoire) : alignement de Kumagai–Oba retenu, alignement à
site unique abandonné ; E_res retiré du ch. 4 ; C14 = ±9,05 meV autour du V_loc aligné ; tab:rcut_M à deux valeurs.
Constats ouverts recopiés dans NOTES_TGAMMA.md §8 : ⟨ΔV⟩_3D anormal pour 18×18, 24×24, 27×27 (dans le vide) ; 7
avertissements `c_bands` du nscf dense 27×27 ; le modèle replié « aligné » de R10 A.3 utilise C_9 à site unique ;
points 1–2 et 5 de R6 encore ouverts.

---

## 7. Les deux systèmes de tests

| Système | Fichiers | Ce qu'ils couvrent | Données | Lancement |
|---|---|---|---|---|
| **pytest** `tests/` (12 fichiers, 199 tests, collecte 0,31 s sans erreur) | `conftest.py` (fixtures jouet + `wannier/27x27` avec sha256 `W90_SHA256`, valeurs de référence `W90_REF`), `test_diagnostics` (16), `test_kgrid` (6), `test_kubo` (27, un test ~13 s), `test_ring` (65), `test_tb_model` (23), `test_velocity_operator` (28), `test_wannier_io` (9) — **EM** ; `test_matrix_io_units` (5, `tmp_path`) — chaîne M ; `test_r8_functions` (8, `disorder_average`, `pole_criterion`, `tb_models`, résolvante directe 6×6), `test_r9_functions` (4, `alignment`, `Mwr_to_Mwk_pairs`, `cluster_ldos`, `ldos_from_eigenpairs`), `test_r10_functions` (8, `ws_images`, `ws_phase`, `defect_mwr`, `recenter_mwr`) — synthétiques | modules EM, `matrix_io`, fonctions R8/R9/R10 de la chaîne T ; **pas** `compute_ML_*`, `non_local`, `wannier_interpolate`, `scattering_rate(_fast)` | aucune donnée externe hors `wannier/27x27` (versionné) | `.venv/bin/python -m pytest tests` ; aucun `sys.path`, aucun import de `scripts/` ni des pilotes |
| **scripts autonomes** `scripts/test_*.py` (8) + `validate_wannier_bands.py` | voir §4.6 | chaîne M (reconstruction KS, zero-padding, interpolation Wannier) et test d'or de la chaîne T | `.save` locaux 5×5/11×11 ou aucune ; 2 exigent le cluster | `.venv/bin/python scripts/test_X.py` ; PASS/FAIL, code 0/1 (sauf `test_local_green_batch` sans code, `test_local_tmatrix_real` toujours 0) |

Le renommage du dossier (quatrième passe) imposait `PYTHONPATH=src` ; le `.pth` est refait depuis la cinquième passe
(§13 étape 0) et `.venv/bin/python -m pytest tests` collecte de nouveau 199 tests (`.venv/bin/pytest` reste mort,
shebang périmé). CLAUDE.md l. 44 (« no pytest ») contredit l. 50 et l'existence de `tests/` ; sa liste l. 139 omet
`test_matrix_io_units`, `test_r8/r9/r10_functions`. `pyproject.toml` n'a pas de `[tool.pytest]` ; `pytest` et
`mpi4py` sont absents de `requirements.txt`.

---

## 8. `figures/` (72 fichiers suivis) et usage par le mémoire

35 noms × {pdf, png} + `hamiltonian_reconstruction.png` (sans producteur, `git mv` de l'étage 6) + `memoire.mplstyle`
(chargé par les 4 `make_figures*`, `_palette`, r4/r5/r8/r9/r10, EM2/EM3/M4/B ; `figure.figsize` = 6,5 × 4,0, pas
6,5 × 3,6 comme le dit CLAUDE.md l. 240).

Producteurs : `make_figures_memoire.py` (6 `_final`/`fig_convergence`/`fig_Ved`), `make_figures.py` (14),
`make_figures_epw.py` (8), `make_figures_em.py` (1). **Orphelines** (aucun script n'écrit ce nom) : les 6
`fig_epw_*_mv0.002` (12 fichiers, anciennes figures degauss 0.002 renommées à la main en P18 ; seule
`fig_epw_vs_ed_mv0.002` est citée, en v1) et `hamiltonian_reconstruction.png`.

Le dépôt du mémoire — **`~/projects/msc-graphene-raman-defects/memoire/`** depuis le 2026-09-30 (HEAD `e07e0d5` ;
chapitres dans `chapitres/` ; historique des emplacements au §0 bis) — inclut **20 figures**, et son `figures/` ne
contient plus que ces 20 PDF (liste et comparaison md5 revérifiées à `e07e0d5`, puis à `4da4d62` : inchangées, 8 diffèrent) :

| Figure (mémoire) | Chapitre | Producteur ici | md5 mémoire vs `figures/` (2026-09-30) |
|---|---|---|---|
| fig_ks_reconstruction, fig_epw_kohn_degauss, fig_epw_phonons (ch. 2), fig_epw_validation, fig_epw_gamma, fig_epw_g_control, fig_epw_phonselfen, fig_epw_decay, fig_Ved_zoom | défauts.tex / eph.tex / théorie.tex | `make_figures.py`, `make_figures_epw.py` | identiques |
| fig_M_map_final, fig_M_scaling_final, fig_Ved, fig_convergence, fig_level2, fig_locality_final, fig_spectral_final, fig_epw_vs_ed | défauts.tex, eph.tex | `make_figures_memoire.py`, `make_figures.py`, `make_figures_epw.py` | **diffèrent** : le mémoire a la version R6 (installée par `install_figures_R6.sh`), `figures/` a la version R10 du 2026-09-30 → **à réinstaller dans le mémoire** |
| KB_projectors_C.pdf, fig_ebands_edos, fig_electron_convergence (théorie.tex) | ch. 2 | **hors dépôt** : `~/projects/qe_pp` (`plot_ebands_edos.py`, `plot_convergence.py`) ; `KB_projectors_C.pdf` probable sortie manuelle de `plotting/plot_psp_radial_proj.py` (copie dans `results/`) | — |

Les figures non incluses qui traînaient dans le `figures/` du mémoire (`fig_band_interp`, `fig_phonons`, `fig_M_map`,
`fig_M_scaling`, `fig_locality`, `fig_spectral`, `fig_plateau`, `fig_rcut`, `fig_Ved_{boundary,map,profile_mean,radial,
radial_masked}`, les 6 `fig_epw_*_mv0.002`) en ont été **retirées** (19 PDF, commit `dfc3180` du dépôt du mémoire) ; celles
produites ici restent dans `figures/` de ce dépôt, `fig_band_interp` et `fig_phonons` ne subsistent que dans l'historique
git du mémoire et dans `~/projects/qe_pp`. `fig_em_coupling` (§2.5) n'y est pas encore (le texte du §2.5 n'est pas
écrit ; le mémoire dit encore que « l'évaluation numérique [du couplage électron-photon] est laissée hors du cadre »,
introduction.tex l. 72, raman.tex l. 297). Le registre `PROVENANCE.md` du dépôt du mémoire reprend ce tableau figure
par figure, avec le label LaTeX et la ligne du `.tex`.

À venir pour le ch. 4 (décidées, pas encore dans le mémoire ni dans `figures/`) : les cinq figures de R8
`campagnes/R/R8_kaasbjerg/fig/{dos_c,spectral_GKM,sigma_K,sensibilites,superposition}.pdf` (produites par
`r8_driver.py fig` et, pour `sigma_K`, `r8_driver.py sigeff` ; noms sans préfixe `fig_`) et, si elle est retenue, `campagnes/M/ch4/alignement_regions.pdf`
(`ch4_chiffres.py regions`). Le texte de `défauts.tex` (dépôt du mémoire, seul exemplaire depuis `6340c2f`) est encore à
l'état v1 (aucune occurrence de « Kumagai » ; 234 nombres appariés sur 606 depuis la partie 5, contre 87 sur 553) : la
réécriture s'appuie sur `campagnes/M/ch4/table_v1_final.md` (table principale + compléments 5.1–5.4) et `defauts_nombres.md`.

**Figures hors de `figures/` (cinquième passe ; Greg : « des figures dispersées dans `article/` et `memoire/` vont dans
le mémoire aussi »).** `figures/` n'est donc pas la seule source de l'étape F. Recensement des figures suivies ailleurs
(PDF + PNG), avec ce qui est déjà décidé ; la colonne de droite est **à arbitrer figure par figure** à l'étape F :

| Emplacement | Figures (noms sans extension) | Producteur | Pour le mémoire ? |
|---|---|---|---|
| `campagnes/R/R8_kaasbjerg/fig/` | `dos_c`, `spectral_GKM`, `sigma_K`, `sensibilites`, `superposition` | `r8_driver.py fig` / `sigeff` | **oui, décidé** (ch. 4 ; `241739e`, `PROVENANCE.md` l. 53) |
| idem | `7a_controle`, `dos_c_rho0_300`, `superposition_rho0_300` | `r8_driver.py` | non (hors portée, README R8) |
| `campagnes/M/ch4/` | `alignement_regions` | `ch4_chiffres.py regions` (rorqual) | « si retenue » : à trancher |
| `campagnes/EM/EM2/` | `em2_A_anneaux`, `em2_B_sigma` | `em2_A_figure.py`, `em2_B_figure.py` | validation du §2.5 : à trancher (annexe ?) ; `fig_em_coupling`, elle, est déjà dans `figures/` |
| `campagnes/EM/B_sigma_complex/` | `b_sigma_complex` | `b_figure.py` | non (perspective B, hors mémoire) |
| `campagnes/R/R10_plateau/fig/` (hors les 28 `fig_*` identiques à `figures/`) | `offset_profiles_13`, `levels_vs_invN`, `eres_vs_grid`, `kaasbjerg_plateau_ws`, `M_map_brut_aligne`, `dV_z_profiles` ; `avant_apres/` (28) | `r10_driver.py` | contrôles de l'alignement : à trancher (candidats naturels si le ch. 4 justifie Kumagai–Oba) |
| `campagnes/R/R9_controles/fig/` | `offset_profiles`, `resonance_vs_nkint{,_R9,_plateau}`, `rcut_aligned`, `folded_vs_R7`, `kaasbjerg_fig3_map` | `r9_driver.py` | contrôles ch. 4, « non installées » : à trancher |
| `campagnes/R/R7_tailles_3m/fig/` | `size_3m`, `localized_3m` (+ `_d1_8pts`, `_d1_reg`) | `r7_driver.py` | prévu pour l'article : à confirmer |
| `campagnes/R/R4_quasi_lie/fig/`, `R5_base_vs_M/fig/`, `R6_production_corrigee/etape2/fig/` | `d1_states`, `d3_alpha_*`, `d4_ladder`, `d5_boundary`, `d6_*` ; `a_delta`, `a_ladder_corrected`, `size`, `delta_vs_n`, `ladder_bands` ; `b3_ladder_M2`, `d3_pole_M2`, `d4_ladder_M2` | `r4_driver.py`, `r5_driver.py`, `r6_d4.py` | prévu pour l'article (R4 en v1) : à confirmer |
| `campagnes/R/R6_production_corrigee/etape3/figures_{v2,fix}/` | versions R6 des figures de `figures/` | `make_figures*` (R6) | non (remplacées par R10) |

Conséquence pour l'étape F : chaque figure retenue hors de `figures/` apporte son producteur (un pilote de campagne, pas
un `make_figures*`) et ses données (`out/`, json, npz de la campagne) ; à inscrire dans `PROVENANCE.md` « figures à
venir », qui ne liste aujourd'hui que les cinq R8, `alignement_regions` et `fig_em_coupling`.

Tables du mémoire (labels `tab:`) et leur source ici : `tests_M` ← `M_tests_summary.csv` ; `L_NL` ← `lnl_frobenius.csv` ;
`rcut_M` ← `m_rcut_convergence.csv` ; `échantillonnage` ← `sampling_table.csv` ; `convergence_gamma`, `tests_T`,
`param_T` ← `level1_summary.csv`, `level2_*.csv`, `nkint_check_9x9.csv`, NOTES_TGAMMA ; `param_M`, `param_dft`,
`param_wannier`, `param_dfpt`, `param_epw` ← `production.json`, `wannier/*/wannier.wout`, NOTES_EPW ; `gamma_ep_conv`,
`gamma_ph_conv`, `gamma_ph_prod`, `kohn_chaines` ← `results/epw/*.npz` (via `epw_*_post.py`) ; §2.5 (à venir) ←
`campagnes/EM/M4_sigma/em_table.tex`, `EM3/em_table.md`. Toutes sont saisies à la main dans le `.tex` (aucun `\input` de csv).

---

## 9. Fichiers de la racine, `config/`, `.claude/`

| Fichier | Rôle | Statut | Points périmés (ligne) |
|---|---|---|---|
| `CLAUDE.md` (25,8 Ko) | guide Claude Code | partiellement périmé ; **corrigé par `90d2022`** : l. 243 (`results_dir` = M2_plateau, `matrices_dir`, `results_dir_frozen`), l. 255 (`sampling_table.csv`), nouvelle ligne « Alignement du potentiel de défaut » (→ `campagnes/M/ch4/`), ligne `campagnes/M/ch4/` du tableau des campagnes ; le reste de la colonne de droite vaut toujours (l. 44, 59, 63, 131, 156, 217, 240 revérifiées sur `656c05c`) | l. 44 « no pytest » ; l. 59/63-65/217 `run.py`, `compute_M_cluster.py` ; l. 91-94 `use_ws_distance` sans la nuance R10 (`ws_images`) ; l. 97-143 carte des modules sans `config`, `matrix_io`, `projwfc_io`, `qe_gamma_io`, `wannier_provenance`, `alignment`, `deltav_pw`, `many_body/*`, `supercell_fold`, `sc_projection`, `electron_phonon`, `plotting` ; l. 131 `single_defect` « orphan » ; l. 139 tests incomplets ; l. 156 « no defective 12×12 » ; l. 157/203 `M_ed.npy` à la racine et « not covered by .gitignore » (faux : `*.npy` ignoré) ; l. 240 figsize ; l. 243 `results/M2` vivant ; l. 302 règle 4 manquante ; l. 314-328 tableau sans R2–R10 ni C ; l. 46-49 nom du dossier : **de nouveau exact** depuis le renommage local en `graphene-raman` (quatrième passe), et l. 38 et 50 (« already done in .venv », « The venv has the editable install ») sont de nouveau vraies depuis l'étape 0 (cinquième passe) |
| `README.md` (4,4 Ko, anglais) | vitrine GitHub | périmé | l. 18-19 M^L « in reciprocal space » (production = noyau réel) ; l. 24-25 « Born fails by 3–16 » (v1 ; v2 : médian 46,5) ; l. 52-59 layout sans `article/`, `memoire/`, avec `notebooks/` inexistant ; l. 56 `results/M` production ; l. 63-65 `data/` en liens (vrais fichiers en local) ; série EM absente |
| `CLEANUP.md` (23,3 Ko) | plan et journal du ménage P13 (rorqual) | clos ; en-tête périmé | l. 6-8 « dry-run… rien supprimé » ; l. 185 « non poussés » ; reste : 9 répertoires vides scratch (l. 235-244), `graphene_scf.save` 1,95 Gio, `phonons/graphene.wfcN` ~0,9 Gio, étage `codes/` éventuel (l. 154-156, nom « étage 7 » réattribué à R1) |
| `INVENTAIRE_2026-09-16.md` | ce fichier (ex-inventaire rorqual) | réécrit 2026-09-30 | — |
| `NOTES_TGAMMA.md` (45,1 Ko depuis `90d2022` ; 45,5 Ko depuis `214ed02` : les trois « non publiés par R10 » du §2 et de C13 — Born/T min et max, Γ_T à c = 0,1 %, ħ/Γ — sont remplis par les compléments 5.4) | chaîne t/Γ : étapes, paramètres, contrôles C1–C18, balayage N_k^int ; **état final R10** | **à jour sur le fond** (mis à jour le 2026-09-30 par la campagne `campagnes/M/ch4/`) : en-tête « Mise à jour R10 — état final », §1 (ligne « Alignement du potentiel » = Kumagai–Oba, 13 C_N, `defect_mwr`, images WS ; coûts R10), §2 et §3 en valeurs finales (`results/M2_plateau`), C14 redéfini, E_res retiré, nouveau §6e (N_k^int final), nouveau §8 (R10 : décision, code, portes, non aligné → final pour Re Σ et Kaasbjerg, Fig. 13–14 de R8, constats ouverts) | reste périmé : l. 6 `ab-initio-defects` ; §0 décrit encore la chaîne sans `defect_mwr` (recentrage seul) ; `fichier:ligne` du 2026-09-18 dérivés (ex. positivité `local_tmatrix.py:269-271` → l. 229/305) ; §4, §5, §7 = états historiques v1/R6 assumés ; §8 cite le hash réécrit `23ee3bb` ; A_k/DOS de R8 résumés au §8 seulement (pas d'étapes `disorder_average` au §0) |
| `NOTES_EPW.md` (53,7 Ko) | état de la chaîne EPW P0–P18, note R6 | fond à jour ; §5-§6 périmés | l. 4 ancien nom ; l. 219-222 `_ph0` restants (supprimés) ; l. 222/228 « commits à pousser » (poussés) ; l. 226 point ouvert : quelle chaîne (0.002 ou 0.02) le chapitre cite |
| `NOTES_EPW_REPERES.md` (4,9 Ko) | repères durables EPW | à jour sauf | l. 8-9 « production degauss 0.002 » (ambigu depuis P18, figures en mv0.02) ; l. 20 grilles restantes (étage 1 fait) |
| `REECRITURE_HISTORIQUE_2026-09-29.md` (38,8 Ko) | compte rendu de la réécriture git **exécutée** (Co-Authored-By retirés, courriels, 238 commits, bundle nearline) | fait et effectif | reste (l. 17-20) : hash dans les sorties json/logs et le code gelé R6/R7/R9 ; `23ee3bb` de R10 non traduit |
| `requirements.txt` | pip freeze de l'ère ABINIT (nov. 2025) | périmé | `abipy`, `netCDF4` (interdits par CLAUDE.md l. 206) ; `-e git+…/ab-initio-electron-defect-interaction` ; `numpy==1.26.4` vs 2.3.5 installé ; `ase`, `spglib` listés mais absents du venv ; `mpi4py`, `pytest` absents |
| `pyproject.toml` | métadonnées minimales (0.1.0, src layout) | fonctionnel | aucune dépendance déclarée, pas de `[tool.pytest]` |
| `.gitignore` | ignore `data/`, `*.npy`, `*.save/`, `.venv`, `cleanup/`, `dev.ipynb`, `results/*` sauf listes blanches (M, M2, M2_plateau : `specwd_*_prod`, `resonance_*`, `mwr_locality`, `M_analysis`, `ks_reconstruction`, `ved_analysis`, `*.csv` ; + `MD5SUMS_*`, `ed_vs_ep_*`, `README.md` pour M2/M2_plateau ; `results/epw/*.npz,*.csv`) ; `figures/` suivi en entier | à jour | `results/M/README.md` non couvert ; `notebooks/` non listé (CLAUDE.md l. 203) |
| `LICENSE` (MIT 2026) | — | à jour | — |
| `dev.ipynb` | brouillon EM, ignoré | jetable | — |
| `.claude/settings.local.json` | 4 permissions Bash, `attribution` vide (pas de trailer co-auteur) | à jour | — |
| `reports/Optique du graphène avec lacunes.md`, `research_notes/…/{ab_initio_methods,experiments,spin_vacancy,theory_models}.md` | revue bibliographique du 2026-09-29 (base du PLAN C) | à jour | — |

Chevauchements : règles de campagne (CLAUDE.md l. 290-313 ↔ CLEANUP l. 12-23, 194, 199) ; conventions EPW écrites
trois fois (CLAUDE.md l. 267-275, NOTES_EPW §3, NOTES_EPW_REPERES) ; conventions de figures (CLAUDE.md l. 224-247,
NOTES_EPW l. 199, REPERES l. 11) ; R6 raconté quatre fois ; unités Ha→eV (NOTES_TGAMMA l. 53, CLAUDE.md l. 246-247,
NOTES_EPW l. 201) ; familles N mod 3 (NOTES_TGAMMA C12, CLAUDE.md l. 249-255) ; description du projet (README ↔
CLAUDE.md l. 5-22).

---

## 10. Ce qui sert au mémoire (à copier plus tard dans `~/projects/msc-graphene-raman-defects`)

Le dépôt du mémoire a maintenant quatre dossiers d'accueil vides (README seulement) : `analyse/` (code de
post-traitement, « copie figée du repo graphene-raman »), `calculs/` (inputs QE, Wannier90, EPW et scripts SLURM, un
dossier par calcul), `donnees/` (données des figures), `scripts_figures/`. Correspondance proposée au §13, étape F.
**À ajouter à la liste du chapitre 4 ci-dessous** : `campagnes/M/ch4/` en entier (pilote, table v1 → final, régions
d'alignement, nombres du chapitre) et les cinq figures R8 de la portée du mémoire (`dos_c`, `spectral_GKM`, `sigma_K`,
`sensibilites`, `superposition`) avec leurs données `campagnes/R/R8_kaasbjerg/out/{dos,spec,sigeff,sens,7a}/`.

**Chapitre 4 (M, matrice T, Γ)** — code : `src/graphene_raman/{config,io/qe_io,io/matrix_io,
io/pseudo_io,io/wannier_io,io/wannier_provenance,utils/fft_utils,utils/lattice,utils/planewaves,wavefunctions/wfk,
wavefunctions/fold_wfk_to_sc,wannier/wannier_hamiltonian,wannier/wannier_interpolation,defects/local_R,
defects/non_local,defects/alignment,defects/deltav_pw,wavefunctions/sc_projection,defects/many_body/local_tmatrix,
defects/many_body/pole_criterion,defects/many_body/single_defect}.py` ; scripts de production `compute_M.py`,
`compute_M_dense_stages.py`, `assemble_M2.py`, `gate_M_normalization.py`, `compute_spectral_wannier.py`,
`resonance_metrics.py`, `resonance_criteria.py`, `rcut_resigma.py`, `nkint_check_post.py`, `m_rcut_convergence.py`,
`mwr_locality_coarse_vs_dense.py`, `analyze_M.py`, `analyze_Ved.py`, `ks_reconstruction_all.py`, `sampling_table.py`,
`lnl_frobenius_all.py`, `level2_families.py`, `finalize_wannier.py` ; figures `make_figures.py`, `make_figures_memoire.py`,
`_palette.py`, `figures/memoire.mplstyle` ; tests `test_ks_reconstruction.py`, `test_zero_pad_dense.py`,
`test_pad_vs_full_supercell.py`, `test_wannier.py`, `test_local_tmatrix{,_real}.py`, `test_local_rcut.py`,
`test_local_green_batch.py`, `tests/test_{matrix_io_units,r9_functions,r10_functions}.py` ; lanceurs
`submit_spectral_wannier_dense.sh`, `submit_rcut_resigma.sh`, `submit_nkint_check.sh`, `submit_lnl_frobenius.sh`,
`submit_post.sh`, `submit_M.sh`, `submit_M_dense.sh`, `campagnes/R/R10_plateau/submit_r10.sh` ; données `config/production.json`,
`results/M2_plateau/*`, `wannier/{24,25,27,28,32}x*` (sans `_nb16`) ; pilotes `campagnes/R/R10_plateau/r10_driver.py` (base
vivante), `campagnes/R/R8_kaasbjerg/r8_driver.py` + `defects/many_body/disorder_average.py` (DOS, A_k), R6 (provenance de v2).
Documentation : NOTES_TGAMMA.md §0–§3, §6e, §8 (état final R10 depuis `90d2022`), `campagnes/M/ch4/table_v1_final.md`, CLAUDE.md l. 72-95 (conventions), l. 249-265
(N mod 3, sous-réseaux), `campagnes/R/R10_plateau/c/table_v2_plateau.md`, `R10_rapport.md`.

**Chapitre 5 (EPW)** — `src/…/electron_phonon/{phself,selfen}.py` ; `scripts/epw_{pp_save,extract_gkk,validate,
selfen_post,phself_post,d2_extract,ring_check,dfpt_path_freq,phdos_extract,ed_vs_ep}.py`, `make_figures_epw.py`,
`submit_epw_p1_post.sh`, `submit_epw_p2_post_mv.sh` ; `results/epw/*.npz` ; NOTES_EPW.md §1c–§1g, §3, §4,
NOTES_EPW_REPERES.md.

**§2.5 (électron-photon)** — `src/…/electron_photon/*.py` ; `tests/test_{tb_model,kgrid,velocity_operator,ring,kubo,
diagnostics,wannier_io}.py` + `conftest.py` ; `campagnes/EM/{EM.md,EM1_tb,EM2,EM3,M4_sigma}` ; `scripts/fig/make_figures_em.py` ;
`wannier/27x27` ; `io/wannier_io.read_w90_tb`.

Hors dépôt mais nécessaires au ch. 2 : `~/projects/qe_pp` (`plot_ebands_edos.py`, `plot_convergence.py`,
`notebooks/{phonons,epw}.ipynb`) pour `fig_ebands_edos`, `fig_electron_convergence`, `fig_phonons`, `fig_band_interp`.

---

## 11. Cassé, obsolète, incohérent (liste exacte)

**Arbitrage du §11 (Greg, 2026-09-30) — réparations faites, non commitées ; les listes plus bas décrivent l'état d'avant.**

| Fichier | Décision | Fait |
|---|---|---|
| `scripts/validation/test_local_tmatrix_real.py` | réparer | lit les M dans `matrices_dir` ; sort en code 0 (PASS) ou 1 (FAIL) ; non rejouable en local (scratch, nommage cluster), syntaxe vérifiée seulement |
| `scripts/m/tag_vacancy_sublattice.py` | réparer | cherche les sidecars dans `matrices_dir` |
| `scripts/slurm/submit_M.sh`, `submit_M_dense.sh` | réparer | écrivent les `M_*` dans `matrices_dir` (`MAT`), plus dans `results_dir` ; `bash -n` OK ; non soumis |
| `scripts/slurm/submit_golden_dense.sh` | réparer | code de sortie = celui du test d'or ; l'appel à `check_M_dense_vs_coarse.py` (obsolète) est retiré, `analyze_M.py` fait ce contrôle |
| `check_ML_coarse_kernel.py`, `check_M_dense_vs_coarse.py`, `check_onsite_and_NL.py` | obsolètes (absorbés par `analyze_M.py`, `mwr_locality_coarse_vs_dense.py`) | non modifiés ; iront dans `scripts/_obsolete/` à l'étape A |
| `_mcheck.py`, `_normtest.py`, `validate_ML_grid_7x7.py`, `check_M_dense_nb20_vs_nb16.py` | obsolètes (non rejouables) | idem |
| `campagnes/R/R9_controles/r9_driver.py`, `campagnes/R/R6_production_corrigee/etape3/r6_compare_v1_v2.py`, dérives de sens R4/R5/R6 | **en suspens** : politique des campagnes figées à arbitrer (archives rejouables à leur commit, pas à HEAD ?) | rien |

**Cassé par le renommage local du dossier (quatrième passe), réparé à la cinquième** : l'installation éditable du venv
(`.pth` → `/Users/gregou/Projects/raman-graphene/src`, inexistant) faisait échouer tout import du paquet sans
`PYTHONPATH=src` (pytest, `scripts/test_*.py`, `make_figures{,_memoire,_epw}.py`, pilotes EM). Étape 0 faite : il ne
reste que les shebangs morts de `.venv/bin/pytest` et `py.test`.

**Cassé par le retrait de `memoire/défauts.tex` (`6340c2f`, cinquième passe), réparé le même jour (non commité)** :
`campagnes/M/ch4/ch4_chiffres.py tex` et `all` plantaient (`FileNotFoundError`, chemin fixe `ROOT/memoire/défauts.tex`). Le
pilote lit maintenant `memoire/chapitres/défauts.tex` dans le dépôt du mémoire (`THESIS` = dépôt voisin
`msc-graphene-raman-defects`, sinon `$MSC_THESIS`, sinon `--tex CHEMIN` ; arrêt avec message si absent). `tex` et `table`
relancés : mêmes comptes (606 ; 234 / 63 / 309 ; 22/207), porte 5.4 PASS ; dans `defauts_nombres.md` seuls l'en-tête et
les numéros de ligne changent (+1 partout : l'original a une ligne vide de plus en tête). README et rapport (addendum)
mis à jour. `notes` relancé aussi : 15/15 valeurs finales retrouvées (il ne réécrit `NOTES_TGAMMA_partie5.diff` que si
`NOTES_TGAMMA.md` a un diff non commité). Reste hors de portée en local : `regions`, donc `all` (rorqual seulement).

**Cassé à l'exécution (chemins R10 : `results_dir` = M2_plateau sans matrices)** :
`scripts/check_ML_coarse_kernel.py:5,8-9` · `check_M_dense_vs_coarse.py:8,15-17` · `check_onsite_and_NL.py:9,16-17,32-33,44` ·
`check_M_dense_nb20_vs_nb16.py:3` (déjà non rejouable) · `_mcheck.py:4,7,11` · `_normtest.py:11` (+ mélange Ha/eV l. 11-12) ·
`test_local_tmatrix_real.py:18,27,29` (+ FAIL sort en code 0, l. 63) · `validate_ML_grid_7x7.py:16` (+ argv l. 9, garde l. 15) ·
`tag_vacancy_sublattice.py:12` (silencieux, 0 sidecar) · `submit_golden_dense.sh:25` ·
`campagnes/R/R9_controles/r9_driver.py:855,903` (rejeu) · `campagnes/R/R6_production_corrigee/etape3/r6_compare_v1_v2.py:14`
(TypeError : `results_dir_frozen` est une liste).
Correction commune : `matrices_dir(cfg)` ou `dense_paths(cfg, S)["mfile"]` à la place de `results_dir(cfg)` pour les
`M_*.npy`. **Écriraient les M au mauvais endroit** : `submit_M.sh:40-42`, `submit_M_dense.sh:25-33`.

**Dérive de sens (lisent M2_plateau au lieu de M2 « tel quel »)** : `r6_level1_gate.py`, `r6_m_rcut_resigma.py`,
`r6_tests_gate_row.py`, `r9_driver.py:607,724,973,995` ; rejeu R4/R5 mélange v1 (`RES = results/M`) et v2
(`dense_paths` → M2) : `r4_driver.py:42,109,216-217`, `r5_driver.py:171`.

**Obsolètes (désactivés ou non rejouables)** : `compute_spectral.py`, `compute_tmatrix.py`, `compute_convergence.py`,
`_eta_scan.py`, `migrate_M_norm.py`, `_old_vs_new_7x7.py`, `check_M_dense_nb20_vs_nb16.py`,
`check_M_dense_vs_coarse_bands.py` (orphelin), `summarize_level1_maps.py`, `submit_spectral.sh`, `submit_tmatrix.sh`,
`submit_spectral_wannier.sh`, `campagnes/R/R9_controles/cloture/{r9_driver_propose.py,submit_r9_propose.sh}`,
`defects/local_G.py` (+ étiquetage v2 douteux via `--kernel G`), `io/pseudo_io.read_psp8`,
`local_tmatrix.scattering_rate_from_wannier`, `utils/interpolation.py`, `plotting/`, `wannier/*_nb16`,
`results/test_recon`, `results/M/M_ed.npy`, `results/epw/ed_vs_ep_*.npz` (v1), `figures/fig_epw_*_mv0.002.*`,
`figures/hamiltonian_reconstruction.png`.

**Incohérences documentaires** (§9) et **valeurs périmées dans les lanceurs** : `submit_post.sh:31` (v1), `:36` (±25 meV),
`make_figures_epw.py` (défauts de tags), `finalize_wannier.py:85-90`, `link_data.sh:2`, docstrings de `test_wannier`,
`validate_wannier_bands`, `matrix_io` (défaut `SUPERCELL`), `local_R.compute_ML_R` et `local_G.prep_reciprocal_inputs`
(`subtract_mean=True` par défaut ≠ production), `scattering_rate_fast` (`ne_per_eta=4` ≠ config 8),
`m_rcut_convergence.py` (csv en ajout), `Mwk_to_Mwr` (commentaire de forme l. 87). Lignes 1-2 avant le shebang :
`compute_M`, `compute_convergence`, `compute_spectral`, `compute_tmatrix`, `epw_ed_vs_ep`, `migrate_M_norm`,
`nkint_check_post`, `make_figures_epw`.

**Chemins cluster codés en dur** (n'empêchent rien sur rorqual, bloquent en local) : `config.dense_paths` (scratch),
`compute_M_dense_stages.py:24`, `check_M_dense_vs_coarse{,_bands}.py`, `check_onsite_and_NL.py`,
`test_local_tmatrix_real.py:26`, `link_data.sh`, `epw_d2_extract.py:28`, `epw_ring_check.py:14`, `finalize_wannier.py:22`,
`submit_epw_p{1,2}_post*.sh`, tous les pilotes `article/` (`PROJ`, `R4DIR`, `R5DIR`), `campagnes/EM/EM2/*` (`PROJECTS`).
Nommage local `defect_unit_cell_<N>.save` (`_paths.py`) ≠ cluster `defect_<N>.save` (`link_data.sh`, scripts de production).

---

## 12. Doublons

| Sujet | Fichiers | Remarque |
|---|---|---|
| Figures de travail vs finales | `make_figures.py` ↔ `make_figures_memoire.py` (helpers `lab`, `famlab`, `panel`, `save`, `load_map`, `bz_vertices`, `FAM3` recopiés ; fig_rcut+fig_plateau ↔ fig_convergence ; fig_locality ↔ `_final` ; fig_spectral ↔ `_final` ; fig_M_map/scaling ↔ `_final` ; fig_Ved_map+radial_masked ↔ fig_Ved) | `memoire` dépend de `make_figures` pour fig_Ved_zoom |
| k communs dense/grossier | `check_M_dense_vs_coarse.py`, `check_M_dense_vs_coarse_bands.py`, `analyze_M.py:173-183` | un seul vivant (analyze_M) |
| noyau dense p=1 vs grossier | `check_ML_coarse_kernel.py`, `analyze_M.py:198-205` | idem |
| nbnd 16→20 | `check_M_dense_nb20_vs_nb16.py`, `analyze_M.py:185-196` | non rejouable |
| rapport L/NL | `lnl_frobenius_all.py`, `_bz_ratio_LNL.py`, `analyze_M.py:116-128` (`mmap_M`, `wannier_V`, `pi_pair`, `box_dirichlet` recopiés) ; `check_onsite_and_NL.py` Q2 | |
| sur-site p_z | `check_onsite_and_NL.py` Q1 ↔ `mwr_locality_coarse_vs_dense.py` | |
| V_ed le long de a₁ et détection de la lacune | `analyze_M.py:132-147` ↔ `analyze_Ved.py` §2 ; argmax dmin recopié 3× alors que `alignment.vacancy_site` existe | |
| assemblage de M | `compute_M.py combine`, `compute_M_dense_stages.py combine`, `assemble_M2.py` | trois voies |
| bandes Wannier vs DFT | `validate_wannier_bands.py`, `compare_bands_qe.py`, `compare_bands_w90_qe.py`, `finalize_wannier.py` (erreur de bandes), `epw_validate.py` (bloc bandes) | |
| tests de la matrice t | `test_local_tmatrix`, `test_local_rcut`, `test_local_tmatrix_real` (référence `compute_T × N_k` recopiée) | |
| parseurs EPW | `read_plot` (d2_extract) ↔ `read_epw_plot` (validate) ; parseur prt ph.out (d2_extract) ↔ `read_dfpt_prt` (validate) ; prtgkk (ring_check ↔ extract_gkk) | |
| Γ^ed vs Γ^ep | `epw_ed_vs_ep.py` ↔ `make_figures_epw.py` fig_epw_vs_ed | |
| `HA2EV` | `config.py:5`, `matrix_io.py:26`, 8 scripts | |
| tables de tailles (p, D) | `config["dense"]` ↔ `PFAC` (`compute_M_dense_stages:23`), `D` (`check_M_dense_vs_coarse:10`, `_bands:6`), `PF` (`check_onsite_and_NL:11`, `test_local_tmatrix_real:23`) | |
| code promu de R5 | `campagnes/R/R5_base_vs_M/r5_{deltav_pw,sc_projection}.py` ↔ `src/…/defects/deltav_pw.py`, `wavefunctions/sc_projection.py` (copies de campagne gardées, encore importées par R6) | |
| données | `results/M/ved_analysis.npz` = `results/M2/ved_analysis.npz` ; `ks_reconstruction.npz` M2 = M2_plateau ; `results/M2/*.csv` = `R6/etape3/csv_v2/*.csv` ; 33 figures `figures/` = `R10_plateau/fig/` ou `R6/etape3/figures_v2/` ; `em_table.csv` M4 ≠ EM3 (contenus différents, même nom) | |
| pilotes proposés vs finaux | `R9_controles/cloture/r9_driver_propose.py` ↔ `r9_driver.py` | |
| chapitre 4 (seconde passe) | `memoire/défauts.tex` ↔ `~/projects/msc-graphene-raman-defects/memoire/chapitres/défauts.tex` (une ligne vide d'écart, revérifié à `e07e0d5`) | la copie du dépôt servait d'entrée à `ch4_chiffres.py tex` ; **doublon résolu à la cinquième passe** : copie supprimée (`6340c2f`), reste à faire lire l'original par le pilote (§11) |
| diffs de session (seconde passe) | `campagnes/M/ch4/NOTES_TGAMMA.diff` ↔ annexe A de `ch4_rapport.md` ↔ `git show 90d2022 -- NOTES_TGAMMA.md` ; `campagnes/M/ch4/CLAUDE_md.diff` ↔ `git show 90d2022 -- CLAUDE.md` (appliqué) ; cinquième passe : `campagnes/M/ch4/NOTES_TGAMMA_partie5.diff` ↔ annexe C du rapport ↔ `git show 214ed02 -- NOTES_TGAMMA.md` | trois copies de chaque diff |
| provenance des figures et tables du mémoire (troisième passe) | §8 et §10 de ce fichier ↔ `~/projects/msc-graphene-raman-defects/PROVENANCE.md` (§1 figures, §2 tables) | même information en deux endroits : ce fichier décrit l'état d'avant ménage, `PROVENANCE.md` deviendra la référence après la copie (colonne « Dans msc- ») |
| tables de chiffres du ch. 4 | `R6/etape3/table_v1_v2.md` (v1 → v2), `R10/c/table_v2_plateau.md` (v2 → final), `campagnes/M/ch4/table_v1_final.md` (fusion, 207 lignes), NOTES_TGAMMA.md §2–§3 | la table de `campagnes/M/ch4/` fait foi pour la réécriture ; 31 écarts de forme entre les deux sources y sont listés |

---

## 13. Proposition de réorganisation (à arbitrer étape par étape ; rien n'est fait)

Principe : ne rien déplacer qui casse un import ou un lanceur gelé sans mettre à jour l'appelant dans le même
commit ; les campagnes `article/` et `memoire/` restent des copies figées (leurs `sys.path` pointent vers `scripts/`
et vers les répertoires rorqual) ; `figures/` reste la cible d'installation.

### Étape 0 — préalable : réparer l'installation éditable (ajout de la quatrième passe) — **FAITE** (Greg, avant 12:31)

Vérifié à la cinquième passe : `.pth` → `/Users/gregou/projects/graphene-raman/src`, `config.ROOT` correct, 199 tests
collectés sans `PYTHONPATH`. Non fait (facultatif) : recréer le venv, seule façon de réparer `.venv/bin/pytest`. Texte
d'origine, pour mémoire :

Depuis le renommage du dossier en `graphene-raman`, le paquet ne s'importe plus en local (§0 bis). Avant toute autre
étape (chaque déplacement des étapes A, C, E se valide par `pytest` et par les scripts autonomes) :
`.venv/bin/python -m pip install -e . --no-deps --no-build-isolation` depuis la racine (setuptools 80.9 est dans le
venv ; réécrit le `.pth` et `direct_url.json`, ne touche aucun fichier suivi), puis `.venv/bin/python -m pytest tests`
(199 tests attendus). Variante : recréer le venv (Python 3.13) et en profiter pour régénérer `requirements.txt`
(étape B), ce qui répare aussi les shebangs morts de `.venv/bin/`.

**Arbitrage (Greg, 2026-09-30)** : une étape à la fois, GO explicite par étape ; périmètre = réorganisation seulement
(tableau « points hors étapes » en fin de section) ; git en lecture seule pour Claude, Greg commite.

### Étape A — EXÉCUTÉE le 2026-09-30 (GO de Greg sur les trois choix : suppression des obsolètes, découpage `m/`/`t/`, `epw_ed_vs_ep` dans `t/`)

Fait : `scripts/{m,t,epw,fig,validation,slurm}/` (7 + 14 + 9 + 4 + 12 + 11 fichiers), helpers `_palette.py`, `_bands.py`, `_paths.py`
gardés à la racine de `scripts/` (les pilotes de `campagnes/` les importent par `sys.path`), index `scripts/README.md` ; **21 scripts
obsolètes supprimés** (liste dans l'index ; dernier commit qui les contient : `d9a3588`). Réparations : racine déduite de `__file__` un
niveau plus bas (`assemble_M2`, `gate_M_normalization`, `finalize_wannier`, `make_figures_em`, trois `epw_*`) ; `sys.path` vers les
helpers dans les 4 `make_figures*` et 5 scripts de validation ; `ks_reconstruction_all` → `validation/test_ks_reconstruction` ;
`run_test_A_batch` ; tous les appels `scripts/<nom>` réécrits dans les 11 lanceurs, `link_data.sh`, `submit_r10.sh`, `r9_driver.py`
(+ `scripts/m` dans son `sys.path`), `em2_B_compare.py`, `src/` (3 docstrings), les README de `results/`, `campagnes/README.md`,
CLAUDE.md, README, les trois NOTES, `EM.md`, `PROVENANCE.md` du mémoire. **Greg a réorganisé `figures/` entre-temps** (`e03e706` :
`electron/`, `electron_defect/`, `electron_phonon/`, `electron_photon/`) : les défauts `--outdir` des quatre `make_figures*` pointent
maintenant sur le bon sous-dossier, CLAUDE.md mis à jour. Vérifié : `bash -n` des lanceurs, `py_compile`, tests autonomes locaux
(PASS), `make_figures_em.py` régénère `fig_em_coupling.png` au md5 identique, `ch4_chiffres.py check`, pytest 199. Reste : sur rorqual,
`sbatch scripts/slurm/<lanceur>.sh` au lieu de `scripts/<lanceur>.sh`.

Proposition d'origine, pour mémoire :

### Étape A — `scripts/` en sous-dossiers (81 → 6 familles)

| Sous-dossier | Contenu | Précautions |
|---|---|---|
| `scripts/m/` (chaîne M) | `compute_M.py`, `compute_M_dense_stages.py`, `assemble_M2.py`, `gate_M_normalization.py`, `tag_vacancy_sublattice.py`, `finalize_wannier.py`, `link_data.sh`, `_diag_mnl_mpi.py`, `_test_mnl_mpi.sh` | `r9_driver` importe `compute_M_dense_stages` par `sys.path` → mettre à jour ou laisser un stub ; lanceurs `submit_M*.sh`, `submit_r6_kernel.sh` |
| `scripts/t/` (chaîne T) | `compute_spectral_wannier.py`, `resonance_metrics.py`, `resonance_criteria.py`, `rcut_resigma.py`, `nkint_check_post.py`, `m_rcut_convergence.py`, `mwr_locality_coarse_vs_dense.py`, `level2_families.py`, `analyze_M.py`, `analyze_Ved.py`, `ks_reconstruction_all.py`, `sampling_table.py`, `lnl_frobenius_all.py`, `epw_ed_vs_ep.py` | `submit_r10.sh`, `submit_post.sh`, `runbook_3.sh` codent `scripts/<nom>` ; `ks_reconstruction_all` importe `test_ks_reconstruction` |
| `scripts/epw/` | les 9 `epw_*.py` (sauf `epw_ed_vs_ep`) | `submit_epw_p{1,2}_post*.sh` |
| `scripts/fig/` | `make_figures.py`, `make_figures_memoire.py`, `make_figures_epw.py`, `make_figures_em.py`, `_palette.py`, `_bands.py` | `_palette` est importé par r4/r5/r8/r9/r10, EM2/EM3/B via `sys.path.insert(…, "scripts")` → soit un stub `scripts/_palette.py`, soit déplacer `_palette` dans `src/…/plotting/` (et `memoire.mplstyle` avec lui) |
| `scripts/validation/` | `test_ks_reconstruction.py`, `test_wannier.py`, `test_zero_pad_dense.py`, `test_pad_vs_full_supercell.py`, `test_local_tmatrix.py`, `test_local_rcut.py`, `test_local_green_batch.py`, `test_local_tmatrix_real.py`, `validate_wannier_bands.py`, `compare_bands_qe.py`, `compare_bands_w90_qe.py`, `run_test_A_batch.py`, `_paths.py` | `test_local_rcut` importe `test_local_tmatrix` ; `r6_kernel_check` importe `test_ks_reconstruction` |
| `scripts/slurm/` | les 14 `submit_*.sh` | `submit_r10.sh` (article) appelle `scripts/submit_*.sh` par chemin |
| `scripts/_obsolete/` (ou suppression après arbitrage) | `compute_spectral.py`, `compute_tmatrix.py`, `compute_convergence.py`, `_eta_scan.py`, `_normtest.py`, `_mcheck.py`, `migrate_M_norm.py`, `_old_vs_new_7x7.py`, `check_M_dense_nb20_vs_nb16.py`, `check_M_dense_vs_coarse_bands.py`, `summarize_level1_maps.py`, `submit_spectral.sh`, `submit_tmatrix.sh`, `submit_spectral_wannier.sh`, `_bz_ratio_LNL.py` | cités par NOTES_TGAMMA, CLEANUP, rapports R6 : garder une ligne de renvoi |
| réparés ou classés (arbitrage du §11, 2026-09-30) | réparés : `test_local_tmatrix_real.py` (→ `validation/`), `tag_vacancy_sublattice.py` (→ `m/`), `submit_golden_dense.sh`, `submit_M*.sh` (→ `slurm/`) ; vers `_obsolete/` en plus de la ligne ci-dessus : `check_ML_coarse_kernel.py`, `check_M_dense_vs_coarse.py`, `check_onsite_and_NL.py`, `validate_ML_grid_7x7.py` | — |

Alternative plus légère : ne créer que `fig/`, `epw/`, `validation/`, `slurm/`, `_obsolete/` et laisser M + T à plat
(les ~25 scripts de production que `submit_r10.sh` appelle).

### Étape B — `.md` de la racine — **EXÉCUTÉE le 2026-09-30** (GO de Greg par lot ; non commitée ; pytest : 199 passés)

| Lot | Fait | Reste / à faire par Greg |
|---|---|---|
| B1 `CLAUDE.md` | corrigé en place (26,7 → 36,8 Ko) : section « Repository layout » ; commandes réelles (`compute_M.py`, `compute_M_dense_stages.py`, pytest, scripts autonomes) ; carte des modules complétée (`config`, `matrix_io`, `wannier_provenance`, `alignment`, `deltav_pw`, `local_tmatrix`, `disorder_average`, `pole_criterion`, `supercell_fold`, `electron_phonon` ; `local_G` historique ; `single_defect` = référence des tests d'or) ; données locales réelles ; nuance `ws_images` ; `figsize` 6,5 × 4,0 ; chemins `notes/`, `admin/` ; règle « Greg commite, git en lecture seule » ; tableau des campagnes avec R2–R10 et C | règle des campagnes figées (pilotes à garder ou non) : en attente d'arbitrage ; règle 4 toujours absente de la numérotation |
| B2 `README.md`, dépendances | README mis à jour (noyau en espace réel, alignement, Born/T ≈ 46, électron-photon, layout réel, installation) ; dépendances déclarées dans `pyproject.toml` (numpy, scipy, h5py, matplotlib, tqdm ; extras `mpi`, `test`, `campaigns`) ; **`requirements.txt` supprimé** | — |
| B3 journaux → `admin/` | `CLEANUP.md` et `REECRITURE_HISTORIQUE_2026-09-29.md` déplacés, en-tête « clos » / « fait » ajouté | ce fichier reste à la racine jusqu'à la fin du ménage, puis `admin/CARTOGRAPHIE_2026-09-30.md` (et `PROVENANCE.md` du mémoire à mettre à jour) |
| B4 notes → `notes/` | les trois NOTES déplacées ; `NOTES_TGAMMA.md` : chemins relatifs, hash `23ee3bb` → `c2bc733`, paragraphe « Depuis R10 » au §0 (`defect_mwr`, `matrices_dir`, `disorder_average`) ; `NOTES_EPW.md` : en-tête, état du §5 ; `NOTES_EPW_REPERES.md` : en-tête, mention de la chaîne 0.02 ; `ch4_chiffres.py notes` lit `notes/NOTES_TGAMMA.md` (15/15) | **rorqual, après `git pull`** : refaire les deux liens `graphene/qe/epw/NOTES_EPW.md` → `notes/NOTES_EPW.md` et `graphene/qe/epw/CLAUDE.md` → `notes/NOTES_EPW_REPERES.md` ; `PROVENANCE.md` du mémoire cite `GR:NOTES_TGAMMA.md` et `GR:NOTES_EPW.md` (l. 70, 79, 80, 92) → `GR:notes/…` ; chaîne EPW citée par le ch. 5 : **0.02, tranché par Greg** (CLAUDE.md, NOTES_EPW §5, REPERES corrigés) ; les mentions par simple nom dans les rapports de campagne (« NOTES_TGAMMA §2 ») restent valables et ne sont pas touchées |
| hors lot | `reports/` et `research_notes/` → `article/C_optique_lacunes/biblio/` (rapport + `research_notes/` à plat) ; `PLAN.md` l. 5 corrigé | `R10_rapport.md` l. 697 cite les anciens dossiers dans une phrase historique (laissée) ; dernière ligne de `PLAN.md` (« ImportError ») toujours périmée |

Proposition d'origine, pour mémoire :

| Fichier | Proposition |
|---|---|
| `CLAUDE.md` | mettre à jour (§9) : commandes réelles, carte des modules complète (dont la chaîne T), `results_dir`/`matrices_dir`, tableau des campagnes R2–R10 + C, tests pytest ; retirer ce qui est dans les NOTES (conventions EPW, règles de ménage → renvoi) |
| `README.md` | réécrire court : layout réel, v2/R10, série EM, `article/`, `memoire/` |
| `NOTES_TGAMMA.md` | **fait pour l'essentiel le 2026-09-30** (`90d2022` : état final R10, §6e, §8) ; reste : mettre le §0 au niveau de `defect_mwr` et des étapes DOS/A_k de R8, rafraîchir les `fichier:ligne`, corriger l. 6 et le hash `23ee3bb` ; §4/§5/§7 → journal ; c'est la source de tab:param_T / tab:tests_T avec `campagnes/M/ch4/table_v1_final.md` |
| `NOTES_EPW.md`, `NOTES_EPW_REPERES.md` | corriger §5-§6 et l. 8-9/20 ; trancher la chaîne citée (0.002 vs 0.02) ; fusionner REPERES dans NOTES_EPW §3 ou l'inverse (une seule copie des conventions) |
| `CLEANUP.md`, `REECRITURE_HISTORIQUE_2026-09-29.md`, `INVENTAIRE_2026-09-16.md` | journaux d'administration : déplacer dans `admin/` (ou `docs/admin/`) avec un en-tête « clos le … » ; ce fichier devient `admin/CARTOGRAPHIE_2026-09-30.md` |
| `requirements.txt` | régénérer depuis le venv (sans abipy/netCDF4, avec mpi4py, pytest, ase/spglib en optionnel) ou déclarer les dépendances dans `pyproject.toml` et supprimer le fichier |
| `reports/`, `research_notes/` | regrouper sous `article/C_optique_lacunes/biblio/` (ils sont la base du PLAN C) |

### Étape C — EXÉCUTÉE le 2026-09-30 (GO de Greg ; non commitée)

Fait : `[tool.pytest.ini_options]` dans `pyproject.toml` (`testpaths = tests`, marqueurs `slow`, `needs_data`, `cluster`) ;
`tests/conftest.py` : `scripts/` et `scripts/validation/` sur `sys.path`, fixture `repo_cwd`, helper `local_data` ;
`tests/test_scripts_tmatrix.py` (golden synthétique, R_cut, `local_green_batch` [slow], test d'or réel [cluster, skip]) et
`tests/test_scripts_M.py` (zero-pad exact, non-régression dense [needs_data, slow], reconstruction KS [needs_data, slow], pipeline
Wannier 11×11 [needs_data], pad vs super-cellule [cluster, skip]) : **208 tests collectés** (199 + 9). Codes de sortie corrigés :
`test_local_green_batch.py` refondu en `main()` (seuil 1e-2 sur fast vs exact, mesuré 5,4e-3) ; `test_local_tmatrix_real.py` déjà
fait au §11. Marqués `slow` : `test_kubo_pauli_edges_real`, `test_sigma_complex_real` (5–8 s, les deux plus longs de la suite EM).
Constat : `test_zero_pad_dense` test 2 (dense 5×5, > 10 min) a planté une fois en local avec `MPIDI_OFI_handle_cq_error … OFI poll
failed (default nic=en12)` : MPICH/libfabric sur l'interface réseau du portable, pas la physique ; à relancer (ou `FI_PROVIDER`).
Les scripts de `scripts/validation/` restent lançables seuls ; CLAUDE.md et `scripts/README.md` mis à jour.

Proposition d'origine, pour mémoire :

### Étape C — deux systèmes de tests

1. Garder **pytest** comme unique lanceur : ajouter `[tool.pytest.ini_options] testpaths = ["tests"]` dans
   `pyproject.toml`, marquer les tests lents (`test_kubo` 13 s).
2. Envelopper les scripts autonomes exécutables sans données (`test_local_tmatrix`, `test_local_rcut`,
   `test_local_green_batch`, `test_zero_pad_dense` test 1) dans `tests/test_tmatrix_golden.py` et
   `tests/test_zero_pad.py` (appel de leurs fonctions, assert sur le résultat) ; corriger d'abord les codes de sortie.
3. Marquer `@pytest.mark.needs_data` (skip si `data/graphene` ou le scratch manque) les scripts qui exigent des
   `.save` : `test_ks_reconstruction`, `test_wannier`, `test_zero_pad_dense` test 2, `test_pad_vs_full_supercell`,
   `test_local_tmatrix_real` (après réparation `matrices_dir`).
4. Les scripts `scripts/validation/` restent comme outils de diagnostic (figures, rapports) mais ne portent plus le
   PASS/FAIL de référence ; CLAUDE.md l. 44 corrigée.

### Étape D — `results/` et données

- Ajouter `results/M/README.md` (copie de `R6/phase0/README_results_M_gele.md`) à la liste blanche ; décider du sort
  de `results/M/ved_analysis.npz` (doublon 13,5 Mio) et de `results/epw/ed_vs_ep_*.npz` (v1) ; supprimer
  `results/test_recon`, `results/M/M_ed.npy` (local) ; `wannier/*_nb16` (4 dirs, aucun lecteur) → nearline ou suppression.
- Réinstaller dans le mémoire les 8 figures R10 qui diffèrent (§8) et y ajouter `fig_em_coupling` quand le §2.5 sera écrit.
- Traduire `23ee3bb` → `c2bc733` dans `campagnes/R/R10_plateau/{README.md,R10_rapport.md}` et `NOTES_TGAMMA.md` §8 (ou noter
  la table de REECRITURE).
- `campagnes/M/ch4/` : corriger le README (« `CLAUDE_md.diff` non appliqué ») ; décider si `NOTES_TGAMMA.diff` et
  `CLAUDE_md.diff` restent (traces redondantes avec git ; de même `NOTES_TGAMMA_partie5.diff` depuis `214ed02`). Le sort
  de `memoire/défauts.tex` est **tranché** (copie supprimée, `6340c2f`) : `ch4_chiffres.py tex` lit
  désormais le chapitre dans le dépôt du mémoire (fait, §11).

### Étape D — EXÉCUTÉE le 2026-09-30 (GO de Greg ; non commitée ; pytest : 199 passés)

| Fait | Détail |
|---|---|
| `results/M/README.md` | créé depuis `R6/phase0/README_results_M_gele.md` (+ mise à jour R10), ajouté à la liste blanche de `.gitignore` |
| doublons supprimés | `results/M/ved_analysis.npz` (= `results/M2/`), `results/M2/ks_reconstruction.npz` (= `results/M2_plateau/`) ; notés dans les README de `results/M` et `results/M2` |
| obsolètes supprimés (suivis) | `results/epw/ed_vs_ep_24k24q{,_mv0.02}.npz` (v1 ; l'écrivain écrit dans `results_dir`), `wannier/{25,27,28,32}x*_nb16/` (24 fichiers, 25 Mo, aucun lecteur), `figures/fig_epw_*_mv0.002.{pdf,png}` (12, chaîne 0.002 abandonnée), `figures/hamiltonian_reconstruction.png` (sans producteur, non inclus) |
| locaux non suivis supprimés | `results/test_recon/`, `results/M/M_ed.npy` (juin, orphelin) |
| `campagnes/M/ch4/` | `NOTES_TGAMMA.diff`, `NOTES_TGAMMA_partie5.diff`, `CLAUDE_md.diff` retirés (redondants avec git et le rapport) ; README corrigé |
| non fait, à décider | (a) les `results/M/M_ed_{5…12}.npy` v1 locaux (238 Mo, non suivis ; sur rorqual `results/M/` sert encore à `assemble_M2.py`) ; (b) réinstaller les 8 figures R10 dans le mémoire : reporté à l'étape F, avec la réécriture du ch. 4 ; (c) `results/KB_projectors_C.pdf` et 4 PNG non suivis (sorties locales de scripts, inoffensifs) ; (d) hash `23ee3bb` dans `R10/README.md` et rapport : archives, laissés |
| constaté | Greg a mis `article/` dans `.gitignore` (`d9a3588`) : `C_optique_lacunes/` (plan et biblio) n'est plus suivi, il reste en local |

### Étape E — EXÉCUTÉE le 2026-09-30 (GO de Greg ; décisions : `results_dir()` pour R8, `plot_psp` → nouveau script, PDF reproductibles ; non commitée)

| Sous-étape | Fait |
|---|---|
| E1 paquet | `src/electron_defect_interaction` → **`src/graphene_raman`** (`pyproject.toml` `name = "graphene_raman"`, réinstallation éditable) ; 104 fichiers réécrits (py, sh, toml, md vivants, `PROVENANCE.md`) ; 4 fichiers d'archives/données touchés par la substitution restaurés depuis HEAD |
| E2 palette et style | `graphene_raman/plotting/palette.py` (ex `scripts/_palette.py` ; + `use_style()`, `save()` reproductible, `STYLE`) et `plotting/memoire.mplstyle` (ex `figures/memoire.mplstyle`, copie unique, `package-data`) ; 19 importeurs réécrits, plus aucun `sys.path` vers `scripts/` ; `_paths.py`, `_bands.py` → `scripts/validation/` |
| E3 résultats par la config | `r8_driver.py` : `md5_listed` lit les `MD5SUMS_*.txt` de `matrices_dir()`, la porte P2 lit `resonance_9x9.npz` dans `results_dir()` (base finale, commentaire : le run de 2026-09-29 comparait M2 tel quel) ; nouvelle clé `epw_results_dir` + `config.epw_dir()` : les 8 scripts `epw_*` et `make_figures_epw.py` n'ont plus de `results/epw` en dur |
| E4 orphelins | supprimés : `utils/interpolation.py`, `plotting/these.mplstyle`, et 12 fonctions mortes (`compute_ML_G_mpi`, `scattering_rate_from_wannier`, `split_counts`, `local_slice`, `read_psp8`, `monkhorst_pack_grid`, `generate_mp_grid`, `write_kpoints`, `make_Cdicts_for_k`, `get_typat`, `read_pdos_m`, `kvec_to_rvec`, constante `RY2EV`), chacune vérifiée sans utilisateur ; `plotting/plot_psp_radial_proj.py` remplacé par **`scripts/fig/make_figures_electron.py`** (fig:KB du ch. 2 : V_PPL vs −Z/r, projecteurs β_{iℓ}(r) avec r_c^(s) = 1,24 et r_c^(p) = 1,30 bohr, depuis `C.upf`) → `figures/electron/fig_kb_pseudo_C` **redessinée, à valider par Greg** (original dans git) |
| E5 défauts | `HA2EV` défini une seule fois (`config.py` ; `matrix_io` et 5 scripts l'importent) ; défauts `load_M_checked(require_bloch_norm=UNIT_CELL)`, `scattering_rate_fast(ne_per_eta=8)`, `compute_ML_R(subtract_mean=False)`, `prep_reciprocal_inputs(subtract_mean=False)` |
| E6 recette | pytest 201 (hors slow) après chaque sous-étape ; `ch4_chiffres.py check` ; les 29 figures régénérables (`make_figures{,_memoire,_epw,_em}.py`) sont **identiques pixel à pixel** aux versions suivies (les octets diffèrent : encodeur PNG et métadonnées de l'ancienne génération) ; **deux exécutions successives donnent des PDF et PNG identiques au bit** (métadonnées fixées dans `palette.save`) — c'est le test de recette de l'étape F ; `results/M2_plateau/*.csv` réécrits à l'identique par `make_figures.py` |
| reste | dossier `src/electron_defect_interaction.egg-info/` (ignoré, périmé) à supprimer à la main ; tags git `memoire-avant-menage` (sur `656c05c`) et `memoire-apres-menage` : Greg |

Proposition d'origine, pour mémoire :

### Étape E — `src/`

- Supprimer ou isoler : `utils/interpolation.py`, `plotting/` (après avoir régénéré `KB_projectors_C.pdf` une fois et
  consigné la commande), `local_G.py` (ou corriger son étiquetage v2 et le garder comme référence des tests de
  zero-padding, qui sont ses seuls utilisateurs), `read_psp8`, `scattering_rate_from_wannier`, fonctions mortes listées §2.
- Ranger `tb_models.py` sous `tests/` (banc synthétique) si R4 n'est plus rejoué.
- Une seule définition de `HA2EV` ; défaut `load_M_checked(require_bloch_norm=UNIT_CELL)` ; défaut `ne_per_eta=8`.

### Étape F0 — mise au style des figures de contrôle — EXÉCUTÉE le 2026-09-30 (décision de Greg : une seule migration)

`scripts/fig/make_figures_controles.py` : les fonctions de tracé de `r7_driver.figs`, `r9_driver.{b_figure,cmd_synth,c_figure}` et
`r10_driver.{a1_figure,a2_figure,c2_figure}` extraites **à contenu identique** (mêmes séries, panneaux, légendes), style du paquet,
lecture des json/npz de `campagnes/R/{R7,R9,R10}` → `figures/electron_defect/fig_{size_3m,localized_3m,resonance_vs_nkint,
resonance_vs_nkint_plateau,rcut_aligned,folded_vs_R7,offset_profiles_13,levels_vs_invN,kaasbjerg_plateau_ws}`. Retirés ensuite :
`r7_driver.py`, `r9_driver.py`, `r4_driver.py`, `r5_{driver,sc_projection,basis_diagnostics,alignment_ext}.py` (gardés jusque-là pour
ces figures), les copies des 9 figures dans `campagnes/R/*/fig/`, et `io/projwfc_io.py` (plus aucun utilisateur). `fig_kb_pseudo_C`
corrigée (trace r·β tel que tabulé, comme l'original ; validée par Greg). Les 5 figures R8 restent produites par `r8_driver.py` dans
`campagnes/R/R8_kaasbjerg/fig/` (données `out/`).

### Étape F — alimenter le dépôt du mémoire `~/projects/msc-graphene-raman-defects` (ajout de la seconde passe, chemin mis à jour à la troisième)

Le dépôt du mémoire a été préparé le 2026-09-30 (`5debac2`, puis README, `requirements.txt` et venv jusqu'à `e07e0d5`)
avec quatre dossiers vides. Correspondance proposée avec le §10 ; à faire **après** les étapes A–E pour copier un état
propre, et après réparation des chemins du §11. À chaque copie, remplir la colonne « Dans msc- » de `PROVENANCE.md`
avec le chemin d'arrivée et le commit source de ce dépôt.

| Dossier du mémoire | Ce qui y va (depuis ce dépôt) | Remarque |
|---|---|---|
| `analyse/` (« copie figée du repo graphene-raman ») | `src/graphene_raman/` (sans les orphelins du §2), `config/production.json`, `tests/`, les scripts de production des chaînes M et T et EPW (§4.1 à §4.3, §4.5 hors obsolètes), `campagnes/R/R10_plateau/r10_driver.py`, `campagnes/R/R8_kaasbjerg/r8_driver.py`, `campagnes/M/ch4/ch4_chiffres.py`, `campagnes/EM/{M4_sigma,EM3,EM2}/*.py` | noter le commit d'origine dans le README ; `pyproject.toml` pour l'installation |
| `calculs/` (un dossier par calcul) | hors dépôt pour l'essentiel : `graphene/qe/defects/{super_cell,unit_cell}` (scf, nscf denses, pp.x, Wannier90), `graphene/qe/epw/24k-24q{,_mv0.02}`, `graphene/qe/electron_photon/{EM1_tb,EM2}` sur rorqual ; depuis ce dépôt : `scripts/submit_*.sh` actifs, `campagnes/R/R10_plateau/submit_r10.sh`, `campagnes/EM/EM1_tb/wannier.win`, `campagnes/EM/EM2/{bands.in,bands_pp.in,submit_*.sh,postw90/}` ; l'annexe D du mémoire (`annexes/paramètres.tex`) liste les paramètres | les inputs QE ne sont pas versionnés ici : à prendre dans les répertoires de travail |
| `donnees/` | `results/M2_plateau/*` (npz, csv, MD5SUMS, README), `results/epw/*.npz` (chaîne `_mv0.02` + celles que `fig_epw_kohn_degauss` lit), `campagnes/EM/M4_sigma/{em_scalars.json,em_sigma_*.npz,em_map_K.npz,em_table.tex}`, `campagnes/EM/EM3/*.npz`, `campagnes/EM/EM2/em2_A.npz`, `campagnes/R/R8_kaasbjerg/out/` (figures R8), `campagnes/M/ch4/{table_v1_final.md,alignement_regions.*}`, `wannier/27x27` (ou les 5 grilles de production) | ~60 Mo de npz ; `results/M2/*.npy` (matrices) restent sur rorqual/nearline |
| `scripts_figures/` | `scripts/fig/make_figures.py`, `make_figures_memoire.py`, `make_figures_epw.py`, `make_figures_em.py`, `_palette.py`, `figures/memoire.mplstyle` ; les sous-commandes `fig` et `sigeff` de `r8_driver.py` ; `ch4_chiffres.py regions` ; hors dépôt : `~/projects/qe_pp` (`plot_ebands_edos.py`, `plot_convergence.py`) pour les figures du ch. 2 | adapter `results_dir` → `donnees/` par un seul argument ou une variable |
| `memoire/figures/` | réinstaller les 8 figures R10 qui diffèrent ; ajouter les 5 figures R8, `fig_em_coupling`, éventuellement `alignement_regions`, et les figures de `article/` et `memoire/` retenues à l'arbitrage (tableau « Figures hors de `figures/` », §8) | voir §8 ; chaque figure retenue entraîne son pilote (`scripts_figures/` ou `analyse/`) et ses données (`donnees/`) |
| `requirements.txt` du mémoire | le régénérer depuis un venv qui fait tourner `analyse/` et `scripts_figures/` (numpy, scipy, h5py, matplotlib, plus `mpi4py` si les noyaux de M sont copiés ; LaTeX pour `text.usetex`) | l'actuel est un gel généraliste (abipy, torch, netCDF4) sans rapport avec ces scripts |

À corriger en même temps : toute note qui cite `~/LaTeX/master_thesis` (vide depuis 11 h 21) ou
`~/LaTeX/msc-graphene-raman-defects` (supprimé) ; le skill `thesis-section-pass` est corrigé (cinquième passe).
`PROVENANCE.md` est commité (`4da4d62`) ; il nomme déjà ce dépôt `graphene-raman`, mais cite ce fichier à HEAD
`656c05c` : à porter au commit du ménage.

### Étape « campagnes » — EXÉCUTÉE le 2026-09-30 (décision de Greg, solution 3 ; non commitée)

Règle de Greg : pour les figures du MÉMOIRE on garde ce qui permet de les régénérer et de les modifier ; pour les figures
HORS MÉMOIRE on supprime tout (figures et pilotes de diagnostic) ; on garde les résultats, rapports, tables, inputs. Les
deux listes de figures (37 noms mémoire, 180 fichiers hors mémoire) sont celles de Greg ; toutes les figures suivies de
`article/` et `memoire/` y figurent (0 inconnue, 0 absente). Fait : **213 fichiers supprimés** (`rm`, git verra des `D`) :
180 figures/PDF hors mémoire (dont les 28 planches `R10/fig/avant_apres/`, les 69 + 6 fichiers de `R6/etape3/figures_{v2,fix}/`,
`R8/ref/kaasbjerg_2020….pdf`) et 33 pilotes/lanceurs (R1 : 4 py ; R2 : 2 py ; R4 : `submit_r4.sh` ; R5 : 2 sh + `r5_deltav_pw.py` ;
R6 : 8 py + 6 sh ; R7 : 3 py + 3 sh ; R9 : `submit_r9.sh`, `cloture/*_propose*`). Gardés : `r8_driver.py`, `submit_r8.sh`,
`r10_driver.py`, `submit_r10.sh` (vivants), `r9_driver.py` (4 figures du mémoire depuis les json), `r7_driver.py` (2 figures
du mémoire depuis `d1_8pts/d1_results.json`) **et, pour cette seule raison, ses imports** `R4/r4_driver.py`,
`R5/r5_{driver,sc_projection,basis_diagnostics,alignment_ext}.py` (à retirer quand les deux figures R7 seront refaites au
style du mémoire à l'étape F). Chaque README de campagne reçoit en tête une note de ménage (retiré / gardé / `git show
32efd9e:…`) et une **fiche** (question, pourquoi, méthode, résultat, interprétation/décision, où sont les données) tirée du
rapport. `campagnes/EM/*` et `campagnes/M/ch4/` ne sont pas touchés (figures toutes au mémoire ou perspective B). pytest : 199 passés.
Conséquences : les rejeus « cassés » de R4/R5/R6/R9 du §11 n'ont plus d'objet (pilotes retirés ou gardés seulement pour des
figures depuis json) ; à l'étape E, `projwfc_io`, `read_w90_HR`, `deltav_pw`, `sc_projection`, `qe_gamma_io`, `supercell_fold`,
`tb_models`, `pole_criterion` n'ont plus que R7/R8/R9/R10, `scripts/` et `tests/` comme utilisateurs (compter avant de supprimer).

### Points hors étapes, à arbitrer aussi (liste de Greg du 2026-09-28, reprise à la quatrième passe)

Ces points figuraient dans la liste d'origine du ménage et ne sont couverts par aucune des étapes A–F :

| Point | État vérifié | Ce que cela engage |
|---|---|---|
| Périmètre : simple réorganisation, ou **extraction du paquet « défauts seulement »** de l'article (`defects/` dont `many_body/`, `wannier/`, `io/`, `utils/`, `wavefunctions/` ; hors périmètre `electron_photon/`, `electron_phonon/`, pilotes et figures) | **tranché (Greg, 2026-09-30) : réorganisation seulement ; extraction du paquet après le mémoire** | les étapes A et E ne renomment ni ne découpent le paquet ; le renommage `defects` → `electron_defects` est reporté avec l'extraction |
| Renommer le sous-paquet `defects` → `electron_defects` | 52 fichiers importent `config`, 19 `local_tmatrix`, etc. (§2) ; les pilotes figés `campagnes/R/R4–R10` importent `defects.*` | casse le rejeu des campagnes figées sauf alias ; à lier au point précédent |
| Wrapper `wannier_io.read_w90_HR` | gardé pour `campagnes/R/R4`, `campagnes/R/R9` (§2) | le supprimer casse ces deux pilotes |
| `campagnes/R/R8_kaasbjerg/ref/kaasbjerg_2020_prb101_045433.pdf` (3,1 Mo, suivi) | article sous droits d'auteur dans un dépôt GitHub | **retiré le 2026-09-30** (étape « campagnes ») ; reste dans l'historique |
| Données de test | les tests ne lisent que 4 fichiers de `wannier/27x27/` (sha256 dans `tests/conftest.py`) | les regrouper sous `tests/data/` ; Zenodo pour les autres wannierisations si le dépôt devient public |
| Historique git | pack de 103 Mo, clone rorqual | **ne pas** réécrire (pas de `filter-repo`) |
| Nom de ce fichier | `INVENTAIRE_2026-09-16.md` pour un contenu du 2026-09-30 ; cité sous ce nom par `PROVENANCE.md` l. 6 | renommage prévu à l'étape B (`admin/CARTOGRAPHIE_2026-09-30.md`) ; mettre `PROVENANCE.md` à jour en même temps |
