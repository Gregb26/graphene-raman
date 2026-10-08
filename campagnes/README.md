# campagnes/ — copies versionnées des campagnes de calcul (règle 5 de CLAUDE.md)

Un sous-dossier par série, un répertoire par campagne, chacun avec son README (fiche : question, pourquoi, méthode, résultat,
interprétation, où sont les données). Depuis le 2026-09-30 (ménage) : ex-`article/R*` → `R/`, ex-`memoire/EM` → `EM/`,
ex-`memoire/ch4` → `M/ch4`. Le chemin n'encode ni le statut ni la destination : ils sont ici et dans les fiches.
Les pilotes de diagnostic des campagnes closes ont été retirés (règle de Greg, 2026-09-30) ; ils restent dans l'historique
git (`git show 32efd9e:article/<campagne>/<fichier>`, `32efd9e:memoire/…` pour EM et ch4).
Renommage local → cluster de la chaîne T (2026-09-30 ; `local_tmatrix` → `cluster_tmatrix`, `extract_V_loc` → `cluster_potential`,
`local_green(_batch)` → `cluster_green(_batch)`, `local_t` → `cluster_t`, `pole_criterion.local_t_cache` → `cluster_t_cache`) : les
pilotes `R/R8_kaasbjerg/r8_driver.py` et `R/R10_plateau/r10_driver.py` gardent les anciens noms (leur copie suit le répertoire de
travail sur rorqual) ; pour les relancer, `git checkout 1d67419` (dernier commit avant le renommage), ou mettre à jour le pilote dans
son répertoire de travail puis resynchroniser la copie (`cp -p`).

| Série / campagne | But | Statut | Destination | Pilote gardé ? |
|---|---|---|---|---|
| `R/R1_vacancy_relaxed` | relaxation de la lacune 9×9, nspin 1/2 ; contrôles 3×3×1 (R1b) | archive (PRODUCTION, miroir md5) | article | non |
| `R/R2_size_series` | relaxations 5×5…12×12 | archive (PRODUCTION, miroir md5) | article | non |
| `R/R4_quasi_lie` | états π/σ de la lacune 9×9, DFT vs chaîne M → Wannier → T (v1) | archive (TEST) | article | non (retiré le 2026-09-30, F0) |
| `R/R5_base_vs_M` | base ou M ? constat A.2 (facteur N_cells) | archive (TEST) | article | non (retirés le 2026-09-30, F0) |
| `R/R6_production_corrigee` | M2 = N_cells·M^L + M^NL, porte A.2, régénération du ch. 4 | archive (PRODUCTION ; produits remplacés par R10) | mémoire (provenance de v2) | non (`scripts/m/assemble_M2.py`) |
| `R/R7_tailles_3m` | famille 3m 15…27, état π et doublet σ vs taille ; R7c relaxations | archive (PRODUCTION, miroir md5) | mémoire (2 figures) + article | non : ses 2 figures → `scripts/fig/make_figures_controles.py` |
| `R/R8_kaasbjerg` | DOS et A_k moyennées sur le désordre (Kaasbjerg Fig. 13/14) | **vivante** (TEST) | mémoire (5 figures) + article (étape 6) | `r8_driver.py`, `submit_r8.sh` |
| `R/R9_controles` | alignement C_N, résonance vs N_k^int, chaîne repliée, Kaasbjerg, audit | archive (TEST, clos 2026-09-29) | mémoire (4 figures de contrôle) | non : ses 4 figures → `scripts/fig/make_figures_controles.py` |
| `R/R10_plateau` | base unique du ch. 4 : Kumagai–Oba 13 tailles, rejeu de la production → `results/M2_plateau` | **vivante** (base du ch. 4, close 2026-09-30) | mémoire | `r10_driver.py`, `submit_r10.sh` (ses 3 figures de contrôle → `make_figures_controles.py`) |
| `EM/EM1_tb` | r(R) de la wannierisation 27×27 (`restart = plot`) | archive (PRODUCTION) | mémoire §2.5 | `em1_check.py` |
| `EM/EM2` | références indépendantes : DFT directe aux anneaux, postw90 `kubo` | archive (PRODUCTION) | mémoire §2.5 (2 figures) | oui (`em2_B_compare.py` importé par EM3) |
| `EM/EM3` | données de la figure et tableau du §2.5 | **vivante** (PRODUCTION, local) | mémoire §2.5 | `make_em3_data.py` |
| `EM/M4_sigma` | σ(ω), cartes et chiffres des anneaux | **vivante** (PRODUCTION, local) | mémoire §2.5 | `m4_prod.py` |
| `EM/B_sigma_complex` | σ(ω) complexe (perspective B) | archive (PRODUCTION, local) | hors mémoire (1 figure gardée) | `b_prod.py`, `b_figure.py` |
| `M/ch4` | chiffres du ch. 4 : table v1 → final, régions d'alignement, nombres de `défauts.tex` | **vivante** (PRODUCTION, local) | mémoire ch. 4 | `ch4_chiffres.py` |

Plan de la série EM : `EM/EM.md`. Les plans des trois articles (A méthode + SCTMA + MoS₂, B spin, C σ(ω)) sont dans le dépôt `../articles/` (`ROADMAP.md`) ;
leurs campagnes viendront ici sous `SCTMA/`, `MOS2/`, `SPIN/`, `OPT/`.
