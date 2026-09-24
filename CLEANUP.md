# CLEANUP — sauvegarde scratch et grand ménage (P13, 2026-09-17)

Suite de l'inventaire P12 (`INVENTAIRE_2026-09-16.md`). Périmètre :
`~/links/projects/rrg-cotemich-ac/gregb26/` et `/scratch/gregb26`.

État au 2026-09-17 : **phase 0 exécutée et vérifiée ; phases 1–2 en DRY-RUN
(manifestes seulement, aucun `rm`) ; phase 3 en attente du GO.** Rien n'a été
supprimé, déplacé ni renommé. Chaque étage attend un GO explicite, séparé.

---

## Règle d'hygiène (à retenir)

> **Rien d'unique ne vit uniquement sur le scratch.** En fin de campagne QE,
> `rsync` des `prefix.save/` de production vers
> `graphene/qe/qe_tmp_backup/` (même arborescence que `qe_tmp/`), sans les
> `prefix.wfcN` en vrac, puis md5 source/destination. Le scratch reste la copie
> de travail (les `outdir` et les liens `ab-initio-defects/data/` y pointent) ;
> la purge Alliance (60 jours sans accès ni modification) ne prévient pas.
> Commande type (rsync ≥ 3.2, `--open-noatime` pour ne pas rafraîchir les
> atimes du scratch ; sans `-g`/`-o` pour que les fichiers prennent le groupe
> `rrg-cotemich-ac` du projet, sinon quota) :
> `rsync -a -r --no-o --no-g --open-noatime --files-from=liste.txt /scratch/gregb26/qe_tmp/ graphene/qe/qe_tmp_backup/`

---

## Phase 0 — EXÉCUTÉE : sauvegarde des données uniques du scratch

Miroir : `graphene/qe/qe_tmp_backup/` (structure de `qe_tmp/` conservée,
mtimes préservés, groupe `rrg-cotemich-ac`). Copié le 2026-09-17 entre 08:36 et 08:47.

| Bloc | Répertoires | Fichiers | Taille (Gio) | md5 |
|---|---|---|---|---|
| `defect_{5..12}x{5..12}_{p,d}/*.save` (SCF supercellules) | 16 | 64 | 38,1 | OK |
| `defect_uc_dense_{24,25,27,28,32}/*.save` (NSCF denses nbnd=20) | 5 | 3 753 | 11,0 | OK |
| `defect_unit_cell_{5,6,7,8,9,10,11,12}x…/*.save` (mailles unitaires grossières) | 8 | 704 | 1,6 | OK |
| `defect_unit_cell_{5x5,7x7,8x8}_wann/*.save` | 3 | 147 | 0,3 | OK |
| `defect_unit_cell_{5x5,7x7,8x8}_wann/` (reste du dossier : wannier.*, wfc en vrac, xml) | 3 | 192 | 0,4 | OK |
| **Total** | **32 `.save` + 3 dossiers** | **5 007** | **51,8** | **5 007 / 5 007 identiques** |

Non copiés, volontairement : `qe_tmp/graphene_scf/` (SCF maille unitaire,
reproductible), `qe_conv/` (convergences, résultats déjà dans le projet), les
`prefix.wfcN` en vrac hors `.save` (copies par processus MPI), les `prefix.xml`
à côté des `.save` (doublons de `data-file-schema.xml`).

Vérification :
- tailles apparentes source = destination pour chacun des 32 `.save`
  (53 626 735 Kio des deux côtés) et des 3 dossiers `_wann` ; nombres de fichiers identiques (table par répertoire ci-dessous) ;
- md5 de chaque fichier, source lue avec `O_NOATIME` (script python, pas
  `md5sum`), destination normale : **0 mismatch, 0 erreur de lecture** ;
- atimes du scratch après copie et md5 : inchangés (`defect_9x9_d.save/wfc1.hdf5`
  2026-06-18, `defect_12x12_p.save/wfc1.hdf5` 2026-06-18, `defect_5x5_d.save/charge-density.hdf5`
  2026-06-17) — aucun contournement de purge ;
- fichier de sommes : `graphene/qe/qe_tmp_backup/MD5SUMS_2026-09-17.txt`
  (5 007 lignes, format `md5  chemin_relatif`, vérifiable par
  `cd graphene/qe/qe_tmp_backup && md5sum -c MD5SUMS_2026-09-17.txt`).

Détail par répertoire (Gio, nombre de fichiers, source → destination) :

```
defect_10x10_d.save   3.14/4 → 3.14/4    defect_10x10_p.save   3.16/4 → 3.16/4
defect_11x11_d.save   4.52/4 → 4.52/4    defect_11x11_p.save   4.54/4 → 4.54/4
defect_12x12_d.save   6.32/4 → 6.32/4    defect_12x12_p.save   6.34/4 → 6.34/4
defect_5x5_d.save     0.26/4 → 0.26/4    defect_5x5_p.save     0.26/4 → 0.26/4
defect_6x6_d.save     0.48/4 → 0.48/4    defect_6x6_p.save     0.48/4 → 0.48/4
defect_7x7_d.save     0.83/4 → 0.83/4    defect_7x7_p.save     0.84/4 → 0.84/4
defect_8x8_d.save     1.36/4 → 1.36/4    defect_8x8_p.save     1.37/4 → 1.37/4
defect_9x9_d.save     2.11/4 → 2.11/4    defect_9x9_p.save     2.12/4 → 2.12/4
defect_uc_dense_24    1.70/579           defect_uc_dense_25    1.85/628
defect_uc_dense_27    2.15/732           defect_uc_dense_28    2.32/787
defect_uc_dense_32    3.02/1027
defect_unit_cell_5x5  0.06/28   6x6 0.09/39   7x7 (defect_7x7.save) 0.12/52   8x8 0.15/67
defect_unit_cell_9x9  0.20/84   10x10 0.24/103   11x11 0.43/184   12x12 0.35/147
defect_unit_cell_{5x5,7x7,8x8}_wann.save  0.06/28  0.12/52  0.15/67 ; dossiers _wann complets 0.12 / 0.25 / 0.32
```

Note : le rsync initial lancé avec `--files-from` sans `-r` explicite n'a créé
que le squelette (la récursion n'est pas implicite avec `--files-from`) ; relancé
avec `-r`. Le miroir ajoute 51,8 Gio au projet.

---

## Phases 1–2 — DRY-RUN : manifestes de suppression (aucune suppression)

Manifestes dans `ab-initio-defects/cleanup/cleanup_manifest_{1..6}.txt`
(format : `chemin ⇥ taille_octets ⇥ mtime ⇥ atime ⇥ justification`, en-tête `#`
avec les totaux par justification). Listings pris le 2026-09-17 (`find`
frais, pas ceux de P12). Non commités.

| Étage | Contenu | Fichiers | Gio récupérés | Manifeste |
|---|---|---|---|---|
| 1 | wfc des `phonons/_ph0` des grilles EPW test 36k-30q (578,9), 30k-24q (262,0), 24k-12q (44,7), 16k-16q (34,7), 16k-12q (19,9), 12k-12q (11,3) | 29 824 | **951,5** | `cleanup_manifest_1.txt` |
| 2 | ère ABINIT : `graphene/abinit/pristine` (403,1), `graphene/abinit/defective` (17,6), `ab-initio-defects/jobs` (27,9) — avec tar de trace nearline (3 490 fichiers texte, 0,5 Gio) | 9 398 | **448,6** | `cleanup_manifest_2.txt` |
| 3 | wfc des runs DFPT prt terminés : `16k-8q/dfpt_g_{G,K}` (68,9 ×2), `24k-24q/dfpt_g_{G,K}` (66,5 ×2) ; **gardés et listés** : 2 874 fichiers, 11,9 Gio | 16 056 | **270,8** | `cleanup_manifest_3.txt` |
| 4 | `results/M` : `*_norm.npy` (28, 19,6), `obsolete_grid_7x7/` (13, 14,7), `*_nb16` (32, 39,5) | 73 | **73,7** | `cleanup_manifest_4.txt` |
| 5 | scratch : `qe_conv/` (9 492, 44,9), `prefix.wfcN` en vrac de `qe_tmp` (758, 15,3), 52 répertoires vides + `compmatphys/` | 10 306 | **60,2** (scratch) | `cleanup_manifest_5.txt` |
| 6 | `graphene/qe/test/` (31,5), `dft/` (1,2), `abinit_assignment/` (2,4, tar nearline), vrac §5.5 (10 entrées) | 12 829 | **35,1** | `cleanup_manifest_6.txt` |
| **Total** | | **78 486** | **1 839,9** (projet 1 779,7 + scratch 60,2) | |

Le projet passerait de ≈ 1 983 Gio (1 931 + 51,8 de miroir) à ≈ 203 Gio.

### Points à arbitrer avant chaque GO

- **Étage 1 — clause « < 100 Mo » inapplicable.** Les wfc de `_ph0` font entre
  0,3 et 71 Mo chacun (30 592 fichiers) ; aucun n'atteint 100 Mo, la clause
  viderait le manifeste. Règle appliquée à la place : nom contenant `.wfc` ET
  chemin sous `phonons/_ph0/`, rien d'autre. Gardés hors manifeste : dvscf
  (239, 3,7 Gio), dyn* (1 918), xml, hdf5 (densités), entrées/sorties, et
  `phonons/graphene.save/wfc*.hdf5` (768 fichiers, 0,9 Gio — supprimables aussi,
  à dire).
- **Étage 2 — tar de trace** : `/nearline/rrg-cotemich-ac/gregb26/abinit_era_trace_2025-11_2026-02.tar.gz`
  (dossier nearline personnel existant, vide). Contenu = entrées `.abi`, sorties
  `.abo*`, logs, `.sbatch/.py/.sh`, EIG/DDB texte, `.win/.wout/.eig/.kpt/.gnu/.dat/.agr/.txt`,
  `.psp8` : tout fichier < 10 Mo hors binaires ; exclus WFK/DEN/POT/UNK, `*.nc`,
  `*.1`, `.mat/.amn/.chk/.mmn`, `.npy`, `.xsf`. La commande est dans l'en-tête du
  manifeste. Étage 2 = supprimer les deux arbres en entier après le tar.
- **Étage 3** : `16k-8q/epw_g_{G,K}` gardés en entier (1,3 Gio chacun, dont
  1,2 Gio de wfc NSCF), conformément à la consigne. Les xml de `_ph0`
  (patterns/dynmat, 1,7–1,8 Gio par dir) sont gardés avec les dyn.
- **Étage 4 — `_nb16`** : deux options. (a) tar nearline
  `results_M_nb16_2026-09.tar` (39,5 Gio, npy incompressibles, 32 fichiers)
  puis suppression ; (b) suppression sèche, le test 16→20 étant consigné dans le
  mémoire (tab:tests_M) et `check_M_dense_nb20_vs_nb16.py` devenant
  inexécutable. Recommandation : (b) si la table du mémoire suffit comme trace,
  sinon (a) — c'est un arbitrage de provenance, pas d'espace. Les `*_norm.npy`
  se régénèrent en minutes (`scripts/migrate_M_norm.py`) ; `obsolete_grid_7x7/`
  n'est lu que par `_old_vs_new_7x7.py`.
- **Étage 5** : `qe_tmp/graphene_scf/*.save` (1,95 Gio, reproductible) n'est
  pas dans le manifeste, seuls ses `.wfcN` en vrac le sont — à ajouter si voulu.
  Les 52 répertoires vides incluent 20 `.save` vides de `qe_conv/`.
- **Étage 6** : `hamiltonian_reconstruction.png` est **suivi par git**
  (commit 2a30600) : `git mv` vers `figures/` ou `git rm`, pas un `rm`. Tar
  `abinit_assignment_gaas_2026-05.tar.gz` (95 fichiers texte, 5 Mo) avant
  suppression des 2,4 Gio.

### Garde-fous respectés

- `abinit_processing/` : intouché ; `abinit_processing/README` (une ligne)
  ajouté : `pseudo/` est le `pseudo_dir` absolu de **193** entrées QE
  (`graphene/qe/**/*.in`) et de `$PSEUDOS` dans `~/.bashrc`.
- `24k-24q/epw1/` et ses liens relatifs `../epw1/*.fmt` : hors de tout manifeste
  (seul le `.epw1.in.swp` et le log OOM y sont listés, étage 6).
- `codes/` : hors manifeste. **Question libwannier90 :** aucun module chargé
  ni variable d'environnement ne pointe vers `codes/wannier90` ; `libwannier90.a`
  (bibliothèque statique) n'est liée par rien à l'exécution. En revanche
  `~/bin/wannier90.x` (en tête du `PATH`, donc prioritaire sur le module
  `wannier90/3.1.0` de la grappe) est une **copie octet pour octet**
  (md5 `e51c0c82…`) du `codes/wannier90/wannier90.x` compilé le 2025-11-24
  (v3.1.0, gfortran), et c'est ce binaire que les jobs Wannier
  (`graphene/qe/defects/unit_cell/*/submit_*.sh`, `srun -n 1 wannier90.x`) exécutent.
  Supprimer `codes/wannier90/` ne casserait donc pas `~/bin/wannier90.x`, mais
  ferait perdre la provenance de ce binaire (`make.inc`, sources). `codes/qe` et
  `codes/abinit` (sources + `.git`, 2,2 Gio) ne sont référencés nulle part.
  Recommandation : garder `codes/wannier90/` (0,6 Gio) tant que la chaîne
  Wannier peut être relancée ; `codes/qe` et `codes/abinit` sont supprimables
  (recompilables depuis les dépôts publics), à décider dans un étage 7 éventuel.
- Scratch : aucun `touch`, aucune lecture sans `O_NOATIME` ; rien déplacé.

---

## Phase 3 — EN ATTENTE DU GO : renommage du dépôt

À faire après confirmation que le dépôt GitHub a été renommé (nom annoncé :
« graphene-raman », à confirmer dans le GO) :

```
cd ~/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
git remote set-url origin https://github.com/Gregb26/graphene-raman.git
git fetch origin && git remote -v
```

Remote actuel : `origin https://github.com/Gregb26/ab-initio-defects.git`.
Le dossier local `ab-initio-defects/` garde son nom ; aucun `sed` sur les
scripts (chemins absolus de §4 de l'inventaire inchangés).

---

## État (mis à jour 2026-09-23, session R0)

- Phase 0 : **exécutée, vérifiée** (2026-09-17 : 51,8 Gio, 5 007 fichiers, md5 OK ; 2026-09-22 : + `vacancy_relaxed/nspin{1,2}` 6,3 Gio, md5 OK, revérifié 2026-09-23).
- Phase 3 : **faite le 2026-09-23** — `origin` → `https://github.com/Gregb26/graphene-raman.git`, `git fetch` OK, `main` à jour avec `origin/main` (78ed94a). Dossier local inchangé.
- Étage 1 : **exécuté le 2026-09-23** (952,4 Gio). Étage 3 : **exécuté le 2026-09-23** (270,8 Gio). Étage 2 : **exécuté le 2026-09-23** (448,6 Gio, tar de trace sur nearline). Étage 6 : **exécuté le 2026-09-23** (35,1 Gio, tar abinit_assignment sur nearline). Étage 4 : **exécuté le 2026-09-23** (73,7 Gio, tar nb16 sur nearline). Étage 5 : **exécuté le 2026-09-24** (60,2 Gio, `vacancy_relaxed/` exclu). Séquence GO du 2026-09-23 : 3 → 2 → 6 → 4, un GO par étage après find frais + diff ; étage 5 **exécuté le 2026-09-24** après R2 (avec `vacancy_relaxed/` exclu).
  Arbitrages reçus : étage 2 tar de trace nearline + md5 + `tar -tzf` avant rm ; étage 6 `git mv` de hamiltonian_reconstruction.png vers figures/ et tar abinit_assignment ; étage 4 option (a) tar `results_M_nb16_2026-09.tar` nearline + md5 + `tar -tf`, `_norm.npy` et `obsolete_grid_7x7/` supprimés directement.
- Étage 7 (R1, scratch `qe_tmp/vacancy_relaxed/k3x3`) : **exécuté le 2026-09-23** (46 fichiers, 146,8 Gio ; voir journal).
- Commits (F) : **faits le 2026-09-23** — 44ceaa4 (CLEANUP, INVENTAIRE, CLAUDE.md, .gitignore `cleanup/`), 11d1695 (article/R1_vacancy_relaxed resynchronisé + .out pw.x), 9fbccd6 (CLAUDE.md règle 5 : répertoire de travail vs copie versionnée). Non poussés.

## Journal

### 2026-09-23 — session R0

- **A. Phase 3** : `git remote set-url origin https://github.com/Gregb26/graphene-raman.git` ; `git fetch origin` sans sortie (rien de nouveau) ; `git remote -v` montre le nouveau remote ; `git status` : `main` à jour, fichiers non suivis inchangés.
- **B. Campagne R1 déplacée** : `graphene/qe/vacancy_relaxed/` → `graphene/qe/defects/super_cell_relaxed/9x9/` (`mv`, contenu intact : nspin1/, nspin2/, k3x3/, projwfc/, scripts, rapports, analyse). Aucun `sed` nécessaire : les scripts (`make_inputs.py`, `analyze_relax.py`, `k3x3/make_k3x3.py`, `k3x3/analyze_k3x3.py`) ne référencent que `graphene/qe/defects/super_cell/9x9/defective/` (inchangé) et le scratch ; les rapports ne citent que `qe_tmp/vacancy_relaxed` (scratch, inchangé) et le miroir. Ligne d'en-tête ajoutée à `R1_rapport.md`. `outdir` scratch et miroir `qe_tmp_backup/vacancy_relaxed/` non touchés.
  Constat : une copie versionnée existait déjà dans le dépôt, `article/R1_vacancy_relaxed/` (commit 78ed94a du 2026-09-22, 27 fichiers). Elle a été resynchronisée depuis le répertoire déplacé (R1_rapport.md, README.md corrigé pour le nouveau chemin) et complétée par les `.out` de pw.x < 5 Mo (nspin{1,2}/relax.out, k3x3/*/scf.out, 6 fichiers, 2,5 Mo). Laissés hors dépôt : `projwfc/` (185 Mo : projwfc.out 22 Mo, `proj_*.projwfc_{up,down}` 75 Mo ×2 et 8 Mo ×2, pdos_*), `slurm-*.{out,err}`, `JOBID`.
- **C. CLAUDE.md** : section « Campagnes de calcul » ajoutée (règle en 4 points + tableau de classement : R1 nspin1/2 = PRODUCTION miroité ; R1b k3x3 = TEST consigné ; chaîne M = PRODUCTION ; EPW 24k-24q = PRODUCTION ; grilles test avril = TEST consigné).
- **D. Manifeste 7 (R1)** : prérequis `md5sum -c` du miroir nspin1 (4 fichiers) et nspin2 (6 fichiers) : tous OK. Scratch `qe_tmp/vacancy_relaxed` : **aucun `.wfcN` en vrac** (wfc collectées) ; manifeste = les 3 `.save` de `k3x3/` (ideal_nspin1 36,7 ; relax1_nspin1 36,7 ; relax2_nspin2 73,5 Gio ; `ideal_nspin2` n'a pas de `.save` sur le scratch). Total 46 fichiers, 146,8 Gio. STOP, GO séparé.
- **E. Étage 1 exécuté.** `find` frais restreint à `{36k-30q,30k-24q,24k-12q,16k-16q,16k-12q,12k-12q}/phonons` (33 902 fichiers) ; règle (a) `_ph0/` + `.wfc` : 29 824 fichiers, 951,52 Gio, **identique au manifeste du 17** (0 chemin en plus, 0 en moins, 0 taille différente) ; règle (b) ajoutée `phonons/graphene.save/wfc*.hdf5` : 358 fichiers, 0,85 Gio. Garde-fous : 0 chemin hors des 6 grilles, 0 lien symbolique, 0 dvscf/dyn/xml dans la liste. Suppression `xargs rm -f` (liste conservée : `cleanup/etage1_supprimes_2026-09-23.lst`) : **30 182 fichiers, 952,4 Gio récupérés, 0 fichier restant**. `du graphene/qe/epw/` : **1,6 T avant → 686 G après**. Grilles après : 36k-30q 2,2 G, 30k-24q 1,4 G, 24k-12q 516 M, 16k-16q 629 M, 16k-12q 426 M, 12k-12q 395 M (dvscf, dyn, patterns, charge-density gardés). Reste non couvert par les deux règles : `phonons/graphene.wfcN` en vrac à la racine de `phonons/` (128 par grille, ≈ 0,9 Gio en tout) — à ajouter à un étage ultérieur si voulu. Quota groupe : 40 → 39 To.
- **F. Préparé, non commité** : voir `git status` ; `cleanup/` reste hors dépôt (à mettre dans `.gitignore` au moment du commit 1 si souhaité).
- **GO étage 7 (R1) exécuté** (09:00) : `find` frais sur `qe_tmp/vacancy_relaxed/k3x3/*/*.save` = 46 fichiers, 146,79 Gio, **identique au manifeste 7** (0 écart de chemin, 0 de taille) ; garde-fous : 0 chemin hors `k3x3/*/k3x3_*.save`, 0 fichier nspin1/2. `xargs rm -f` (liste : `cleanup/etage7_supprimes_2026-09-23.lst`) puis `rmdir` des 3 `.save` vides : **46 fichiers, 146,8 Gio récupérés, 0 restant**. Scratch `qe_tmp/vacancy_relaxed` : 154 G → 6,4 G (nspin1/2 `.save` + 5 `prefix.xml`). Miroir nspin1/nspin2 revérifié après coup : md5 OK.
- Règle 5 ajoutée à CLAUDE.md (commit 9fbccd6) : répertoire de travail complet hors dépôt (`graphene/qe/defects/super_cell_relaxed/9x9/`), copie versionnée dans `article/<campagne>/` (inputs, submit, scripts, rapports, analyses, README, `.out` pw.x < 5 Mo).
- **GO étage 3 exécuté** (09:03) : `find` frais sur les 4 `dfpt_g_{G,K}` (18 332 fichiers) ; à supprimer (nom contenant `wfc`) = 16 056 fichiers, 270,79 Gio, **identique au manifeste du 17** (0 écart de chemin, 0 de taille) ; garde-fous 0/0. `xargs rm -f` (liste : `cleanup/etage3_supprimes_2026-09-23.lst`) : **16 056 fichiers, 270,8 Gio récupérés, 0 restant**. `du` avant → après : 16k-8q/dfpt_g_G 72 G → 2,4 G, dfpt_g_K 72 G → 2,4 G, 24k-24q/dfpt_g_G 69 G → 2,4 G, dfpt_g_K 69 G → 2,4 G. Gardés : 2 276 fichiers (ph.out, nscf.out, scf.out, 124 dvscf, dyn, xml, densités, logs) ; `16k-8q/epw_g_{G,K}` non touchés. `graphene/qe/epw/` : 686 G → 415 G.
- **GO étage 2 exécuté** (09:07–09:09) : `find` frais sur `graphene/abinit` + `ab-initio-defects/jobs` = 9 398 fichiers, 448,55 Gio, **identique au manifeste du 17** (0 écart de chemin, de taille ou de classe). Tar de trace : `/nearline/rrg-cotemich-ac/gregb26/abinit_era_trace_2025-11_2026-02.tar.gz` (3 490 fichiers, 0,51 Gio → 131 Mo compressé, chemins relatifs à `gregb26/`), md5 `9d5566fc5fb745bb6ff0c31b1e25e93d` (fichier `.md5` à côté), `tar -tzf` = 3 490 entrées identiques à la liste, extraction réelle de 3 fichiers échantillon comparée `cmp` OK. Puis `rm` des 9 398 fichiers (liste : `cleanup/etage2_supprimes_2026-09-23.lst`, 0 fichier hors liste dans les deux arbres) et `rm -rf` des deux arbres : **448,6 Gio récupérés**. `du` avant → après : `graphene/abinit` 421 G → supprimé, `jobs/` 28 G → supprimé ; `graphene/` 526 G, `ab-initio-defects/` 133 G. Quota groupe 39 To, fichiers 4 710 k → 4 682 k.
- **GO étage 6 exécuté** (09:11–09:13) : `find` frais sur `graphene/qe/test`, `dft`, `abinit_assignment` = 12 819 fichiers, 35,06 Gio, **identique au manifeste du 17** (0 écart). Tar `/nearline/rrg-cotemich-ac/gregb26/abinit_assignment_gaas_2026-05.tar.gz` (95 fichiers texte, 1,05 Mo, md5 `d5c3415e2afdb2634feed566aeeea25f`, `tar -tzf` = 95 entrées, extraction test OK). `rm` des 12 819 fichiers (liste : `cleanup/etage6_supprimes_2026-09-23.lst`) puis `rm -rf` des trois arbres : **35,1 Gio récupérés** (`graphene/qe/test` 32 G, `dft` 1,2 G, `abinit_assignment` 2,5 G → supprimés). Vrac : `.epw1.in.swp`, `epw1_oom_20637122.out`, lien cassé `gaas.win`, `.pytest_cache/` supprimés ; `hamiltonian_reconstruction.png` → `figures/` par `git mv` (indexé, commit avec ce journal). `graphene/` 526 G → 495 G.
- **GO étage 4 exécuté** (09:18–09:21) : `find` frais sur `results/M` = 73 fichiers, 73,72 Gio, **identique au manifeste du 17** (0 écart). Option (a) : tar `/nearline/rrg-cotemich-ac/gregb26/results_M_nb16_2026-09.tar` (32 fichiers `*_nb16.{npy,json}`, 39,45 Gio, non compressé, 2 min 23 s), md5 `82336114a3cb045d483c55361a7f2ff6` (fichier `.md5` à côté), `tar -tf` = 32 entrées identiques à la liste, extraction test d'un json et de `M_dense_5x5_nb16.npy` comparés `cmp` OK. Puis `rm` des 73 fichiers (liste : `cleanup/etage4_supprimes_2026-09-23.lst`) et `rmdir obsolete_grid_7x7/` : **73,7 Gio récupérés** (nb16 39,5 ; `_norm` 19,6 ; obsolète 14,7). `results/M` 133 G → 59 G (130 entrées). Attention : `compute_spectral.py`, `compute_tmatrix.py`, `_eta_scan.py`, `_mcheck.py` lisent `M_ed_{N}_norm.npy` (chaîne grossière, hors production) — relancer `scripts/migrate_M_norm.py` avant de les réutiliser ; `check_M_dense_nb20_vs_nb16.py` et `_old_vs_new_7x7.py` ne sont plus exécutables (données sur nearline / supprimées).

### Bilan des suppressions du 2026-09-23

| Étage | Fichiers | Gio | Archive nearline |
|---|---|---|---|
| 1 | 30 182 | 952,4 | — |
| 3 | 16 056 | 270,8 | — |
| 2 | 9 398 | 448,6 | `abinit_era_trace_2025-11_2026-02.tar.gz` (131 Mo) |
| 6 | 12 819 + 7 vrac | 35,1 | `abinit_assignment_gaas_2026-05.tar.gz` (1 Mo) |
| 4 | 73 | 73,7 | `results_M_nb16_2026-09.tar` (39,5 Gio) |
| 7 (scratch) | 46 | 146,8 | — |
| 5 (scratch, 2026-09-24) | 10 306 | 60,2 | — |
| **Total** | **78 880** | **1 987,6** (projet 1 780,6 + scratch 207,0) | 39,6 Gio sur nearline |

Projet gregb26 après : `graphene/` 495 G (dont `qe/epw` 415 G, `qe_tmp_backup` 58 G, `qe/defects` ≈ 22 G), `ab-initio-defects/` ≈ 60 G (`results/M` 59 G), `codes/` 2,8 G ; ≈ 560 Gio au total contre 1 983 avant.

### 2026-09-24 — R2 phase 2 : miroir des `.save` de la série (ajout, aucune suppression)

`graphene/qe/qe_tmp_backup/vacancy_relaxed/series/` : 84 fichiers, 54,45 Go (51 Gio), `MD5SUMS_series_2026-09-24.txt` vérifié 84/84 contre le scratch (O_NOATIME).
`qe_tmp_backup` passe de 58 G à ≈ 109 G. Piège rencontré : `rsync -a` du scratch recrée les dossiers sans setgid → fichiers au groupe `gregb26` ; corrigé (`chgrp -R` + `chmod g+s`), désormais `--chmod=Dg+s`.
Étage 5 (scratch `qe_conv`, wfc vrac, dossiers vides) toujours à faire, avec `vacancy_relaxed/` (R1 + `series/`) exclu.

### 2026-09-24 — Étage 5 EXÉCUTÉ (scratch, GO de Greg)

Manifeste frais `cleanup/cleanup_manifest_5_2026-09-24.txt` (find du jour, `vacancy_relaxed/` exclu : R1 + `series/` R2) : 10 306 entrées, 60,21 Gio,
diff contre le manifeste du 17 = 0 chemin en plus, 0 en moins, tailles de fichiers identiques (seuls les 52 répertoires vides passent de 0 à 25 600 octets, bloc Lustre).
Les 4 `compmatphys_vide` sont des journaux Wannier90 non vides (7–109 Ko, nom avec U+00A0), repris par nom explicite.
Supprimé : 10 254 fichiers (60,21 Gio : `qe_conv/` 9 492 = 44,86 Gio ; `prefix.wfcN` en vrac de 11 dossiers `qe_tmp` 758 = 15,34 Gio ; 4 journaux),
52 répertoires vides listés, puis 236 répertoires de `qe_conv/` vidés par l'étape 1 et `qe_conv/` lui-même. Liste : `cleanup/etage5_supprimes_2026-09-24.lst` (10 542 lignes). 0 ignoré.
Non touché : `qe_tmp/graphene_scf/graphene_scf.save` (1,95 Gio, hors manifeste), `qe_tmp/vacancy_relaxed/` (k3x3, nspin1, nspin2, series).
Après : scratch 113 Go / 5 751 fichiers (avant 166 Go / 16 000) ; `qe_tmp` 121 G → 109 G ; `qe_conv` 38 G → supprimé.
Restent 9 répertoires devenus vides par l'étape 2 (parents des vides listés, hors manifeste, non supprimés) :
  /scratch/gregb26/compmatphys/basic/Al
  /scratch/gregb26/compmatphys/basic/AlFe
  /scratch/gregb26/compmatphys/basic/Si/wannier/test
  /scratch/gregb26/graphene/supercell/centered
  /scratch/gregb26/graphene/unitcell
  /scratch/gregb26/jobs/10x10x1_5x5_Sc/5wann
  /scratch/gregb26/jobs/10x10x1_5x5_Sc/8wann
  /scratch/gregb26/jobs/9x9
  /scratch/gregb26/qe_tmp/graphene_scf/_ph0

Bilan cumulé des étages 1–7 : 78 880 fichiers, 1 987,6 Gio (projet 1 780,6 + scratch 207,0). Tous les étages du plan P13 sont exécutés.
