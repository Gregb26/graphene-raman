# R7 — Super-cellules 15×15 et 18×18, lacune non relaxée, sans spin (DFT seulement) ; famille 3m à cinq points (article)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 4 figures hors mémoire (`size_3m`, `localized_3m`, `*_d1_reg`), `make_inputs_r7.py`, `make_inputs_r7_relax.py`, `phase0_inventory.py`, `submit_r7.sh`, `submit_mirror_wfc*.sh`. `r7_driver.py` est gardé : sa sous-commande `tables` redessine `fig/size_3m_d1_8pts` et `fig/localized_3m_d1_8pts` (mémoire) depuis `d1_8pts/d1_results.json` ; il importe `../R4_quasi_lie/r4_driver.py` et `../R5_base_vs_M/r5_*.py`, gardés pour cette seule raison. Gardés : inputs et `scf.out` des super-cellules 15…27, `relax.in`/`relax.out` de R7c, rapport, tables, json/npz de `d1/ d1_reg/ d1_8pts/`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R7_tailles_3m/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R7 — Famille 3m de 6×6 à 27×27 : état π quasi-lié et doublet σ en fonction de la taille (+ R7c)
- **Question** : vers quoi tendent l'état π quasi-lié et le doublet σ de la lacune (non relaxée, sans spin, DFT seule) quand on prolonge la famille 3m (6, 9, 12) à 15, 18, puis 21, 24 et 27 ? Quelle loi, 1/N ou 1/N², décrit la limite diluée ? R7c : relaxer en nspin 1 les lacunes de 15 à 27.
- **Pourquoi** : la référence R5 n'avait que trois points (6, 9, 12) et donnait des ε_∞ très différents selon la loi (π : −0,046 en 1/N, −0,409 en 1/N² ; σ : +0,055 / +0,087). Pas de M ni de matrice T dans cette campagne.
- **Méthode** : inputs dérivés de la 12×12 (régression `--check` byte à byte PASS ; lacune sur le sous-réseau A, la plus proche du centre ; nbnd = N_occ + 30 pour 15 et 18, N_occ + ⌈30 (N/12)²⌉ pour 21, 24 et 27). scf pw.x, puis pp.x (V_KS), parfaite et lacune. D1 selon le protocole R4/R5, route quadruplet : E_D, alignement Lu 1,0 Å, parité miroir, w₂ sur un disque de 2 Å, seuil 3⟨w₂⟩_P, fenêtre [−3, +1] eV ; porte de régression 6/9/12 contre R5 C. Ajustements ε_∞ + a/N^p. Jobs : scf 21818658/685/687/689, pp.x 21818684–690, D1 à 5 points 21819066 ; GO 2 : scf 21833299–309, pp.x 21833300–310, D1 à 8 points 21833311 ; miroirs wfc 21833312 et 21850526 ; R7c : relax 21852173–21852177, reprise de la 24×24 21908867.
- **Résultat** :
  - π quasi-lié (ε − E_D) : −1.065, −0.737, −0.551, −0.459, −0.388, −0.344, −0.307, −0.279 eV pour N = 6 → 27 ; w₂ 0,284 → 0,166.
  - Sur huit points, loi 1/N : ε_∞ = −0,0523 eV, a = −6,089 eV, rms 4,5 meV (résidu max 8,9 meV). Loi 1/N² : ε_∞ −0,2963 eV, rms 46,5 meV. Le terme 1/N² ajouté en supplément est nul (b = +0,06 eV, rms inchangé).
  - Doublet σ : dégénéré à < 0,1 meV, w₂ 0,711–0,713 quel que soit N. Il vaut +0,156 (6), +0,101 (9), +0,114 (12), puis reste sur un plateau à +0,102–0,107 eV de 15 à 27 (+0,1024 à 27). ε_∞ de +0,084 à +0,122 eV selon la loi (rms 6–10 meV).
  - Porte 6/9/12 contre R5 C : PASS, écart 0 sur 42 grandeurs. E_D −4,23886 / −4,23888 / −4,23890 eV pour 21 / 24 / 27 ; E_F(p) − E_D = +19,0 meV à toutes les tailles.
  - R7c : la 15×15 relaxe en 28 pas BFGS (ΔE_relax −365,62 meV, paire 2,14997 Å), la 18×18 en 33 pas (−362,97 meV, paire 2,14773 Å) ; référence 12×12 R2 : −367,45 meV, 2,15980 Å. 21×21 et 24×24 : non consigné (en cours à la fin du rapport). 27×27 : en hold.
- **Interprétation / décision** : pour le π, la loi 1/N l'emporte nettement : à cinq points, son rms est 6,5 fois plus petit que celui de 1/N², et le coefficient a est stable depuis trois points (−6,14 → −6,09). ε_∞ = −0,052 eV. Pour le σ, les points ne départagent pas les lois (les ajustements sont pilotés par 6 et 9), mais la position est sur un plateau depuis N = 15. Tous les `.save` de 15 à 27 sont miroités dans `qe_tmp_backup/` (md5 OK). R7c : le départ transplanté depuis la 12×12 n'apporte qu'un gain marginal (28 pas contre 31), la relaxation étant collective. La 27×27 attend le GO de Greg, avec transplant depuis la 24×24 relaxée. pp.x + D1 sur les cellules relaxées seulement sur demande.
- **Où sont les données** :
  - Campagne : `R7_rapport.md` ; `d1_8pts/D1_tables.md`, `d1_8pts/d1_results.json`, `d1_8pts/window_{6x6,…,27x27}.npz` ; `d1/` (5 points) et `d1_reg/` (régression), même format ; `fig/size_3m_d1_8pts.{pdf,png}`, `fig/localized_3m_d1_8pts.{pdf,png}` (figures du mémoire, à refaire au style `memoire.mplstyle` ; les variantes à cinq points et de régression ont été retirées) ; `super_cell/{15x15,18x18,21x21,24x24,27x27}/{defective,pristine}/{scf.in,scf.out,pp.in,submit.scf,submit.pp}` ; `super_cell_relaxed/series/NxN/nspin1/relax.in` (+ `relax.out` pour 15×15 et 18×18), `super_cell_relaxed/series/NxN/transplant.log` ; `relax_phase0.txt`, `phase0_analysis.txt`, `inputs_diff_vs_12x12.txt`, `inputs_diff_21_24_27.txt`, `r7_log.txt`, `JOBID`.
  - Hors dépôt : les `.save` sont sur le scratch `qe_tmp/defect_NxN_{d,p}/` et miroités dans `qe_tmp_backup/defect_NxN_{d,p}/` (`MD5SUMS_R7_2026-09-25.txt`, `MD5SUMS_R7_2026-09-26.txt`). Les `Vks_NxN_{d,p}` sont dans `graphene/qe/defects/super_cell/`.

---

- But : prolonger la famille 3m (6, 9, 12) pour la position de l'état π quasi-lié et du doublet σ en fonction de la taille,
  vers la limite diluée (ajustements 1/N et 1/N² à cinq points). Pas de M ni de matrice T dans cette campagne.
- Prompt d'origine : R7 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) scf ×4 → pp.x ×4 → D1 → rapport → STOP.
- Statut : **PRODUCTION** (règle 3 de CLAUDE.md : miroir `qe_tmp_backup/` avec md5 en fin de campagne ; rien d'unique sur le
  scratch ; aucune suppression dans cette campagne, tout nettoyage passe par un manifeste et un GO séparé).
- Dates : phase 0 le 2026-09-25 (inputs écrits) ; GO reçu le 2026-09-25 (ressources A ; pp.x chaînés en afterok) ; jobs soumis vers 15 h 30,
  tous COMPLETED à 17 h 51 (scf, pp.x, D1 cinq points) ; GO 2 ≈ 18 h (miroir wfc 15/18 fait, md5 4/4) ; 21/24/27 : scf, pp.x et D1 huit points
  COMPLETED dans la nuit du 25 au 26 (dernier job 05 h 49, `d1_8pts/`) ; rapport complété le 2026-09-26 ; pas de commit (consigne).
- Emplacements (règle 5 de CLAUDE.md) : les calculs QE vivent comme les tailles 5…12 dans
  `graphene/qe/defects/super_cell/{15x15,18x18}/{defective,pristine}/` (`scf.in`, `submit.scf`, `pp.in`, `submit.pp`, puis
  `scf.out`, `pp.out`, `Vks_NxN_{d,p}`) avec `outdir` sur le scratch (`qe_tmp/defect_NxN_{d,p}/`) ; ce répertoire-ci
  (`graphene/qe/defects/R7_tailles_3m/`, à côté de R4/R5 qu'il prolonge) porte le générateur d'inputs, le pilote D1, le rapport
  et les résultats (`d1/`, `fig/`). Copie versionnée prévue : `graphene-raman/campagnes/R/R7_tailles_3m/` (inputs, submit, `scf.out`
  < 5 Mo, générateur, pilote, rapport, tables, figures, extrait npz de la fenêtre ; jamais `.save`, `Vks_*`, `pp.out`, `slurm-*`, `JOBID`).
- Contenu : `R7_rapport.md` (phase 0 : inputs, règle de position, règle nbnd, ressources, disque, manifeste proposé, plan ; phase GO : jobs, D1, résultats),
  `r7_driver.py` + `submit_r7.sh` (D1 : protocole R4/R5 route quadruplet, porte de régression sur R5 C, ajustements 1/N et 1/N², extraits npz, tables, figures),
  `make_inputs_r7.py` (dérive 15×15/18×18 du 12×12 ; `--check` régénère 9×9 et 12×12 byte à byte), `inputs_diff_vs_12x12.txt`
  (diffs des huit fichiers d'entrée, positions exclues), `phase0_inventory.py` + `phase0_analysis.txt` (inventaire de phase 0 : règle de position, nbnd, XML, temps, mémoire, tailles).
- Inputs dérivés du 12×12 : seuls nat, CELL_PARAMETERS, ATOMIC_POSITIONS, nbnd, prefix, outdir et ressources changent.
  Lacune sur le sous-réseau A, atome A le plus proche du centre (i = j = ⌊N/2⌋ : 15×15 atome 225, 18×18 atome 343 de la parfaite).
- Extension 21×21 / 24×24 / 27×27 (demande de Greg, 2026-09-25) : inputs écrits dans `super_cell/{21x21,24x24,27x27}/`, **GO reçu ≈ 18 h,
  soumis** (scf 21833299/301/303/305/307/309, pp.x afterok 21833300–310, D1 huit points 21833311 → `d1_8pts/`, miroir wfc 15/18 21833312) ; deux règles adaptées (rapport §0b) : nbnd = N_occ + ⌈30 (N/12)²⌉ (+92/+120/+152, couverture ≈ 1,96 eV comme la 12×12) et
  nœuds entiers à 192 rangs (2 × 96, 3 × 64, 6 × 32 ; `--mem=0` ; 6/8/12 h) ; lacune A atomes 441 / 601 / 729 ; `inputs_diff_21_24_27.txt`.
- Jobs (fichier `JOBID`) : scf 15×15 p 21818658, d 21818685 ; 18×18 p 21818687, d 21818689 ; pp.x chaînés afterok 21818684 / 21818686 /
  21818688 / 21818690 ; D1 régression 6/9/12 21818859 (`d1_reg/`) ; D1 complet (6, 9, 12, 15, 18) à soumettre en afterok des quatre pp.x
  une fois la régression PASS. Surveillance sans relance. Pilote `r7_driver.py` (d1, tables) + `submit_r7.sh` (16 cœurs, 64 G, 1 h).
- R7c (demande de Greg, 2026-09-26) : relaxation BFGS nspin 1 des lacunes 15/18/21/24/27 depuis la géométrie relaxée de la 12×12 (R2 nspin1)
  transplantée par inversion (B → A, R_cut 14,5 Å) ; inputs écrits dans `super_cell_relaxed/series/{15x15,…,27x27}/nspin1/` (`make_inputs_r7_relax.py`,
  régressions PASS, `relax_phase0.txt`, `transplant.log` par taille) ; GO reçu le 2026-09-26 (« GO en parallèle »), cinq relax soumises (JOBID) ; rapport §R7c.
