# R7 — Super-cellules 15×15 et 18×18, lacune non relaxée, sans spin (DFT seulement) ; famille 3m à cinq points (article)

- But : prolonger la famille 3m (6, 9, 12) pour la position de l'état π quasi-lié et du doublet σ en fonction de la taille,
  vers la limite diluée (ajustements 1/N et 1/N² à cinq points). Pas de M ni de matrice T dans cette campagne.
- Prompt d'origine : R7 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) scf ×4 → pp.x ×4 → D1 → rapport → STOP.
- Statut : **PRODUCTION** (règle 3 de CLAUDE.md : miroir `qe_tmp_backup/` avec md5 en fin de campagne ; rien d'unique sur le
  scratch ; aucune suppression dans cette campagne, tout nettoyage passe par un manifeste et un GO séparé).
- Dates : phase 0 le 2026-09-25 (inputs écrits) ; GO reçu le 2026-09-25 (ressources A ; pp.x chaînés en afterok) ; jobs soumis vers 15 h 30,
  tous COMPLETED à 17 h 51 (scf, pp.x, D1 cinq points) ; rapport et copie article le 2026-09-25 (STOP : miroir des wfc, GO 21/24/27, commit).
- Emplacements (règle 5 de CLAUDE.md) : les calculs QE vivent comme les tailles 5…12 dans
  `graphene/qe/defects/super_cell/{15x15,18x18}/{defective,pristine}/` (`scf.in`, `submit.scf`, `pp.in`, `submit.pp`, puis
  `scf.out`, `pp.out`, `Vks_NxN_{d,p}`) avec `outdir` sur le scratch (`qe_tmp/defect_NxN_{d,p}/`) ; ce répertoire-ci
  (`graphene/qe/defects/R7_tailles_3m/`, à côté de R4/R5 qu'il prolonge) porte le générateur d'inputs, le pilote D1, le rapport
  et les résultats (`d1/`, `fig/`). Copie versionnée prévue : `graphene-raman/article/R7_tailles_3m/` (inputs, submit, `scf.out`
  < 5 Mo, générateur, pilote, rapport, tables, figures, extrait npz de la fenêtre ; jamais `.save`, `Vks_*`, `pp.out`, `slurm-*`, `JOBID`).
- Contenu : `R7_rapport.md` (phase 0 : inputs, règle de position, règle nbnd, ressources, disque, manifeste proposé, plan ; phase GO : jobs, D1, résultats),
  `r7_driver.py` + `submit_r7.sh` (D1 : protocole R4/R5 route quadruplet, porte de régression sur R5 C, ajustements 1/N et 1/N², extraits npz, tables, figures),
  `make_inputs_r7.py` (dérive 15×15/18×18 du 12×12 ; `--check` régénère 9×9 et 12×12 byte à byte), `inputs_diff_vs_12x12.txt`
  (diffs des huit fichiers d'entrée, positions exclues), `phase0_inventory.py` + `phase0_analysis.txt` (inventaire de phase 0 : règle de position, nbnd, XML, temps, mémoire, tailles).
- Inputs dérivés du 12×12 : seuls nat, CELL_PARAMETERS, ATOMIC_POSITIONS, nbnd, prefix, outdir et ressources changent.
  Lacune sur le sous-réseau A, atome A le plus proche du centre (i = j = ⌊N/2⌋ : 15×15 atome 225, 18×18 atome 343 de la parfaite).
- Extension 21×21 / 24×24 / 27×27 (demande de Greg, 2026-09-25) : inputs écrits dans `super_cell/{21x21,24x24,27x27}/`, **non soumis, GO
  attendu** ; deux règles adaptées (rapport §0b) : nbnd = N_occ + ⌈30 (N/12)²⌉ (+92/+120/+152, couverture ≈ 1,96 eV comme la 12×12) et
  nœuds entiers à 192 rangs (2 × 96, 3 × 64, 4 × 48 ; `--mem=0` ; 6/8/12 h) ; lacune A atomes 441 / 601 / 729 ; `inputs_diff_21_24_27.txt`.
- Jobs (fichier `JOBID`) : scf 15×15 p 21818658, d 21818685 ; 18×18 p 21818687, d 21818689 ; pp.x chaînés afterok 21818684 / 21818686 /
  21818688 / 21818690 ; D1 régression 6/9/12 21818859 (`d1_reg/`) ; D1 complet (6, 9, 12, 15, 18) à soumettre en afterok des quatre pp.x
  une fois la régression PASS. Surveillance sans relance. Pilote `r7_driver.py` (d1, tables) + `submit_r7.sh` (16 cœurs, 64 G, 1 h).
