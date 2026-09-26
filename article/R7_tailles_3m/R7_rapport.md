# R7 — Super-cellules 15×15 et 18×18, lacune non relaxée, sans spin (DFT seulement) : rapport de campagne

Prompt R7 (Greg, 2026-09-25). Statut PRODUCTION. Ordre : phase 0 → STOP → (GO) scf ×4 → pp.x ×4 → D1 → rapport → STOP.
But : prolonger la famille 3m (6, 9, 12) pour la position de l'état π quasi-lié et du doublet σ en fonction de la taille, vers la
limite diluée. Pas de M ni de matrice T. Aucune suppression.

## Phase 0 — Inventaire et préparation (2026-09-25 ; rien soumis)

### 0.1 Inputs des tailles existantes et construction des positions

**Inputs 5…12** (`graphene/qe/defects/super_cell/NxN/{defective,pristine}/`) : `scf.in`, `submit.scf`, `pp.in`, `submit.pp`, puis
`scf.out`, `pp.out`, `Vks_NxN_{d,p}` ; `outdir` sur le scratch `qe_tmp/defect_NxN_{d,p}/` (`.save` + `prefix.xml`, aucun `.wfcN`
résiduel). Le diff 9×9 ↔ 12×12 de `scf.in` ne touche que les deux lignes de commentaire, `prefix`, `outdir`, `nat`, `nbnd`, les deux
premières lignes de `CELL_PARAMETERS` et les positions ; tout le reste est commun : `ecutwfc 100` Ry, `occupations 'smearing'`,
`mv`, `degauss 0.01`, `assume_isolated '2D'`, `electron_maxstep 400`, `conv_thr 1e-10`, `mixing_beta 0.3`, `C 12.011 C.upf`,
`pseudo_dir .../abinit_processing/pseudo`, `K_POINTS gamma`, c = 30 bohr. `submit.scf` : 64 tâches × 1 cœur, 4 G/cœur, sortie SLURM
`/dev/null`, `module restore qe`, `srun pw.x < scf.in > scf.out` ; seule la limite de temps varie (1 j pour la 9×9, 4 j pour la
12×12). `pp.in` : `plot_num = 1`, `filplot = 'Vks_NxN_{d,p}'`, `&plot iflag 3, output_format 6` sans `fileout` → le cube part sur
la sortie standard (`pp.out` : 184 Mo pour la 9×9, 328 Mo pour la 12×12). `submit.pp` : 1 tâche, 32 G, 1 h.

**Script de construction** : aucun dans le projet (`abinit_processing/`, `codes/`, `graphene/qe/convergence`, `unit_cell`,
`graphene-raman/{scripts,notebooks,src}` : rien). Les positions des huit tailles sont reproduites exactement (écart ≤ 2e-9, format
`%.10f`) par la formule A = ((i + 1/3)/N, (j + 1/3)/N), B = ((i + 2/3)/N, (j + 2/3)/N), i extérieur, j intérieur, A puis B ; la
cellule est la 12×12 × N/12 (la 9×9 existante en est exactement les 3/4) : a = 4,6597846766 bohr (2,465852 Å), γ = 60°. Le
générateur `make_inputs_r7.py` code cette règle ; `--check` régénère la 9×9 (A, i = j = 4) et la 12×12 (B, i = j = 5) et exige
l'identité byte à byte avec les `scf.in` et `pp.in` existants : **PASS**.

**Règle de position de la lacune** (fichier défectueux = fichier parfait moins une ligne, ordre conservé, pour les huit tailles) :

| N | nat p/d | atome retiré (1-based) | (i, j) | sous-réseau | position réduite | écart au centre (0,5, 0,5) |
|---|---|---|---|---|---|---|
| 5 | 50/49 | 25 | (2, 2) | A | (0.466667, 0.466667) | −0.0333 |
| 6 | 72/71 | 30 | (2, 2) | B | (0.444444, 0.444444) | −0.0556 |
| 7 | 98/97 | 49 | (3, 3) | A | (0.476190, 0.476190) | −0.0238 |
| 8 | 128/127 | 73 | (4, 4) | A | (0.541667, 0.541667) | +0.0417 |
| 9 | 162/161 | 81 | (4, 4) | A | (0.481481, 0.481481) | −0.0185 |
| 10 | 200/199 | 90 | (4, 4) | B | (0.466667, 0.466667) | −0.0333 |
| 11 | 242/241 | 121 | (5, 5) | A | (0.484848, 0.484848) | −0.0152 |
| 12 | 288/287 | 132 | (5, 5) | B | (0.472222, 0.472222) | −0.0278 |

Règle : l'atome du sous-réseau choisi le plus proche du centre de la cellule ; sous-réseau A pour 5, 7, 8, 9, 11 (i = j = ⌊N/2⌋),
B pour 6, 10, 12 (i = j = N/2 − 1). Pour R7 (A, prompt) : **15×15 → i = j = 7, atome 225, (0.4888888889, 0.4888888889)**
(écart −0.0111) ; **18×18 → i = j = 9, atome 343, (0.5185185185, 0.5185185185)** (écart +0.0185, du même côté que la 8×8 ;
sans conséquence, la cellule est périodique). Vérification indépendante des fichiers écrits (détection par image minimale) :
atome 225 / 343, sous-réseau A, trois premiers voisins à 1,4237 Å, nat 450/449 et 648/647.

**Inputs R7 écrits** (`super_cell/15x15/…`, `super_cell/18x18/…`, diffs complets dans `inputs_diff_vs_12x12.txt`) : seuls changent
les commentaires (« 15x15 »), `prefix`, `outdir`, `nat`, `nbnd`, la cellule (60.5323785915 / 34.9483850745 et 72.6388543098 /
41.9380620894 bohr), les positions, et dans `submit.scf` les ressources (§0.3) ; `pp.in` : `prefix`, `outdir`, `filplot` ;
`submit.pp` identique. Versions : `module restore qe` = PWSCF 7.5 (runs de juin 2026 et relaxations R2 de septembre) ;
`pseudo_dir` existe, `C.upf` md5 `34a24e64…` identique à la copie des `.save` de production. Les `outdir` 15×15/18×18 n'existent
pas encore sur le scratch.

### 0.2 Règle de nbnd

nbnd = N_occ + 30 pour les tailles 6…12, parfaite et lacune (5×5 : +30 parfaite, +32 lacune). Appliquée : **15×15 : 930 (p) /
928 (d) ; 18×18 : 1326 (p) / 1324 (d)**. Couverture de la fenêtre D1 [−3, +1] eV lue dans les XML de production :

| cas | nbnd | E_F (eV) | E_D quadruplet (étalement) | ε(nbnd) − E_D | états ≤ E_D + 1 eV | bandes de marge au-dessus de E_D + 1 |
|---|---|---|---|---|---|---|
| 6×6 p / d | 174 / 172 | −4.2185 / −4.2390 | −4.2375 (0.0000) | +4.92 / +4.87 | 146 / 146 | 28 / 26 |
| 9×9 p / d | 354 / 352 | −4.2195 / −4.2219 | −4.2385 (0.0000) | +3.33 / +3.30 | 326 / 326 | 28 / 26 |
| 12×12 p / d | 606 / 604 | −4.2197 / −4.2172 | −4.2387 (0.0000) | +1.96 / +1.95 | 578 / 578 | 28 / 26 |

Le nombre d'états dans [E_D, E_D + 1 eV] croît en N² (2 à 4 états vides à la 12×12), donc la dernière bande attendue reste à
≈ +1,6 eV (15×15) et ≈ +1,3 eV (18×18) au-dessus de E_D (loi n ∝ N²E², étalonnée sur la 12×12) : la fenêtre est couverte avec
≈ 20–25 bandes de marge. La règle est gardée telle quelle (même famille, même convergence Davidson des bandes vides sans
`diago_full_acc`, comme 5…12) ; l'alternative +30 × (N/12)² (47 / 68) coûterait ∝ nbnd sans rien ajouter dans la fenêtre.

### 0.3 Ressources

Murs scf 5…12 (64 rangs, `scf.out`) et mémoire :

| N | mur d (s) / it | mur p (s) / it | RAM QE estimée par rang | FFT | G denses |
|---|---|---|---|---|---|
| 5 | 12.9 / 22 | 12.9 / 16 | 43 MB | 150×150×192 | 0.95 M |
| 6 | 29.4 / 23 | 32.9 / 26 | 76 MB | 180×180×192 | 1.37 M |
| 7 | 50.4 / 20 | 41.0 / 16 | 131 MB | 216×216×192 | 1.87 M |
| 8 | 118.2 / 24 | 98.9 / 17 | 209 MB | 240×240×192 | 2.44 M |
| 9 | 199.8 / 24 | 160.9 / 19 | 321 MB | 270×270×192 | 3.09 M |
| 10 | 300.8 / 23 | 331.9 / 18 | 475 MB | 300×300×192 | 3.81 M |
| 11 | 721.0 / 26 | 472.7 / 18 | 684 MB | 360×360×192 | 4.61 M |
| 12 | 907.1 / 27 | 673.3 / 17 | 946 MB | 360×360×192 | 5.49 M |

Ajustements t ∝ N^p : lacune p = 4,96 (N = 5…12) / 5,28 (N = 8…12) → **15×15 : 44–51 min, 18×18 : 110–134 min** ; parfaite
p = 4,59 / 4,87 → 31–35 / 71–85 min ; par itération (× 27 itérations comme la 12×12) : 45–49 / 106–122 min. Le nombre
d'itérations peut croître avec la longueur de cellule (12×12 : 27 pour la lacune, 17 pour la parfaite) ; 40 itérations à la 18×18
feraient ≈ 3 h.

Mémoire par rang : estimation QE (borne inférieure) × (N/12)⁴ → 2,3 G (15×15), 4,7 G (18×18) ; sacct de juin : 12×12 AveRSS
1,24–1,27 G par rang, **MaxRSS d'un rang 5,2 G (parfaite) et 7,3 G (lacune)** (11×11 : 5,0–5,3 G) — le rang 0 porte la base plus
la fonction d'onde complète (6,6 G) au moment de l'écriture de `wfc1.hdf5`. Extrapolation × (N/12)⁴ : moyenne 3,1 / 6,4 G, rang
maximal ≈ 18 / 37 G (= base + wfc 16 / 34 G). En juin les 64 rangs étaient répartis sur 5–10 nœuds à 4 G/cœur (cgroup ≈ 28 G par
nœud, suffisant pour 7,3 G) ; ce n'est plus suffisant pour 18–37 G sur un rang.

**Proposition A (écrite dans les `submit.scf`)** : 64 rangs, `--ntasks-per-node=8` (8 nœuds), `--mem-per-cpu=8G` pour la 15×15
(64 G par nœud ≥ 19 + 7 × 3,1 = 41 G) et `16G` pour la 18×18 (128 G par nœud ≥ 40 + 7 × 6,4 = 85 G), `--time=03:00:00` /
`06:00:00` (≥ 2,2 × l'ajustement 8…12 ; partitions bycore b1/b2, 314 nœuds). Facturation ≈ 128 / 256 cœurs-équivalents
(4 G par cœur sur les nœuds 192 cœurs / 768 G). Sortie SLURM `/dev/null` conservée comme 5…12 (suivi par `scf.out` et sacct,
où un OUT_OF_MEMORY est visible).
**Option B (non écrite)** : un nœud exclusif, 192 rangs, `--mem=0` : même coût que A pour la 18×18, probablement 2–3 × plus rapide
si la parallélisation R&G tient à 192 rangs sur 12 M de vecteurs G (non mesuré dans cette famille) ; A garde des murs comparables
à 5…12. À décider au GO.

pp.x (`plot_num = 1`, 1 tâche) : 12×12 en 1 min 25, MaxRSS 5,3 G sous 32 G → 15×15 ≈ 8 G, 2–3 min ; 18×18 ≈ 12 G, 4–5 min :
`submit.pp` inchangé (32 G, 1 h).

### 0.4 Espace disque

Tailles (12×12 mesurées ; `wfc1.hdf5` = npw × nbnd × 16 o, npw ∝ N², nbnd ∝ N²) :

| fichier | 12×12 | 15×15 | 18×18 |
|---|---|---|---|
| `wfc1.hdf5` (par cellule) | 6.64 Go (685 905 × 604) | 16.2 Go | 33.6 Go |
| `charge-density.hdf5` | 154 Mo | 240 Mo | 346 Mo |
| `data-file-schema.xml` + `C.upf` | 0.3 Mo | ≈ 0.5 Mo | ≈ 0.7 Mo |
| `Vks_NxN_{d,p}` (filplot, /project) | 428 Mo | 669 Mo | 963 Mo |
| `pp.out` (cube sur stdout, /project) | 328 Mo | 512 Mo | 737 Mo |

Scratch (`outdir`) : quatre `.save` ≈ 101 Go ; si pw.x écrit aussi les 64 `.wfcN` distribués en fin de run (vu pour le nscf R5),
le pic est ≤ 200 Go. Quota utilisateur : 116 Go / 20 To, 5 902 / 1 M fichiers → OK. Purge Alliance à 60 jours → miroir
obligatoire (règle 3). Projet : répertoires de travail `super_cell/{15x15,18x18}` (Vks + pp.out + scf.out) ≈ 5,8 Go ; miroir
`qe_tmp_backup/defect_NxN_{d,p}/` : XML + charge-density + C.upf ≈ 1,2 Go, avec les wfc ≈ 101 Go ; quota du groupe
rrg-cotemich-ac 39 / 180 To, 4,67 M / 15 M fichiers → OK (nearline rrg : 311 Go / 100 To, possible pour les wfc). Copie versionnée :
inputs (22–31 Ko), submit, `scf.out` (≈ 50–100 Ko), générateur, pilote, rapport, tables, figures, extrait npz (< 1 Mo).

### 0.5 Manifeste proposé (à décider au GO)

| élément | emplacement | taille (4 cellules) | proposition |
|---|---|---|---|
| `data-file-schema.xml`, `charge-density.hdf5`, `C.upf` | `qe_tmp_backup/defect_NxN_{d,p}/defect_NxN_{d,p}.save/` (même chemin que 5…12) + md5 | ≈ 1,2 Go | miroir oui (rsync -a --no-o --no-g --open-noatime, `md5sum -c`) |
| `Vks_NxN_{d,p}` (pp.x) | répertoires de travail (/project), comme 5…12 | 3,3 Go | déjà hors scratch ; pas de miroir |
| `pp.out` (cube) | répertoires de travail (/project) | 2,5 Go | gardé comme 5…12, non versionné (> 5 Mo) ; listé |
| extrait npz de la fenêtre (valeurs propres, alignement, parité, w₂, w₁ ; bandes de [−3, +1] eV, p et d) | `R7_tailles_3m/d1/window_NxN_{p,d}.npz` + `article/R7_tailles_3m/` | < 1 Mo | oui (versionné) |
| `wfc1.hdf5` complets | `qe_tmp_backup/…` ou nearline | 101 Go | à décider : sans eux, parité / w₂ / w₁ ne sont plus recalculables après purge ; l'extrait couvre les grandeurs rapportées |
| `.wfcN` distribués (s'ils existent) | scratch | ≤ 100 Go | listés, non supprimés (aucune suppression dans R7) |

### 0.6 D1 (protocole R4/R5) et plan d'exécution

D1 = `cmd_c` de R5 restreint à la route « quadruplet » (pas de `.save` de maille pour N = 15, 18 ; inutile pour 3m) : E_D = moyenne
des quatre états de la parfaite les plus proches de E_F (étalement 0 pour 6/9/12) ; alignement Lu = `alignment.far_atom_alignment`
(moyennes sphériques 0,5 et 1,0 Å autour de l'atome le plus loin de la lacune ; décalage à 1,0 Å appliqué au spectre de la lacune) ;
fenêtre [−3, +1] eV ; `qe_gamma_io.read_wfc_gamma` + `mirror_parity_z` + `density_2d` → parité, w₂ (disque 2 Å), w₁ (1 Å) ; seuil
3 ⟨w₂⟩ de la parfaite de même N ; π = état impair localisé de plus grand w₂, doublet σ = deux états pairs de plus grand w₂ ;
comptages (pairs/impairs dans la fenêtre, localisés) ; `size_fits` 1/N et 1/N² sur cinq points (6, 9, 12, 15, 18) pour π et σ
(moyenne du doublet) ; figure ε − E_D contre 1/N (`memoire.mplstyle`, texte français, marine σ / orange π). Porte de régression :
6/9/12 doivent redonner R5 C à 1e-6 (π −1.065 / −0.737 / −0.551 ; σ +0.156 / +0.101 / +0.114 ; Lu −98.62 / −24.65 / −29.24 meV ;
seuils 0.1831 / 0.0812 / 0.0457). Référence trois points (R5) : π 1/N ε_∞ −0.046 (rms 0.007), 1/N² −0.409 (rms 0.024) ; σ 1/N
+0.055 (0.012), 1/N² +0.087 (0.011).

Code prévu (répertoire de travail, rien dans `src/`) : `r7_driver.py` (sous-commande `d1` ; réutilise `r5_driver.window_states_gamma`,
`load_pot_eV`, le schéma `sc_paths`, `r5_alignment_ext.{dirac_quadruplet,size_fits}` et les modules `qe_gamma_io`, `alignment`
commis dans 75c8656 ; routines de production intouchées) + `submit_r7.sh` (patron `submit_r5.sh`). Coût : R5 C sur huit tailles
2 min 38, MaxRSS 2,6 G (12×12 : 55/57 états, lecture 12 s) ; 18×18 ≈ 125 états × 1,54 M npw × 16 o = 3,1 G de coefficients +
FFT 540×540×192 → ≈ 8–10 G, quelques minutes par taille → job 16 cœurs, 64 G, 1 h.

Ordre après GO : (1) `sbatch submit.scf` dans `15x15/pristine`, `15x15/defective`, `18x18/pristine`, `18x18/defective` (JOBID
consigné ici) ; (2) surveillance sans relance (sacct, `scf.out` : « convergence has been achieved », « JOB DONE », « running with
the 2D cutoff ») ; (3) `sbatch submit.pp` ×4 après fin des scf (ou `--dependency=afterok` dès le GO si Greg le souhaite) ;
contrôles : E_F, étalement du quadruplet, taille des `Vks` ; (4) J1 `d1` ; (5) rapport (tables, figure, manifeste), copie
`article/R7_tailles_3m/` ; STOP.

**STOP — attente du GO** (quatre scf ; choix A/B des ressources ; pp.x chaînés ou non ; sort des wfc dans le manifeste).

## Phase GO (GO de Greg le 2026-09-25 : « GO ! oui fait le pp.x directement après »)

Décisions : ressources A (telles qu'écrites dans les `submit.scf`) ; pp.x chaînés en `--dependency=afterok` sur leur scf ; sort des
wfc dans le manifeste : à décider au rapport (non tranché par le GO).

### Jobs soumis le 2026-09-25 (≈ 15 h 30 ; `JOBID`)

| cellule | scf (pw.x, 64 rangs, 8 par nœud) | pp.x (afterok, 1 tâche, 32 G) |
|---|---|---|
| 15×15 parfaite | 21818658 (8G/cœur, 3 h) | 21818684 |
| 15×15 lacune | 21818685 (8G/cœur, 3 h) | 21818686 |
| 18×18 parfaite | 21818687 (16G/cœur, 6 h) | 21818688 |
| 18×18 lacune | 21818689 (16G/cœur, 6 h) | 21818690 |

D1 : pilote `r7_driver.py` (sous-commandes `d1`, `tables`) + `submit_r7.sh` (16 cœurs, 64 G, 1 h) écrits après le GO ;
régression 6/9/12 contre R5 C (tolérance 1e-6, `d1_reg/`) : job 21818859 ; D1 complet (6, 9, 12, 15, 18, `d1/`) à soumettre en
afterok des quatre pp.x une fois la régression PASS. Surveillance des jobs par sacct sans relance.


## Phase 0b — Extension 21×21, 24×24, 27×27 (demande de Greg le 2026-09-25 ; inputs préparés, **rien soumis, GO attendu**)

Inputs écrits par le même générateur (`make_inputs_r7.py --sizes 21,24,27`, diffs dans `inputs_diff_21_24_27.txt`) dans
`super_cell/{21x21,24x24,27x27}/{defective,pristine}/` ; régression `--check` toujours PASS (9×9 et 12×12 byte à byte, et les huit
fichiers déjà écrits des 15×15/18×18 reproduits à l'identique). Lacune A la plus proche du centre : 21×21 (10, 10) atome 441 ;
24×24 (12, 12) atome 601 ; 27×27 (13, 13) atome 729 ; vérification indépendante des fichiers : nat 882/881, 1152/1151, 1458/1457,
trois voisins à 1,4237 Å, a = 4,6597846766 bohr. Deux règles changent au-delà de 18×18, consignées ici :

**nbnd.** La règle +30 couvre 1,96 eV au-dessus de E_D à la 12×12 et, mesuré sur la 15×15 terminée, 1,77 eV (parfaite : 902 états
≤ E_D + 1 eV sur 930, marge 28 bandes ; lacune 26). La couverture décroît avec N (le nombre d'états par eV croît en N²) : à 21–27
elle tomberait vers 1,0–1,4 eV, à la limite de la fenêtre D1 [−3, +1]. Règle adoptée : **nbnd = N_occ + ⌈30 (N/12)²⌉** = +92 (21),
+120 (24), +152 (27), soit la couverture de la 12×12 (≈ 1,96 eV) ; surcoût 3,5 / 3,9 / 4,1 % sur nbnd (temps, mémoire, wfc).
nbnd : 21×21 1856/1854, 24×24 2424/2422, 27×27 3068/3066.

**Ressources.** Nœuds entiers (`--mem=0`) et 192 rangs pour les trois tailles : 192 = nr3, le nombre de plans FFT selon z
(c = 30 bohr), borne de la parallélisation R&G sans groupes de tâches. Extrapolation depuis la 15×15 mesurée (100 s par
itération à 64 rangs, moyenne des deux cellules ; MaxRSS rang 0 13,0 G = moyenne 2,65 G + 0,65 × wfc ; wfc 15,96 Go),
t ∝ N⁵ × nbnd, mémoire ∝ N⁴ × nbnd, efficacité 0,85 supposée de 64 à 192 rangs (non mesurée) :

| N | `submit.scf` | s/it à 64 rangs | s/it à 192 | mur 27 it (h) | limite | RAM totale | par nœud (rang 0 inclus) | wfc/cellule | 2 cellules | Vks | D1 (RAM) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 21 | 2 nœuds × 96, `--mem=0`, 6 h | 556 | 218 | 1,6 | 6 h | 0,66 To | 378 G | 63 Go | 127 Go | 1,3 Go | ≈ 57 G |
| 24 | 3 nœuds × 64, `--mem=0`, 8 h | 1089 | 427 | 3,2 | 8 h | 1,13 To | 455 G | 109 Go | 217 Go | 1,7 Go | ≈ 98 G |
| 27 | 4 nœuds × 48, `--mem=0`, 12 h | 1968 | 772 | 5,8 | 12 h | 1,81 To | 577 G | 174 Go | 349 Go | 2,2 Go | ≈ 158 G |

Les nœuds font 768 G ; à efficacité 0,6 le 27×27 prendrait 8,2 h (limite 12 h). Facturation : 2 / 3 / 4 nœuds pleins pendant le mur.
Un seul nœud ne suffit pas dès 21×21 (0,66 To + rang 0). `submit.pp` : 32 G conservé pour 21 et 24 (pp.x ≈ 16 / 21 G attendus,
12×12 : 5,3 G), 48 G pour 27 (≈ 27 G). D1 : `submit_r7.sh` avec `--mem` 128 / 192 / 256 G selon la taille (1 nœud), 2 h.

**Disque.** Scratch : 693 Go de `.save` pour les trois tailles (≈ 800 Go avec 15 et 18), quota 20 To ; projet : Vks + pp.out ≈ 10 Go.
Miroir : XML + charge-density ≈ 3 Go pour les trois ; wfc complets 693 Go supplémentaires (question du manifeste, non tranchée).

Ordre proposé après GO : scf 21 (p, d) → pp.x chaînés → puis 24 → puis 27 (ou les trois en parallèle : 9 nœuds), D1 à huit points
(6, 9, 12, 15, 18, 21, 24, 27) en afterok. Rien n'est soumis.

## Résultats (2026-09-25, 15 h 30 → 17 h 51 ; tous les jobs COMPLETED, aucune relance)

### Runs QE (64 rangs, 8 par nœud ; `scf.out`, sacct)

| cellule | job scf | it | mur | E_F (eV) | E_tot (Ry) | MaxRSS rang 0 / moyenne | `wfc1.hdf5` | job pp.x | mur pp.x / MaxRSS | `Vks` / `pp.out` |
|---|---|---|---|---|---|---|---|---|---|---|
| 15×15 parfaite | 21818658 | 19 | 36 min 20 s | −4.2198 | −5421.51446966 | 13.06 / 2.72 G | 15.96 Go | 21818684 | 2 min 26 / 8.2 G | 669 / 512 Mo |
| 15×15 lacune | 21818685 | 29 | 42 min 46 s | −4.2174 | −5408.86552118 | 13.03 / 2.65 G | 15.93 Go | 21818686 | 2 min 04 / 8.1 G | 669 / 512 Mo |
| 18×18 parfaite | 21818687 | 18 | 77 min 00 s | −4.2199 | −7807.01354247 | 30.93 / 5.47 G | 32.76 Go | 21818688 | 4 min 15 / 11.7 G | 963 / 737 Mo |
| 18×18 lacune | 21818689 | 30 | 105 min 06 s | −4.2188 | −7794.36212492 | 30.94 / 5.45 G | 32.71 Go | 21818690 | 3 min 13 / 11.5 G | 963 / 737 Mo |

Tous : PWSCF 7.5, « running with the 2D cutoff », « convergence has been achieved », « JOB DONE », grilles 450×450×192 et
540×540×192, npw (demi-sphère Γ) 1 071 681 et 1 542 x ; `charge-density.hdf5` 240 / 346 Mo ; aucun `.wfcN` distribué dans les
`outdir` ; scratch occupé par les quatre `.save` : 93 Go. Prévisions de phase 0 : murs 44–51 / 110–134 min (lacune) et 31–35 /
71–85 min (parfaite) → mesurés 43 / 105 et 36 / 77 min ; MaxRSS rang 0 ≈ 18 / 37 G → mesurés 13 / 31 G ; wfc 16,2 / 33,6 Go →
15,96 / 32,7 Go. Le rang 0 porte la moyenne + ≈ 0,65 × wfc au moment de l'écriture. Dernière bande − E_F : +1.750 / +1.735 (15×15),
+1.695 / +1.587 (18×18) : la règle +30 a couvert la fenêtre [−3, +1] avec ≥ 0,6 eV de marge.

### D1 (job 21819066, 4 min 26, MaxRSS 26,1 G ; lecture des états de la fenêtre 18×18 en 136 s ; `d1/`)

Porte de régression 6/9/12 contre R5 C : **PASS**, écart 0 sur les 42 grandeurs (déjà PASS en `d1_reg/`, job 21818859).
Tables complètes dans `d1/D1_tables.md` ; extraits npz `d1/window_NxN.npz` (valeurs propres complètes p et d, états de la fenêtre :
bandes, ε, ε alignée, ε − E_D, parité, w₂, w₁ ; E_D, E_F, décalages Lu, seuil, s_vac, grille, cellule).

| N | E_D quadruplet (eV) | E_F(p) − E_D | Lu 1,0 / 0,5 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ (meV) | seuil w₂ | fenêtre P / D (σ, π) | π quasi-lié : ε − E_D ; w₂ ; w₁ ; bande | doublet σ : ε − E_D ; w₂ ; w₁ | localisés σ / π | 2ᵉ w₂ impair |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | −4.23745 | +19.0 meV | −98.62 / −125.83 | +20.83 | 0.1831 | 19 (0, 19) / 21 (2, 19) | −1.0654 ; 0.284 ; 0.053 ; 140 | +0.1560 ; 0.711 ; 0.280 | 2 / 1 | — |
| 9 | −4.23847 | +19.0 | −24.65 / −21.59 | +9.09 | 0.0812 | 28 (0, 28) / 30 (2, 28) | −0.7373 ; 0.246 ; 0.047 ; 320 | +0.1014 ; 0.713 ; 0.281 | 2 / 4 | 0.114 |
| 12 | −4.23869 | +19.0 | −29.24 / −36.07 | +5.06 | 0.0457 | 55 (0, 55) / 57 (2, 55) | −0.5508 ; 0.223 ; 0.044 ; 572 | +0.1140 ; 0.713 ; 0.281 | 2 / 9 | 0.079 |
| 15 | −4.23878 | +19.0 | −15.64 / −14.54 | +1.42 | 0.0292 | 76 (0, 76) / 78 (2, 76) | −0.4593 ; 0.207 ; 0.041 ; 896 | +0.1033 ; 0.713 ; 0.281 | 2 / 10 | 0.099 |
| 18 | −4.23884 | +19.0 | −18.18 / −20.34 | −8.48 | 0.0203 | 127 (0, 127) / 129 (2, 127) | −0.3881 ; 0.194 ; 0.039 ; 1292 | +0.1070 ; 0.712 ; 0.281 | 2 / 15 | 0.053 |

Faits bruts :
- Quadruplet dégénéré à ≤ 3,2e-7 eV pour les cinq tailles ; E_D varie de −4,23745 (6) à −4,23884 eV (18), 0,06 meV entre 15 et 18 ;
  E_F(p) − E_D = +19,0 meV à toutes les tailles (lissage mv 0,01 Ry sur le quadruplet demi-rempli).
- Les états de la fenêtre de la parfaite sont tous impairs (π) ; la lacune ajoute exactement deux états pairs, le doublet σ,
  dégénéré à < 0,1 meV à toutes les tailles, avec w₂ = 0,712–0,713 et w₁ = 0,280–0,281 indépendants de N. Position :
  +0,156 (6), puis +0,101 / +0,114 / +0,103 / +0,107 pour 9 → 18 (dispersion 13 meV pour N ≥ 9).
- État π quasi-lié (impair de plus grand w₂) : −1,065 / −0,737 / −0,551 / −0,459 / −0,388 eV ; w₂ décroît 0,284 → 0,194 et w₁
  0,053 → 0,039 ; il reste l'impair de plus grand w₂ avec un facteur ≥ 2 sur le suivant (0,194 contre 0,053 à 18×18). Le nombre
  d'impairs au-dessus du seuil croît (1, 4, 9, 10, 15) avec le seuil 3 ⟨w₂⟩_P ∝ 1/N² ; les autres impairs localisés ont w₂ ≤ 0,114.
- États de la parfaite dans ]E_D, E_D + 1 eV] : 2, 2, 3, 2, 15 ; la première couronne de la grille repliée autour de K (|b|/N,
  0,163 Å⁻¹ à 18×18, ≈ 1 eV) entre dans la fenêtre à 18×18 seulement.
- Décalage Lu (1,0 Å) : −98,6 / −24,7 / −29,2 / −15,6 / −18,2 meV ; ⟨V_d⟩ − ⟨V_p⟩ 3D : +20,8 / +9,1 / +5,1 / +1,4 / −8,5 meV.

Ajustements ε(N) = ε_∞ + a/N^p sur la famille 3m à cinq points (eV ; `fig/size_3m`) :

| état | loi | ε_∞ | a | rms | max résidu | référence R5 à trois points (6, 9, 12) |
|---|---|---|---|---|---|---|
| π | 1/N | −0.0493 | −6.113 | 0.0055 | 0.0088 | −0.0458 (rms 0.0069) |
| π | 1/N² | −0.3490 | −26.665 | 0.0355 | 0.0591 | −0.4089 (rms 0.0241) |
| π | 1/N + 1/N² (supplément, hors prompt) | −0.0325 | −6.461 ; b = +1.55 | 0.0052 | — | — |
| σ (moyenne du doublet) | 1/N | +0.0750 | +0.427 | 0.0111 | 0.0211 | +0.0549 (rms 0.0123) |
| σ | 1/N² | +0.0944 | +2.008 | 0.0092 | 0.0178 | +0.0869 (rms 0.0107) |
| σ | 1/N + 1/N² (supplément) | +0.1570 | −1.278 ; b = +7.59 | 0.0060 | — | — |

Pour le π, le rms de la loi 1/N (5,5 meV, résidu max 8,8 meV) est 6,5 fois plus petit que celui de la loi 1/N² (35,5 meV), et
l'ajout du terme 1/N² ne réduit le rms que de 0,3 meV (b = +1,55 eV contre a = −6,46 eV) ; les deux points nouveaux (15, 18)
tombent sur la droite 1/N ajustée à trois points à 10 meV près. Pour le σ, les rms (6–11 meV) sont du même ordre que la
dispersion des quatre tailles N ≥ 9 (13 meV) et les trois lois donnent des ε_∞ de +0,075 à +0,157 eV : les cinq points ne
départagent pas les lois ; le point 6×6 (+0,156) est le seul hors de la bande +0,10–0,11. Figures : `fig/size_3m` (a) π,
(b) σ, points annotés par N, droites 1/N (trait plein) et 1/N² (tirets) ; `fig/localized_3m` : tous les états localisés contre 1/N.

### Fichiers de la campagne et manifeste

- Répertoire de travail `graphene/qe/defects/R7_tailles_3m/` : `README.md`, ce rapport, `make_inputs_r7.py`, `phase0_inventory.py`,
  `phase0_analysis.txt`, `inputs_diff_vs_12x12.txt`, `inputs_diff_21_24_27.txt`, `r7_driver.py`, `submit_r7.sh`, `r7_log.txt`,
  `JOBID`, `slurm-r7-*`, `d1_reg/` (régression), `d1/` (json, tables, 5 npz), `fig/`.
- Répertoires QE `super_cell/{15x15,18x18}/{defective,pristine}/` : inputs, `scf.out` (56–58 Ko), `pp.out` (512 / 737 Mo, non
  versionné), `Vks_NxN_{d,p}` (669 / 963 Mo), `pp.x_*.{out,err}` ; `super_cell/{21x21,24x24,27x27}/` : inputs seulement (§0b).
- Scratch `qe_tmp/defect_NxN_{d,p}/` : `.save` complets (93 Go) + `prefix.xml` ; rien d'autre.
- Miroir `qe_tmp_backup/defect_NxN_{d,p}/` (2026-09-25, 17 h 55) : `data-file-schema.xml`, `charge-density.hdf5`, `C.upf`,
  `prefix.xml` (rsync `-a -r --no-o --no-g --open-noatime --chmod=Dg+s`, `chgrp -R rrg-cotemich-ac`, 16 fichiers, 1,2 Go,
  `MD5SUMS_R7_2026-09-25.txt` ; `md5sum -c` 16/16 OK, md5 source = miroir vérifié fichier par fichier). **Les `wfc1.hdf5`
  (97,4 Go) ne sont pas miroités** : décision de Greg (miroir `qe_tmp_backup` ou nearline, ou extrait npz seulement).
- Copie versionnée `graphene-raman/article/R7_tailles_3m/` : rapport, README, scripts, tables, npz, figures, inputs et `scf.out`
  des cinq répertoires QE (jamais `.save`, `Vks`, `pp.out`, `slurm-*`, `JOBID`).
- Aucune suppression. Rien de commis (Greg).

**STOP.** Décisions restantes : (1) miroir des wfc 15/18 (97 Go) ; (2) GO ou non pour 21/24/27 (§0b : cascade ou parallèle) ;
(3) commit de `article/R7_tailles_3m/`.

## Phase GO 2 (GO de Greg le 2026-09-25, ≈ 18 h : « Miroir vers qe_tmp_backup et GO pour les plus grosses supercellules. Ne commit rien pour l'instant. »)

Amendement au §0b avant soumission : l'excès du rang 0 mesuré à la 18×18 est 0,78 × wfc (25,4 G pour 32,7 Go) et non 0,65 ;
avec 4 nœuds le 27×27 aurait ≈ 600–630 G par nœud sur 768 (marge 18 %). **27×27 passé à 6 nœuds × 32 rangs** (toujours 192 rangs,
≈ 450 G par nœud) ; `make_inputs_r7.py` mis à jour, `submit.scf` 27×27 régénérés, `--check` PASS. 21×21 (2 × 96, ≈ 390 G/nœud)
et 24×24 (3 × 64, ≈ 480 G/nœud) inchangés. Les six scf sont soumis en parallèle (11 nœuds au total), pp.x chaînés en afterok.

| cellule | scf | nœuds × rangs, limite | pp.x (afterok) |
|---|---|---|---|
| 21×21 parfaite | 21833299 | 2 × 96, 6 h | 21833300 |
| 21×21 lacune | 21833301 | 2 × 96, 6 h | 21833302 |
| 24×24 parfaite | 21833303 | 3 × 64, 8 h | 21833304 |
| 24×24 lacune | 21833305 | 3 × 64, 8 h | 21833306 |
| 27×27 parfaite | 21833307 | 6 × 32, 12 h | 21833308 |
| 27×27 lacune | 21833309 | 6 × 32, 12 h | 21833310 |

D1 à huit points (6, 9, 12, 15, 18, 21, 24, 27) : job 21833311, afterok des six pp.x, 16 cœurs, 256 G, 3 h, sorties `d1_8pts/`
(le `d1/` à cinq points est conservé). Miroir des `wfc1.hdf5` 15×15/18×18 (97,4 Go) vers `qe_tmp_backup/` : job 21833312
(`submit_mirror_wfc.sh` : rsync, chgrp, md5 source et miroir comparés, ajout à `MD5SUMS_R7_2026-09-25.txt`). Rien de commis
(consigne). Surveillance sans relance.

Miroir wfc 15×15/18×18 fait (job 21833312, 22 h 39 → 23 h 02 : rsync 97,4 Go en 8 min ≈ 200 Mo/s, md5 source 8 min, md5 miroir 5 min) :
**md5 source = miroir 4/4 OK**, `qe_tmp_backup/defect_{15x15,18x18}_{d,p}/` = 93 Go, `MD5SUMS_R7_2026-09-25.txt` 20 lignes (16 + 4 wfc),
groupe rrg-cotemich-ac. Les `.save` 15/18 sont donc entièrement miroités (règle 3) ; rien supprimé sur le scratch.

Premiers temps à 192 rangs (23 h 05) : 21×21 ≈ 330 s par itération (3 it en 18 min 40), 24×24 ≈ 700 s (1 it en 16 min 30, mise en
place comprise), 27×27 aucune itération après 21 min : efficacité ≈ 0,55 par rapport à 64 rangs, contre 0,85 supposé au §0b.
Projection : 21×21 lacune ≈ 2,5 h (limite 6 h), 24×24 lacune ≈ 5 h (8 h), 27×27 lacune ≈ 10 h pour 30 itérations (12 h) :
à surveiller sur les premières itérations (décision : laisser courir ou annuler et resoumettre avec 24 h).

## Résultats GO 2 — 21×21, 24×24, 27×27 (nuit du 25 au 26 septembre 2026 ; tous COMPLETED, aucune relance)

### Runs QE (192 rangs, nœuds entiers ; `scf.out`, sacct)

| cellule | job scf | nœuds × rangs | it | mur | E_F (eV) | E_tot (Ry) | MaxRSS rang 0 / moyenne | `wfc1.hdf5` | pp.x : mur / MaxRSS |
|---|---|---|---|---|---|---|---|---|---|
| 21×21 parfaite | 21833299 | 2 × 96 | 18 | 76 min 29 s | −4.2199 | −10626.23505960 | 49.3 / 3.6 G | 62.40 Go | 6 min 15 / 15.5 G |
| 21×21 lacune | 21833301 | 2 × 96 | 30 | 105 min 44 s | −4.2205 | −10613.58213655 | 52.8 / 3.6 G | 62.34 Go | 5 min 09 / 15.7 G |
| 24×24 parfaite | 21833303 | 3 × 64 | 20 | 163 min 17 s | −4.2199 | −13879.18001899 | 29.4 / 5.7 G | 106.44 Go | 6 min 15 / 20.8 G |
| 24×24 lacune | 21833305 | 3 × 64 | 33 | 245 min 40 s | −4.2222 | −13866.52609300 | 29.7 / 5.7 G | 106.35 Go | 9 min 11 / 20.8 G |
| 27×27 parfaite | 21833307 | 6 × 32 | 19 | 281 min 40 s | −4.2199 | −17565.84899959 | 142.2 / 9.5 G | 170.49 Go | 14 min 30 / 25.8 G |
| 27×27 lacune | 21833309 | 6 × 32 | 30 | 388 min 33 s | −4.2237 | −17553.19434454 | 167.9 / 9.6 G | 170.38 Go | 9 min 08 / 26.2 G |

Tous : PWSCF 7.5, 2D cutoff, convergence atteinte (ΔE ≤ 8e-11 Ry), JOB DONE ; grilles 630², 720², 810² × 192 ; `charge-density`
471 / 615 / 778 Mo ; aucun `.wfcN` ; scratch total R7 : 877 Go. Murs contre les limites : 1,8 / 6 h, 4,1 / 8 h, 6,5 / 12 h.
Efficacité 64 → 192 rangs (mur par itération rapporté à l'extrapolation N⁵ de la 15×15) ≈ 0,55–0,6, comme relevé à 23 h 05.
Mémoire : rang 0 = moyenne + 0,8–0,95 × wfc (27×27 : 168 G, nœud à ≈ 475 G sur 768) ; le passage à 6 nœuds était utile mais 4
auraient tenu (≈ 630 G). Dernière bande − E_F : +1.92 / +1.93, +2.00 / +1.94, +2.04 / +1.92 eV : la règle nbnd = N_occ + ⌈30 (N/12)²⌉
a donné la couverture visée (≈ 2 eV). pp.x : `Vks` 1.29 / 1.71 / 2.17 Go, `pp.out` 0.99 / 1.31 / 1.66 Go par cellule.

### D1 à huit points (job 21833311, 24 min 13, MaxRSS 109 G ; `d1_8pts/`, porte 6/9/12 contre R5 : PASS, écart 0)

| N | E_D quadruplet | Lu 1,0 / 0,5 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ (meV) | seuil w₂ | fenêtre P / D | π : ε − E_D ; w₂ ; w₁ ; bande | doublet σ : ε − E_D ; w₂ | localisés σ / π |
|---|---|---|---|---|---|---|---|---|
| 21 | −4.23886 | −13.77 / −13.44 | +1.03 | 0.0149 | 154 (0, 154) / 168 (2, 166) | −0.3437 ; 0.183 ; 0.036 ; 1760 | +0.1030 ; 0.711 | 2 / 23 |
| 24 | −4.23888 | −13.12 / −12.94 | −62.74 | 0.0114 | 217 (0, 217) / 220 (2, 218) | −0.3074 ; 0.174 ; 0.035 ; 2300 | +0.1026 ; 0.712 | 2 / 29 |
| 27 | −4.23890 | −12.77 / −12.69 | −99.91 | 0.0090 | 262 (0, 262) / 269 (2, 267) | −0.2785 ; 0.166 ; 0.034 ; 2912 | +0.1024 ; 0.712 | 2 / 34 |

Ajustements sur les huit points (6 → 27 ; `fig/size_3m_d1_8pts`) :

| état | loi | ε_∞ (eV) | a | rms (eV) | max résidu | trois points R5 (6, 9, 12) | cinq points (6 → 18) |
|---|---|---|---|---|---|---|---|
| π | 1/N | −0.0523 | −6.089 | 0.0045 | 0.0089 | −0.0458 | −0.0493 |
| π | 1/N² | −0.2963 | −29.345 | 0.0465 | 0.0786 | −0.4089 | −0.3490 |
| π | 1/N + 1/N² (supplément) | −0.0519 | −6.100 ; b = +0.06 | 0.0045 | — | — | −0.0325 |
| σ | 1/N | +0.0843 | +0.352 | 0.0095 | 0.0221 | +0.0549 | +0.0750 |
| σ | 1/N² | +0.0972 | +1.864 | 0.0075 | 0.0189 | +0.0869 | +0.0944 |
| σ | 1/N + 1/N² (supplément) | +0.1221 | −0.621 ; b = +4.86 | 0.0059 | — | — | +0.1570 |

Faits bruts :
- π quasi-lié (ε − E_D) : −1.065, −0.737, −0.551, −0.459, −0.388, −0.344, −0.307, −0.279 eV pour N = 6 → 27. La loi 1/N tient sur
  huit points avec un rms de 4,5 meV (résidu max 8,9 meV) ; le coefficient a = −6,09 eV est stable depuis trois points (−6,14) ;
  le terme 1/N² ajusté en supplément est nul (b = +0,06 eV, rms inchangé). La loi 1/N² seule a un rms de 46 meV. ε_∞ = −0,052 eV.
- Doublet σ : +0.1030, +0.1026, +0.1024 eV pour 21, 24, 27 (plateau +0,102–0,107 depuis N = 15 ; 6×6 +0,156 et 9×9 +0,101 hors
  plateau) ; dégénéré à < 0,1 meV ; w₂ 0,711–0,712, w₁ 0,281 à toutes les tailles. Les ajustements restent pilotés par 6 et 9
  (ε_∞ de +0,084 à +0,122 selon la loi, rms 6–10 meV).
- w₂ du π : 0,284 → 0,166 (6 → 27), w₁ 0,053 → 0,034 ; le π reste l'impair de plus grand w₂ (suivant ≤ 0,064 à 21×21).
- E_D : −4,23886 / −4,23888 / −4,23890 eV (0,04 meV entre 21 et 27) ; E_F(p) − E_D = +19,0 meV partout ; E_F de la lacune décroît
  avec N (−4,2205 / −4,2222 / −4,2237 contre −4,2199 pour la parfaite).
- Décalage Lu (1,0 Å) : −13,8 / −13,1 / −12,8 meV (0,5 Å à ≤ 0,4 meV près) ; la moyenne 3D ⟨V_d⟩ − ⟨V_p⟩ vaut +1,0 (21), −62,7 (24),
  −99,9 meV (27), sans lien avec le décalage local ; rapporté tel quel.
- Fenêtre [−3, +1] : 154 / 217 / 262 états (parfaite, tous impairs), 168 / 220 / 269 (lacune, dont exactement 2 pairs) ; impairs
  localisés au-dessus du seuil 3 ⟨w₂⟩_P : 23 / 29 / 34.

### Fichiers et manifeste (21/24/27)

- Répertoires QE `super_cell/{21x21,24x24,27x27}/{defective,pristine}/` : inputs, `scf.out`, `pp.out` (1–1,7 Go, non versionné),
  `Vks_NxN_{d,p}` (1,3–2,2 Go), `pp.x_*.{out,err}`. Scratch : six `.save` (678 Go de wfc + 3,7 Go), aucun `.wfcN`.
- `d1_8pts/` : `d1_results.json`, `D1_tables.md`, huit `window_NxN.npz` ; `fig/size_3m_d1_8pts`, `fig/localized_3m_d1_8pts`.
- Miroir : partie légère (XML, charge-density, C.upf, prefix.xml) et wfc (678 Go) vers `qe_tmp_backup/defect_NxN_{d,p}/`,
  selon la décision « miroir vers qe_tmp_backup » du GO 2 ; voir les lignes datées ci-dessous. Aucune suppression ; rien de commis.
- 2026-09-26, 09 h 05 : miroir léger 21/24/27 fait (XML, charge-density, C.upf, prefix.xml : 24 fichiers, `MD5SUMS_R7_2026-09-26.txt`,
  `md5sum -c` 24/24, source = miroir sur les 18 fichiers de données, groupe rrg-cotemich-ac). Miroir des six `wfc1.hdf5` (678 Go) :
  job 21850526 (`submit_mirror_wfc_21_24_27.sh`, 10 h ; rsync, chgrp, md5 source et miroir, ajout à `MD5SUMS_R7_2026-09-26.txt`).
