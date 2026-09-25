# EM1 — Couplage électron-photon : inventaire et éléments de position r(R)

Date : 2026-09-25. Série EM (couplage électron-photon, §2.5 du mémoire), distincte des P (mémoire) et des R (article).
Statut : **PRODUCTION**. Répertoire de travail : `graphene/qe/electron_photon/EM1_tb/` (hors dépôt).
Copie versionnée : `memoire/EM1_tb/` dans le dépôt (créée le 2026-09-25 sur GO de Greg ; règle 5 et tableau de CLAUDE.md complétés). Commit fait, non poussé.
Référence (LECTURE SEULE, rien écrit) : `graphene/qe/defects/unit_cell/27x27/`. Aucun push (Greg pousse lui-même).

Chemins abrégés : `$P` = `/home/gregb26/links/projects/rrg-cotemich-ac/gregb26`, `$S` = `/home/gregb26/links/scratch`.

## Résumé

- La wannierisation de référence du ch. 4 (maille unitaire, 27×27×1, nbnd 20, chaîne dense du 2026-09-04) contenait déjà
  `write_tb = true` : le `wannier_tb.dat` de référence existe depuis le 2026-09-04 (md5 `4a3f65cb…`, identique à la copie
  suivie par git `wannier/27x27/wannier_tb.dat`, sha256 = manifeste).
- `restart = plot` sur une copie du `.chk` (binaire Wannier90 3.1.0 du module `quantumespresso/7.5`, celui de la référence)
  reproduit `_tb.dat`, `_hr.dat`, `_wsvec.dat`, `_u.mat`, `_u_dis.mat` **à l'octet près hors ligne de date**. Aucun STOP.
- H(R) : même liste de R (741), mêmes ndegen, écart 0. Centres ⟨0i|r|0i⟩ = centres du .wout à 4×10⁻⁷ Å.
- r(R) **n'est pas hermitien** : max |r_ij(R) − r_ji(−R)*| = 2,5×10⁻³ Å (x), 1,7×10⁻³ Å (y), 10⁻¹⁴ Å (z), alors que H l'est à
  10⁻¹⁵ eV. C'est une propriété de la formule Eq. (44) de Wang-Yates-Souza-Vanderbilt 2006 employée par Wannier90 ; postw90
  prend explicitement la partie hermitienne de r(k) (commentaire dans `get_oper.F90`). Fait rapporté, décision à Greg.
- Binaire `~/bin/wannier90.x` (ifort) : SIGSEGV en `restart = plot`. Ne pas l'utiliser pour cela.

## Étape 0 — Inventaire (lecture seule)

### 0.1 Chemins

Chaîne d'origine : `$P/graphene/qe/defects/unit_cell/27x27/`, job SLURM 20212561 (`submit_dense_nb20.sh`, 2026-09-04 :
`pw.x < nscf.in` → `wannier90.x -pp wannier` → `pw2wannier90.x < pw2wann.in` → `wannier90.x wannier` → `finalize_wannier.py 27x27`).
Tous les fichiers ci-dessous sont présents (rien n'a disparu au ménage P13).

| Fichier (dans `unit_cell/27x27/`) | Taille (o) | mtime 2026-09-04 | md5 |
|---|---|---|---|
| `wannier.win` | 24 791 | 12:05:33 | `86d737657df7012d5930697c4efdc57a` |
| `wannier.chk` | 3 870 101 | 12:19:35 | `42311f8201958cb418f465e2b879c705` |
| `wannier.mmn` | 86 611 131 | 12:19:17 | `0326c0a0186d5cc8a09b75c134d27daf` |
| `wannier.amn` | 4 884 399 | 12:18:20 | `fe3484f0e18112780721bce88038378c` |
| `wannier.eig` | 568 620 | 12:19:17 | `9cb7e96a07f4d275121c09fccbb9bcb9` |
| `wannier_hr.dat` | 930 066 | 12:19:36 | `b47bd2ff89f5c806e70ed1cace7684be` |
| `wannier_tb.dat` (déjà présent) | 2 882 097 | 12:19:46 | `4a3f65cb962acbcbac90892867545108` (sha256 `baa17b88…c98be`) |
| `wannier_wsvec.dat` | 896 817 | 12:19:46 | `16ea5c13aa1fe930140df6e84880657c` |
| `wannier_u.mat` / `wannier_u_dis.mat` | 599 310 / 2 294 235 | 12:19:46 | — |
| `wannier.nnkp` / `wannier.bvec` / `wannier.wout` | 195 764 / 285 828 / 198 224 | 12:09 / 12:19 / 12:19 | `03ca117b…` / `4c85f040…` / — |
| `pw2wann.in` | 216 | 06:30 | prefix `defect_uc_dense_27`, outdir `$S/qe_tmp/defect_uc_dense_27`, seedname `wannier`, `write_mmn`/`write_amn` = .true., `write_unk` = .false. |
| `nscf.in` / `nscf.out` | 32 843 / 8 707 | 12:05 / 12:09 | 729 k explicites (grille 27×27×1 complète, `nosym`, `noinv`), nbnd 20, `diago_full_acc`, conv_thr 1e-14, ecutwfc 100 Ry, MV 0.01 Ry, `assume_isolated='2D'` ; PWSCF 7.5, 207 s, E_F = −4,2199 eV |

`.save` du nscf :
- scratch : `$S/qe_tmp/defect_uc_dense_27/{defect_uc_dense_27.save, defect_uc_dense_27.xml}` — 733 fichiers, 2,1 Go, fichiers datés du 2026-09-04 ;
- miroir : `$P/graphene/qe/qe_tmp_backup/defect_uc_dense_27/defect_uc_dense_27.save` — 732 fichiers, 2,2 Go, md5 dans
  `qe_tmp_backup/MD5SUMS_2026-09-17.txt` (732 entrées) ; 5 fichiers revérifiés aujourd'hui (xml, charge-density, C.upf, wfc1, wfc496) : OK.

Copie suivie par git dans le dépôt : `wannier/27x27/{wannier_tb.dat, wannier.eig, wannier.wout, wannier_u.mat, wannier_u_dis.mat,
wannier_manifest.json}` (déposée par `finalize_wannier.py`) ; md5 de `wannier_tb.dat` et `wannier.eig` identiques à la référence.

### 0.2 Version de wannier90.x

- `.wout` de référence : « Release: 3.1.0, 5th March 2020 », « Using CODATA 2006 constant values », « Running in serial (with serial executable) ».
- Deux binaires 3.1.0 disponibles aujourd'hui, bannières identiques :
  1. module `quantumespresso/7.5` (StdEnv/2023, `module restore qe`) : `/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/wannier90.x`,
     md5 `9abd4c5c9aec9715634d5a89a6e29cec`, gfortran + flexiblas ; c'est ce que `which wannier90.x` renvoie sous `module restore qe` ;
  2. `~/bin/wannier90.x` = `$P/codes/wannier90/wannier90.x` (même md5 `e51c0c821e7a293b646ec4317a303a91`, ifort + openblas, 2025-11-24).
- Discrimination par l'expérience (scratchpad, copies de .win/.chk/.eig) : le binaire 1 reproduit en `restart = plot` les cinq fichiers
  de référence à l'octet près hors ligne de date (§ a) ; le binaire 2 **plante (SIGSEGV)** après avoir écrit un `_hr.dat` différent (mise en
  forme ifort) et un `_u_dis.mat` vide. L'en-tête « list-directed » du `_tb.dat` de référence (`2.1354903421575360`, `0.0000000000000000`)
  est d'ailleurs du style gfortran. **Origine = binaire 1 ; EM1 utilise le binaire 1.**
- Note : `$P/codes/wannier90/src` est un checkout *develop* (merge du 2025-04-30, CODATA 2022) étiqueté `w90_version = '3.1.0'` ; il ne
  correspond exactement à aucun des deux binaires. Il n'a servi qu'à lire les routines d'écriture (`hamiltonian_write_tb`, `ws_write_vec`),
  dont les formats sont ceux observés dans les fichiers de référence.

### 0.3 Paramètres du `.win`

| Paramètre | Valeur | Remarque |
|---|---|---|
| `num_wann` | 5 | |
| `num_bands` | 20 | |
| `projections` | `C1:sp2;pz` / `C2:pz` | WF 1–3 = sp² de C1, WF 4 = p_z C1, WF 5 = p_z C2 |
| `dis_win_min` / `dis_win_max` | −25,0 / 15,0 eV | |
| `dis_froz_min` | **absent** → = `dis_win_min` = −25,0 eV | doc : « If dis_froz_max is given, then the default for dis_froz_min is dis_win_min » ; .wout l. 1013 : « Inner: −25.00000 to −1.74000 (eV) » |
| `dis_froz_max` | −1,74 eV | |
| `mp_grid` | 27 27 1 | 729 k listés dans le .win, ordre = `nscf.in` (k₂ le plus rapide) |
| `use_ws_distance` | **absent** → défaut `.true.` (3.x) | confirmé par l'en-tête de `wannier_wsvec.dat` : « with use_ws_distance=.true. » |
| `write_*` présents | `write_bvec = .true.`, `write_tb = true`, `write_u_matrices = true`, `write_hr = .true.` | `write_tb` et `write_hr` étaient donc déjà actifs |
| autres | `iprint 2`, `num_iter 5000`, `conv_tol 1E-12`, `conv_window 4`, `dis_num_iter 5000`, `dis_conv_tol 1E-12`, `num_print_cycles 10`, `dis_mix_ratio 0.5`, `guiding_centres = true` | |
| absents | `bands_plot`, `wannier_plot`, `fermi_surface_plot`, `restart` | rien à désactiver |
| cellule | `unit_cell_cart` en bohr : (4,0354919061, −2,3298923383, 0), (4,0354919061, 2,3298923383, 0), (0, 0, 30) ; `atoms_frac` C1 (1/3, 1/3, 0), C2 (2/3, 2/3, 0) | |

### 0.4 Centres et étalements finaux (.wout, « Final State », cycle 150, « Wannierisation convergence criteria satisfied », Δ < 10⁻¹² sur 4 itérations)

| WF | centre (Å) | étalement (Å²) | en regard |
|---|---|---|---|
| 1 | (1,067745, 0,616462, −0,000000) | 0,61094820 | milieu de liaison C1–C2(−a₂) = (1,067745, 0,616463, 0) |
| 2 | (1,067745, −0,616462, −0,000000) | 0,61094820 | milieu de liaison C1–C2(−a₁) = (1,067745, −0,616463, 0) |
| 3 | (2,135489, 0,000000, −0,000000) | 0,61094848 | milieu de liaison C1–C2 = (2,135490, 0, 0) |
| 4 | (1,423660, 0,000000, 0,000000) | 0,96695368 | atome C1 = (1,423660, 0, 0) |
| 5 | (2,847320, −0,000000, 0,000000) | 0,96695465 | atome C2 = (2,847320, 0, 0) |

Somme des étalements 3,76675321 Å² ; Ω_I 3,024611, Ω_D 0,013641, Ω_OD 0,728499, Ω_tot 3,766750 Å². Réseau (Å) : a₁ = (2,135490, −1,232926, 0),
a₂ = (2,135490, 1,232926, 0), a₃ = (0, 0, 15,875316). Distance C–C 1,42366 Å.

### 0.5 E_D au point K

K = (2/3, 1/3, 0) = k n° 496 (base 1) de la liste du .win (écart 4,7×10⁻⁹). `wannier.eig` (eV) à K, bandes 1…20 :
−16,864233, −16,864233, −14,829462, **−4,238895, −4,238895**, 6,420368, 8,662661, 8,662661, 9,265238, 10,691041, 10,691042, 11,261801,
11,261801, 11,568570, 11,909682, 12,916933, 12,916933, 13,071169, 13,071170, 14,017232.

| Grandeur | Valeur |
|---|---|
| E_D = (E₄ + E₅)/2 à K | **−4,238895 eV** (E₄ − E₅ = −1,8×10⁻⁷ eV) |
| dis_froz_max − E_D | **+2,498895 eV** |
| dis_froz_min − E_D | **−20,761105 eV** (dis_froz_min effectif = dis_win_min) |
| dis_win_min − E_D / dis_win_max − E_D | −20,761105 / +19,238895 eV |
| E_F (nscf, MV 0,01 Ry) − E_D | −4,2199 − (−4,238895) = +19,0 meV |
| bandes ≤ dis_froz_max par k | 4 à 5 (5 près de K : π* sous −1,74 eV) |

Cohérent avec `efermi_read = E_D = −4,2389 eV` de NOTES_EPW et avec E_D(SC 9×9 parfaite) = −4,23847 eV (R4).

## Étape 1 — Génération de `wannier_tb.dat` par `restart = plot`

1. Répertoire `EM1_tb/` : copies `wannier.win`, `wannier.chk`, `wannier.eig` (`cp -p`) ; `ref/` = liens symboliques vers les sorties de
   référence (`_hr.dat`, `_tb.dat`, `_wsvec.dat`, `_u.mat`, `_u_dis.mat`, `.wout`, `.bvec`) + copie `wannier.win.ref`. Rien écrit dans la
   référence (mtimes tous du 2026-09-04, vérifiés après le run).
2. Diff du `.win` (seule modification ; `write_tb`/`write_hr` étaient déjà présents et sont laissés tels quels ; aucun `*_plot` à désactiver) :

```diff
--- ref/wannier.win.ref
+++ wannier.win
@@ -22,0 +23 @@
+ restart = plot
```

3. Exécution sur le nœud de connexion (`module restore qe`, série, 1 cœur) :

```
binary: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/wannier90.x
Fri 25 Sep 2026 12:53:29 PM EDT   wannier90.x wannier   rc=0   real 0m10.313s
.wout : « Reading restart information from file wannier.chk » / « Restarting Wannier90 from plotting routines ... »
        « Time for plotting 9.967 (sec) » / « Total Execution Time 10.180 (sec) »
```

   Sorties : `wannier_tb.dat` (2 882 097 o), `wannier_hr.dat` (930 066), `wannier_wsvec.dat` (896 817), `wannier_u.mat` (599 310),
   `wannier_u_dis.mat` (2 294 235), `wannier.bvec` (285 828), `wannier.wout` (77 980), `run.log`. Le `.chk` n'est pas réécrit (octets identiques).
4. Pas d'échec, donc pas de repli. Pour mémoire : le repli aurait été une wannierisation complète depuis `.mmn`/`.amn` avec un `.win` identique
   (sans `restart`), non exécutée.

## Vérifications

Script : `em1_check.py` (numpy, `.venv` du dépôt) ; sorties `em1_check_2026-09-25.txt` et `em1_check_c_shells_2026-09-25.txt`.

### a. H(R) — jauge

| Comparaison | Résultat |
|---|---|
| liste de R (741, ordre inclus) : nouveau `_hr.dat` / bloc H du `_tb.dat` vs `_hr.dat` de référence | identique |
| ndegen (717 × 1, 24 × 2 ; Σ 1/ndegen = 729 = 27·27·1) | identique |
| max \|ΔH\| nouveau `_hr.dat` vs référence `_hr.dat` (F12.6 des deux côtés) | **0** |
| max \|ΔH\| bloc H du `_tb.dat` vs référence `_hr.dat` | 5,0×10⁻⁷ eV (= arrondi F12.6 du `_hr.dat` ; le `_tb.dat` est en E15.8) |
| max \|ΔH\| bloc H du `_tb.dat` vs référence `_tb.dat` | **0** |
| max \|Δr\| bloc r du `_tb.dat` vs référence `_tb.dat` | **0** |
| md5 hors ligne de date (nouveau = référence) | `_tb.dat` `0c9248e4…`, `_hr.dat` `bd6a870b…`, `_wsvec.dat` `6447a116…`, `_u.mat` `a3fc88c0…`, `_u_dis.mat` `6554cbec…`, `.bvec` `dc299427…` |

Même jauge que le ch. 4 (c'est le même `.chk`, relu). Énergies sur site ⟨0i|H|0i⟩ (eV) : −15,084664, −15,084664, −15,084667, −3,912484, −3,912484.

### b. Diagonale ⟨0i|r|0i⟩ (bloc R = 0 du `_tb.dat`) vs centres du .wout

| WF | ⟨0i|r|0i⟩ (Å) | .wout (Å) | écart (Å) |
|---|---|---|---|
| 1 | (1,0677454, 0,6164618, −9×10⁻¹²) | (1,067745, 0,616462, −0,000000) | (4,0×10⁻⁷, −1,9×10⁻⁷, −9×10⁻¹²) |
| 2 | (1,0677454, −0,6164618, −9×10⁻¹²) | (1,067745, −0,616462, −0,000000) | (4,0×10⁻⁷, 1,9×10⁻⁷, −9×10⁻¹²) |
| 3 | (2,1354891, 7×10⁻¹³, −1×10⁻¹¹) | (2,135489, 0, −0,000000) | (1,0×10⁻⁷, 7×10⁻¹³, −1×10⁻¹¹) |
| 4 | (1,4236600, 4×10⁻¹¹, 6×10⁻¹¹) | (1,423660, 0, 0) | (0, 4×10⁻¹¹, 6×10⁻¹¹) |
| 5 | (2,8473200, −3×10⁻¹¹, 6×10⁻¹¹) | (2,847320, −0,000000, 0) | (0, −3×10⁻¹¹, 6×10⁻¹¹) |

Écart max par composante : **(4,0×10⁻⁷, 1,9×10⁻⁷, 5,9×10⁻¹¹) Å** ; le .wout n'imprime que 6 décimales (±5×10⁻⁷). Partie imaginaire de la diagonale : 0.

### c. Hermiticité : max |r_ij(R) − r_ji(−R)*| sur tous (i, j, R), (i, R) ≠ (j, 0)

18 520 paires ; tout R a son −R dans la liste.

| Grandeur | x | y | z |
|---|---|---|---|
| max \|r_ij(R) − r_ji(−R)*\| (Å) | **2,542×10⁻³** | **1,652×10⁻³** | 1,0×10⁻¹⁴ |
| max \|H_ij(R) − H_ji(−R)*\| (eV) | 1,0×10⁻¹⁵ | | |

- Maximum atteint à R = (−1, 0, 0), i = 1, j = 3 : r_13(R) = (−0,046858, −0,086006, 0) vs r_31(−R)* = (−0,049400, −0,085563, 0).
- Bloc R = 0 seul (i ≠ j) : 3,4×10⁻¹⁵ Å (hermitien). Le défaut n'apparaît qu'entre cellules différentes.
- Norme de Frobenius sur tous les R, ‖partie anti-hermitienne‖/‖partie hermitienne‖ : **r 1,21×10⁻³**, H 7,3×10⁻¹⁷.
- Dans l'espace k, r(k) = Σ_R e^{ik·R} r(R)/ndegen : max |r(k) − r(k)†| = 3,3×10⁻³ (Γ), 4,9×10⁻³ (K), 4,3×10⁻³ Å (M) ; H(k) : ≤ 2×10⁻¹⁴ eV.

Par couche |R| (écart max, toute paire (i, j) ; dernière colonne = écart / max |r_ij| de la couche) :

| \|R\| (Å) | n_R | x | y | max \|r_ij(R)\| (Å) | ratio |
|---|---|---|---|---|---|
| 0 | 1 | 2,9×10⁻¹⁵ | 3,4×10⁻¹⁵ | 2,20×10⁻¹ | 0 |
| 2,4659 | 6 | 2,54×10⁻³ | 1,37×10⁻³ | 2,20×10⁻¹ | 0,012 |
| 4,2710 | 6 | 7,5×10⁻⁴ | 2,4×10⁻⁴ | 8,64×10⁻³ | 0,087 |
| 4,9317 | 6 | 1,68×10⁻³ | 1,65×10⁻³ | 7,59×10⁻³ | 0,22 |
| 6,5240 | 12 | 5,1×10⁻⁴ | 8,0×10⁻⁴ | 2,55×10⁻³ | 0,31 |
| 7,3976 | 6 | 1,02×10⁻³ | 1,11×10⁻³ | 1,45×10⁻³ | 0,77 |
| 8,5420 | 6 | 2,8×10⁻⁴ | 2,7×10⁻⁴ | 3,3×10⁻⁴ | 0,85 |
| 9,8634 | 6 | 5,2×10⁻⁴ | 5,8×10⁻⁴ | 4,2×10⁻⁴ | 1,41 |
| 36,4913 | 8 | 8,3×10⁻⁵ | 1,5×10⁻⁵ | 4,7×10⁻⁵ | 1,75 |
| 38,4388 | 2 | 7,9×10⁻⁵ | 5×10⁻¹¹ | 4,2×10⁻⁵ | 1,89 |

Contexte (source Wannier90, `postw90/get_oper.F90`, checkout develop, l. 715–725) : « Since Eq.(44) WYSV06 does not preserve the Hermiticity
of the Berry potential matrix, take Hermitean part » — postw90 remplace r(k) par ½[r(k) + r(k)†] avant sa transformée vers R. Le `_tb.dat`
contient la quantité brute d'Eq. (44) (formule aux différences finies i Σ_b w_b b ⟨u_jk|u_i,k+b⟩), sans symétrisation. Décision à Greg
(symétriser r(k) comme postw90, ou r(R) ↔ r(−R)†, ou garder brut) ; le rapport ne tranche pas.

### d. Décroissance : max_{i≠j} |r_ij(R)| et max |H_ij(R)| selon |R|

| \|R\| (Å) | n_R | max \|r_ij\|, i ≠ j (Å) | max \|r_ii\| (Å) | max \|H_ij\|, i ≠ j (eV) | max \|H_ij\|, tous (eV) |
|---|---|---|---|---|---|
| 0 | 1 | 2,204×10⁻¹ | 2,847 (centre) | 2,909 | 15,08 |
| 2,4659 | 6 | 2,204×10⁻¹ | 3,35×10⁻² | 2,909 | 2,909 |
| 4,2710 | 6 | 8,64×10⁻³ | 3,9×10⁻⁷ | 2,69×10⁻¹ | 2,69×10⁻¹ |
| 4,9317 | 6 | 7,59×10⁻³ | 2,43×10⁻³ | 1,47×10⁻¹ | 1,47×10⁻¹ |
| 6,5240 | 12 | 2,55×10⁻³ | 1,22×10⁻³ | 1,99×10⁻² | 1,99×10⁻² |
| 7,3976 | 6 | 1,45×10⁻³ | 2,5×10⁻⁴ | 1,98×10⁻² | 1,98×10⁻² |
| 8,5420 | 6 | 3,3×10⁻⁴ | 4×10⁻⁷ | 3,75×10⁻³ | 3,75×10⁻³ |
| 8,8908 | 12 | 3,8×10⁻⁴ | 2,1×10⁻⁵ | 1,00×10⁻³ | 1,00×10⁻³ |
| 9,8634 | 6 | 4,2×10⁻⁴ | 2,1×10⁻⁴ | 3,55×10⁻³ | 3,55×10⁻³ |
| 10,7484 | 12 | 1,9×10⁻⁴ | 6,3×10⁻⁵ | 6,8×10⁻⁴ | 6,8×10⁻⁴ |
| 12,3293 | 6 | 1,4×10⁻⁴ | 3,2×10⁻⁵ | 7,4×10⁻⁴ | 7,4×10⁻⁴ |
| 13,0481 | 12 | 1,1×10⁻⁴ | 3,8×10⁻⁵ | 4,1×10⁻⁴ | 4,1×10⁻⁴ |
| … (70 couches) | | | | | |
| 36,3243 | 12 | 4,3×10⁻⁵ | 6,1×10⁻⁶ | 1,4×10⁻⁵ | 1,4×10⁻⁵ |
| 38,4388 (bord WS) | 2 | 4,2×10⁻⁵ | 5×10⁻¹² | 5,8×10⁻⁶ | 9,0×10⁻⁶ |

Bloc R = 0, |r_ij| max sur les composantes (Å) : sp²–sp² 0,099 (entre WF 1 et 2) et 0,086 (1–3, 2–3) ; sp²–p_z(C1) 0,220 ; sp²–p_z(C2) 0,024 (WF 1, 2) et 0,220 (WF 3) ; p_z–p_z 0,0018.
Rapport bord/centre : r 4×10⁻⁵/0,22 ≈ 2×10⁻⁴ ; H 1,4×10⁻⁵/2,9 ≈ 5×10⁻⁶.

### e. Format exact de `wannier_tb.dat` (Wannier90 3.1.0, `hamiltonian_write_tb`)

Unités : Å (réseau, r) et eV (H). 40 070 lignes = 56 lignes d'en-tête + 741 × 27 (bloc H) + 741 × 27 (bloc r).

1. Ligne 1 : date (list-directed, espace initial) : ` written on 25Sep2026 at 12:53:29`.
2. Lignes 2–4 : a₁, a₂, a₃ en Å (list-directed) :
   ```
      2.1354903421575360       -1.2329259238968218        0.0000000000000000
      2.1354903421575360        1.2329259238968218        0.0000000000000000
      0.0000000000000000        0.0000000000000000        15.875316257699998
   ```
3. Ligne 5 : `num_wann` (`           5`) ; ligne 6 : `nrpts` (`         741`).
4. Lignes 7–56 : `ndegen(1..741)`, 15 entiers par ligne (format `15I5`), 49 lignes pleines + 1 de 6.
5. **Bloc H** (lignes 57–20063), pour chaque R dans l'ordre de la liste : une ligne vide, la ligne `R₁ R₂ R₃` (`3I5`), puis 25 lignes
   `j i Re Im` (format `2I5,3x,2(E15.8,1x)`), **j (1ʳᵉ colonne) varie le plus vite** ; valeur = ⟨0j|H|Ri⟩ = `ham_r(j,i,R)` (WF j dans la
   cellule 0, WF i dans la cellule R). Premières lignes :
   ```
   (ligne 57 vide)
     -17    8    0
       1    1    0.12754893E-08 -0.99226183E-14
       2    1    0.32645942E-09 -0.55719318E-14
       3    1   -0.10511129E-09 -0.19365586E-13
       4    1    0.58306842E-13 -0.10866184E-13
   ```
   La ligne « 0 0 0 » du bloc H (R = 0) est la ligne 10 048 (ligne vide en 10 047). Même contenu que `_hr.dat` (dont les lignes sont `R₁ R₂ R₃ j i Re Im`, `F12.6`, sans ligne vide).
6. **Bloc r** (lignes 20064–40070), même structure : ligne vide, `R₁ R₂ R₃`, puis 25 lignes
   `j i Re(x) Im(x) Re(y) Im(y) Re(z) Im(z)` (format `2I5,3x,6(E15.8,1x)`), valeur = ⟨0j|r|Ri⟩ en Å, même convention d'indices que H.
   Premières lignes :
   ```
   (ligne 20064 vide)
     -17    8    0
       1    1    0.80297721E-13 -0.19815367E-13 -0.13060570E-14 -0.71805246E-13  0.34461386E-12  0.10396751E-12
       2    1   -0.20627583E-10 -0.78507922E-13  0.11805305E-10 -0.10052321E-12  0.11516010E-12 -0.12737978E-12
       3    1    0.75379952E-11  0.15181538E-12  0.14914294E-10  0.38676512E-13 -0.14242180E-12 -0.12249651E-12
   ```
   La ligne « 0 0 0 » du bloc r (R = 0) est la ligne 30 055 ; sa diagonale = centres (§ b).
7. **ndegen** : ni H ni r ne sont divisés par ndegen dans le fichier. Source : `ham_r` est accumulé avec `fac = exp(−i k·R)/num_kpts`
   (`hamiltonian.F90` l. 412–413) et `pos_r` avec le même `fac` (l. 949–950) ; `_hr.dat` écrit le même `ham_r` (vérifié : bloc H du
   `_tb.dat` = `_hr.dat` à l'arrondi F12.6 près). L'interpolation doit donc appliquer **la même règle aux deux blocs** :
   X(k) = Σ_R e^{i k·R} X(R)/ndegen(R), comme `read_w90_HR` + `Hwr_to_Hwk` le font pour H dans le dépôt.
8. Formules (source, `hamiltonian_write_tb`) : i = j → ⟨0i|r|Ri⟩ = −Σ_k (e^{−ik·R}/N_k) Σ_b w_b **b** Im ln M_ii(k,b) (Eq. 32 de
   Marzari-Vanderbilt 1997 pour R = 0, Eq. 44 WYSV06 modifiée par Eqs. 27, 29 de MV97 sinon) ; i ≠ j → Eq. 44 WYSV06 :
   i Σ_k (e^{−ik·R}/N_k) Σ_b w_b **b** M_ji(k,b), avec M_ji(k,b) = ⟨u_jk|u_i,k+b⟩ (jauge Wannier). Non hermitien par construction (§ c).
9. **`wannier_wsvec.dat`** (18 525 entrées = 741 × 5 × 5) : ligne 1 `## written on … with use_ws_distance=.true.` ; puis, pour chaque R
   (ordre de la liste), i, j (i puis j en boucle interne) : `R₁ R₂ R₃ i j` (`5I5`), `ndeg` (`I5`), puis `ndeg` lignes `T₁ T₂ T₃` (vecteurs
   de super-réseau à AJOUTER à R). 462 entrées avec ndeg > 1, 806 entrées avec au moins un T ≠ 0. Sur ces 806 entrées : max |r_ij(R)| =
   6,8×10⁻⁵ Å, max |H_ij(R)| = 1,4×10⁻⁵ eV (bord de la cellule de Wigner-Seitz 27×27).
   **S'applique-t-il à r ?** Doc (`files.md`, section `seedname_wsvec.dat`) : fichier écrit « if write_hr = true or write_rmn = true or
   write_tb = true », les T « should be added to the R vector to obtain the correct centre of the Wannier function that underlies a given
   matrix element (e.g. the Hamiltonian matrix elements in seedname_hr.dat) » ; docstring de `ws_write_vec` : « to be added to R vector in
   seedname_hr.dat, seedname_rmn.dat, etc. ». Donc oui, par construction il est prévu pour tout ⟨0i|O|Rj⟩, r compris. Dans le source
   *develop* sur disque, postw90 applique les mêmes décalages à AA_R (`get_oper.F90` l. 746, `operator_wigner_setup`) ; le source exact de
   la 3.1.0 n'est pas sur disque, cette dernière affirmation n'est donc vérifiée que sur le checkout develop. Rappel : la chaîne du ch. 4
   n'applique **pas** `use_ws_distance` à H (CLAUDE.md) ; l'amplitude sur les entrées concernées est ≤ 6,8×10⁻⁵ Å.

## md5 (fichiers produits, `MD5SUMS_EM1_2026-09-25.txt`)

| Fichier | md5 complet | md5 hors ligne de date (= référence) |
|---|---|---|
| `wannier_tb.dat` | `b3f71b9b368c50deda00c87542cba138` (sha256 `4460e83d…0bd8`) | `0c9248e4f6b11e94c50597028c60780a` |
| `wannier_wsvec.dat` | `5cd3d8cf7f63096d4db6dfc52733d9c7` | `6447a116494f4ca7aae1bef62226b06c` |
| `wannier_hr.dat` | `a596193ec439a3913fdbc3a31b56feb5` | `bd6a870bc5a236e31d2085f18bdacb13` |
| `wannier_u.mat` / `wannier_u_dis.mat` | `79689ea6…` / `2b394569…` | `a3fc88c0…` / `6554cbec…` |
| `wannier.win` (copie + restart) | `7ff04c67d6d577096508ed521ea4e2c2` | (référence `86d73765…`) |
| `wannier.chk` / `wannier.eig` (copies) | `42311f82…` / `9cb7e96a…` | = référence |

## Emplacements, statut, miroir

- **Répertoire de travail** (règle 5) : `graphene/qe/electron_photon/EM1_tb/`, à côté de `graphene/qe/epw/` (chapitre 5) — la série EM est
  un nouveau volet du mémoire, comme EPW ; il contient tout (inputs, `.chk`/`.eig` copiés, sorties, script, rapport, md5).
- **Copie versionnée** : CLAUDE.md règle 5 prévoit `article/<campagne>/` (article) ou « le répertoire du chapitre (mémoire) », mais aucun
  répertoire de chapitre n'existe dans le dépôt (ch. 4 = `results/M`, `wannier/`, `NOTES_TGAMMA.md` ; ch. 5 = `results/epw`, `scripts/epw_*`,
  `NOTES_EPW*.md`). **Proposition** : `memoire/EM1_tb/` (symétrique de `article/`), avec `README.md`, `EM1_rapport.md`, `em1_check.py`,
  `em1_check_*.txt`, `wannier.win` + `wannier.win.diff`, `run.log`, `MD5SUMS_EM1_2026-09-25.txt`, `wannier_wsvec.dat` (0,9 Mo, pas encore
  versionné ailleurs) ; sans `wannier_tb.dat` (déjà suivi, identique, sous `wannier/27x27/`), sans `_hr.dat` (redondant), sans `.chk`, `.eig`
  (suivi), `.wout`, `.bvec`, `_u*.mat` (suivis). **Fait le 2026-09-25** (GO de Greg) : répertoire créé, règle 5 de CLAUDE.md complétée (« `memoire/<campagne>/` (mémoire) ») et ligne EM1 ajoutée au tableau des campagnes.
- **Commit** (fait le 2026-09-25, non poussé) :
  ```
  EM1 : éléments de position r(R) de la wannierisation de référence 27×27 (restart = plot) — inventaire, vérifications a–e, format du _tb.dat, rapport
  ```
  (`memoire/EM1_tb/` + CLAUDE.md ; copie faite par `cp -p` depuis le répertoire de travail, `wannier.win.diff` = `diff -u` référence → copie).
- **Statut PRODUCTION** : aucun `.save` produit, rien sur le scratch ; toutes les sorties (≤ 3,9 Mo) sont dans le répertoire de travail sur
  `/project`. Le `.save` du nscf d'origine est déjà miroité (`qe_tmp_backup/defect_uc_dense_27/`, md5 2026-09-17, 5 fichiers revérifiés).
- Dépôt `graphene-raman` : propre avant EM1 (HEAD `75c8656`, R4 terminé) ; seul ajout = `memoire/EM1_tb/` + CLAUDE.md ; rien du mémoire (ch. 4/5) touché.

## Fichiers du répertoire de travail

`README.md`, `EM1_rapport.md`, `wannier.win` (modifié), `wannier.chk`, `wannier.eig` (copies), `ref/` (liens + `wannier.win.ref`),
`wannier_tb.dat`, `wannier_hr.dat`, `wannier_wsvec.dat`, `wannier_u.mat`, `wannier_u_dis.mat`, `wannier.bvec`, `wannier.wout`,
`wannier90_restart_plot.log` (vide : W90 n'écrit que dans le .wout), `run.log`, `em1_check.py`, `em1_check_2026-09-25.txt`,
`em1_check_c_shells_2026-09-25.txt`, `MD5SUMS_EM1_2026-09-25.txt`.

STOP après ce rapport : pas de code d'interpolation, pas de nscf ni de bands.x (EM2), pas de figure, rien dans le mémoire.
