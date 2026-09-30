# Évaluation de t / Γ (matrice t locale de Wannier) — paramètres réels et contrôles implémentés

Transcription pour le mémoire (chapitre 4, §4.1.5 et section spectrale), faite le 2026-09-18 à partir des fichiers réels :
provenances `fichier:ligne` re-grep ce jour ; valeurs de production lues dans `results/M/*.npz`, `results/M/*.csv` et
`results/M/logs/*.out` (jamais recalculées ici, sauf le test `test_local_green_batch.py`, 10 s, et le comptage d'états
de la grille interne, indiqué comme tel). Chemins relatifs à la racine du dépôt (`graphene-raman`, ex `ab-initio-defects` ; ce fichier est
dans `notes/` depuis le 2026-09-30, il était à la racine). Même rôle que `NOTES_EPW.md` §1f pour
§5.4 : une ligne par paramètre, une ligne par contrôle, pour confronter la liste des dix contrôles du mémoire à ce que le
code fait réellement. **Les jugements sont laissés au lecteur ; ce document constate.**

**Mise à jour R6 (2026-09-27).** La campagne R6 a corrigé la normalisation de la partie locale de M : M^L de production était en
norme super-cellule, trop petit d'un facteur N_cells = N² (constat R5-A.2). Les M corrigés (v2 : M2 = N_cells·M^L + M^NL, sidecar
`M_normalization = v2`) et tous les produits sont dans `results/M2/` (`config/production.json` : `results_dir`) ; `results/M/` est gelé
(README). Après R6, les §1, §2, §3 et §6 donnaient les valeurs v2 ; les §4 et §5 sont l'état
du 2026-09-18 (v1), inchangés ; le §7 résume R6. Rapport : `graphene/qe/defects/R6_production_corrigee/R6_rapport.md` (copie
`campagnes/R/R6_production_corrigee/`) ; table complète ancien → nouveau : `campagnes/R/R6_production_corrigee/etape3/table_v1_v2.md`. Le dépôt
s'appelle `graphene-raman` depuis le 2026-09-24 (ex `ab-initio-defects`) ; les provenances `fichier:ligne` datent du 2026-09-18.

**Mise à jour R10 (2026-09-30) — état final.** La production du chapitre 4 est `results/M2_plateau/` (`config/production.json` :
`results_dir`) : matrices M2 inchangées (`results/M2/`, clé `matrices_dir` ; `results/M/` et `results/M2/` gelés, `results_dir_frozen`) et
**alignement de Kumagai–Oba** du potentiel de défaut (bloc `alignment`, symbole du mémoire ΔV_PA^(N) = C_N ; définition au §1, ligne
« Alignement du potentiel », et au §8). Nomenclature : « non aligné » (= « tel quel » des archives, `results/M2`), « alignement à site unique »
(= « Lu »), « alignement de Kumagai–Oba » (= « plateau (i) »), « final » = v2 + Kumagai–Oba. Dans ce document, **la valeur en tête est la
valeur finale** ; entre crochets, seulement « [v1 : …] » (ce que dit le mémoire actuel) ; la valeur non alignée n'apparaît qu'au §8
(alignement), pour Re Σ et pour la grandeur de Kaasbjerg (Fig. 3) près de K. Table complète v1 → final, chiffres nouveaux et nomenclature :
`campagnes/M/ch4/table_v1_final.md` (pilote `campagnes/M/ch4/ch4_chiffres.py`) ; régions d'alignement : `campagnes/M/ch4/alignement_regions.md`.
Rapport R10 : `graphene/qe/defects/R10_plateau/R10_rapport.md` (copie `campagnes/R/R10_plateau/`) ; table v2 → final :
`campagnes/R/R10_plateau/c/table_v2_plateau.md`.

## 0. Ce que le script fait, dans l'ordre (chaîne de production, taille de référence 9×9)

Scripts : `scripts/t/compute_spectral_wannier.py` (carte niveau 1 : R_cut × grille × η), `scripts/t/resonance_metrics.py`
(Γ(ε), Born, δρ, T̄ en K), `scripts/t/resonance_criteria.py` (det, valeurs propres, règle de somme, Γ à c = 0,1 %),
`scripts/t/rcut_resigma.py` (Re/Im Σ par R_cut). Tous chargent `config/production.json` (`config.py:10–19`, clés
obligatoires l. 13) et passent la porte de jauge avant toute lecture.

1. **Config figée** (`config/production.json`, 2026-09-05) → 2. **porte de jauge** sha256 (`wannier_provenance.py:42–68`) →
3. **chargement de M** avec sidecar `bloch_norm = unit_cell` et `units = hartree`, conversion Ha → eV **une seule fois**
(`matrix_io.py:50–78`, `HA2EV = 27.211386245988` l. 24) → 4. **rotation de jauge** M_W = V† M V, V = U_dis·U
(`wannier_interpolation.py:11–25`), puis **double TF** vers M_W(R, R′) sur la boîte MP-duale 27×27 (l. 56–70 ; pas de
`ndegen`, pas de `use_ws_distance`, cohérent avec H) → 5. **recentrage** du défaut à R₀ = 0 (`local_tmatrix.py:79–94`,
seul le label R bouge) → 6. **garde-fou a** : poids ‖M_W(R, 0)‖ maximal en R = 0, sinon `AssertionError` (l. 36–51) →
7. **V_loc = P M_W P** sur les mailles |R| ≤ R_cut (`compute_spectral_wannier.py:101`, `resonance_metrics.py:36`) et
**garde-fou b** : hermiticité de V_loc, résidu > 1e-10 ⇒ symétrisation et message (`local_tmatrix.py:54–76`) →
8. **E_D** = milieu du gap minimal π/π* sur une grille 90×90 (`resonance_metrics.py:47–48`, `compute_spectral_wannier.py:92–94`) →
9. **g₀(ε)** = bloc local de la fonction de Green du réseau, somme sur la grille interne N_k^int × N_k^int
(`local_tmatrix.py:122–156`, `local_green_batch`, restructuration exacte : différences R − R′ + base propre de H(k)),
sur une grille d'énergie de pas η/ne_per_eta (`resonance_metrics.py:55`) →
10. **t(ε) = V_loc [1 − g₀(ε) V_loc]⁻¹** (`local_tmatrix.py:159–162` ; `resonance_metrics.py:59`) →
11. **Γ_nk = −2 Im ⟨nk| t(ε_nk) |nk⟩** avec ⟨wR|nk⟩ = U_wn(k) e^{2πik·R} (`local_tmatrix.py:262–268`,
`resonance_metrics.py:66–70`), état sur couche évalué au point d'énergie le plus proche (`local_tmatrix.py:266` ;
`resonance_metrics.py:67`, `np.rint`), états hors fenêtre ±e_window → NaN (l. 250–255) → 12. **positivité** Γ ≥ −1e-8 sinon
`AssertionError` (`local_tmatrix.py:269–271` ; `rcut_resigma.py:55`) → 13. **agrégats** : médiane de |Γ|·N_cells sur les
états de la fenêtre (`compute_spectral_wannier.py:109`), position de résonance = argmax de Γ sur les états à |ε − E_D| ≤ 1,5 eV
(l. 112), courbes Γ(ε) = moyenne lorentzienne des états (`resonance_metrics.py:75–79`).

**Depuis R10 (2026-09-30)** : les étapes 4 à 6 sont faites par un seul appel, `local_tmatrix.defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N)` :
rotation de jauge, double TF, **soustraction de C_N** (alignement de Kumagai–Oba, `config.alignment_C`) sur la diagonale M_W(R, R) des N² mailles de
la boîte (étiquettes brutes), puis recentrage et garde-fou a. Les matrices sont lues dans `matrices_dir` (`results/M2`), les produits écrits dans
`results_dir` (`results/M2_plateau`). La DOS moyennée sur le désordre et A_k (Kaasbjerg, Fig. 13–14) ne passent pas par ces scripts : module
`defects/many_body/disorder_average.py`, pilote `campagnes/R/R8_kaasbjerg/r8_driver.py` (§8). Les `fichier:ligne` ci-dessus datent du 2026-09-18 et ont
dérivé de quelques lignes ; carte à jour de la chaîne, fonction par fonction : inventaire du 2026-09-30, §3.

## 1. Paramètres réels (taille de référence 9×9)

| paramètre | valeur | fichier:ligne / preuve |
|---|---|---|
| Super-cellule de référence, M | 9×9 (N_cells = 81), lacune sur le sous-réseau A, s_red = (13/27, 13/27) ; M dense zero-padded p = 3, D = 27, forme (20, 729, 20, 729), Hartree, norme `unit_cell` pour M^L et M^NL (`M_normalization = v2`, R6 : M2 = N_cells·M^L + M^NL) ; matrices brutes, l'alignement de Kumagai–Oba est appliqué en base de Wannier (ligne « Alignement du potentiel », §8) | `config/production.json` (`matrices_dir`, `alignment`) ; `results/M2/M_dense_9x9.json` (sidecar) [v1 : `results/M/M_dense_9x9.json`, M^L en norme super-cellule] |
| Cellule primitive dense (nscf) | grille 27×27×1 = 729 k, `nosym` + `noinv` (grille complète), nbnd = 20, `diago_full_acc`, ecutwfc 50 (unités Ha du XML, soit 100 Ry), ecutrho 200, smearing m-v degauss 5,0e-3 (Ha du XML), conv_thr 5e-15, `assume_isolated = 2D` | `scratch/qe_tmp/defect_uc_dense_27/defect_uc_dense_27.save/data-file-schema.xml` (`<nk>`, `<nbnd>`, `<nosym>`, `<noinv>`, `<ecutwfc>`, `<smearing>`, `<assume_isolated>`) ; `config/production.json:16` (nbnd_dense 20) |
| Wannierisation 27×27 | num_wann = 5 (C1 : sp², p_z ; C2 : p_z), 20 bandes, fenêtre externe [−25, 15] eV, fenêtre gelée [−25, −1,74] eV (= E_D + 2,5 eV), Ω_I = 3,0246, Ω_D = 0,0136, Ω_OD = 0,7285, Ω_tot = 3,7668 Å² ; étalements 0,611 (×3, sp²) et 0,967 Å² (×2, p_z) ; run_id `baa17b88b1b69e51` | `wannier/27x27/wannier.wout` l. 125 (grille), 863 (num_wann), 1012–1013 (fenêtres), 2553–2555 (Ω), « Final State » ; `config/production.json` clé `wannier` ; `wannier/27x27/wannier_manifest.json` (sha256 de tb / u / u_dis) |
| Maille | a₁ = (2,135490, −1,232926, 0), a₂ = (2,135490, 1,232926, 0) Å → a = 2,4659 Å, \|b\| = 4π/(√3 a) = 2,9423 Å⁻¹ ; c = 15,875 Å | `wannier/27x27/wannier.wout` l. 89–92 |
| Unités | M stocké en Ha, multiplié par 27,211386245988 **une fois** au chargement ; H(R) Wannier et toutes les énergies de la matrice t en eV ; Γ rapporté en eV (×10³ = meV) | `config.py:5`, `matrix_io.py:24, 77–78` ; `config/production.json` clé `units` |
| R_cut | 3, sur la **norme euclidienne des indices réduits** (i² + j² ≤ 9, `np.linalg.norm(R_mwr)`), soit n_L = 29 mailles et dim V_loc = 29 × 5 = 145 | `config/production.json:3` ; `compute_spectral_wannier.py:101` ; `resonance_metrics.py:36` ; log `resonance_9x9_20421370.out` « [cluster] R_cut=3: nL=29 sites, dim=145 » |
| Grille de sortie N_k^out | 240×240 = 57 600 k (Γ-centrée, non décalée, K = (2/3, 1/3) sur la grille : indice 38480) × 5 bandes = 288 000 états ; **41 266 états** dans la fenêtre ±3 eV | `config/production.json:4` ; `local_tmatrix.py:24–28` ; `resonance_metrics.py:49, 51` ; log « 41266 on-shell states in +-3.0 eV on 240x240 » ; « k_out[38480]=[0.667 0.333 0.] » |
| **Grille interne N_k^int (pour g₀)** | **300×300 = 90 000 k**, Δk = \|b\|/300 = 0,00981 Å⁻¹ ; découplée de la grille de sortie ; valeur figée ET valeur par défaut codée en dur (`k_int=None` → 300) ; **jamais variée** : les 53 occurrences dans `results/M/logs/*.out`, `scripts/*.py`, `scripts/*.sh` valent 300 (grep 2026-09-18) ; balayée depuis : §6 (P13 en v1, R6 en v2) | `config/production.json:13` ; `compute_spectral_wannier.py:49, 89` ; `resonance_metrics.py:22, 46` ; `resonance_criteria.py:23, 32` ; `rcut_resigma.py:22, 36` ; défaut : `local_tmatrix.py:171–172, 240–241` |
| Échantillonnage de la grille interne près de E_D (compté ce jour sur les bandes Wannier 27×27) | états à ‖ε − E_D‖ ≤ η = 0,02 eV : **4** (les deux points K × π, π*) ; ≤ 0,1 eV : 52 ; ≤ 0,5 eV : 1 300 ; ≤ 3 eV : 64 528. Pour comparaison : 240² → 4 / 28 / 832 / 41 266 ; 600² → 4 / 220 / 5 140 / 258 196 | calcul direct (`lt.Hwr_to_Hwk` sur `mp_grid(n, n, 1)`, `wannier/27x27/wannier_tb.dat`), non versionné |
| η | 0,02 eV (lorentzienne ε + iη dans g₀ ; même η pour la moyenne lorentzienne des courbes et pour ρ₀) | `config/production.json:5` ; `local_tmatrix.py:150` ; `resonance_metrics.py:75` |
| Fenêtre d'énergie | ±3,0 eV autour de E_D (états sur couche retenus) ; grille d'énergie de g₀ : [E_D − 3 − η, E_D + 3 + η] | `config/production.json:14` ; `resonance_metrics.py:51, 55` ; `local_tmatrix.py:258` |
| Pas de la grille d'énergie | η/ne_per_eta = 0,02/8 = **2,5 meV**, 2 417 énergies (g₀ et dg₀/dε) ; courbes énergétiques sur `egrid[::2]` = 1 209 points (5 meV) | `config/production.json:15` ; `resonance_metrics.py:55, 74` ; log « [g0] 2417 energies, spacing 2.50 meV » |
| E_D (point de Dirac, Wannier) | **−4,238896 eV** (gap minimal 0,0 meV sur 90×90 : K est sur la grille) ; égal à `efermi_read` d'EPW (chapitre 5, −4,2389 eV) | `resonance_metrics.py:47–48` ; log « [dirac] E_D = -4.2389 eV » ; npz `E_D` |
| Grille de ρ₀ | 900×900 (`--rho0-grid`), ρ₀ par cellule et par spin ; `rho0_240` (sur la grille de sortie) aussi stocké | `resonance_metrics.py:18, 80–82` |
| Concentration | c = 1 % par cellule pour ρ_dis = ρ₀ + c δρ ; c = 0,1 % pour la comparaison Kaasbjerg | `config/production.json` clé `defect_concentration_for_dos` ; `resonance_metrics.py:83` ; `resonance_criteria.py:19` |
| Born | t_B = V + V g₀ V (deux premiers termes) → Γ_Born = −2 Im ⟨V g₀ V⟩ (le terme du 1er ordre est réel) | `resonance_metrics.py:60` |
| Alignement du potentiel (ΔV_PA^(N) = C_N) | **Alignement de Kumagai–Oba** (Kumagai et Oba, PRB 89, 195205, 2014 ; région élargie, `campagnes/M/ch4/alignement_regions.md`) : C_N = moyenne des potentiels de site ΔV (sphères de 1,0 Å, `alignment.atom_sphere_shifts`) sur les atomes à distance vraie ≥ 0,75 r_max de la lacune (image minimale) ; rms de ces potentiels = incertitude ; **9×9 : C_9 = −25,1437 meV, rms 9,05 meV, 53 atomes** (13 tailles : 5×5 −57,5626 ; 6×6 −50,5091 ; 7×7 −26,9143 ; 8×8 −14,4103 ; 9×9 −25,1437 ; 10×10 −9,3417 ; 11×11 −9,5507 ; 12×12 −18,6896 ; 15×15 −16,2185 ; 18×18 −15,0621 ; 21×21 −14,2343 ; 24×24 −13,4815 ; 27×27 −13,0058 meV). Application en base de Wannier : M_W(R, R) − C_N·𝕀₅ sur les N² mailles de la boîte (approximation (i), exacte à 0,06 meV à 9×9, R9 A.2 ; `local_tmatrix.defect_mwr`) ; étiquettes hors grille par images de Wigner-Seitz (`ws_images`, `ws_phase`). **La production n'a jamais soustrait la moyenne de ΔV** (`subtract_mean=False`, `scripts/m/compute_M.py:118`). C14 = ±rms (±9,05 meV) autour du V_loc aligné (§3) [v1 : décalage G ≈ 0 = moyenne diagonale de M^L, 67,0 meV, Γ_T recalculé avec M − ⟨M^L⟩·1] | `config/production.json` bloc `alignment` (`C_N_eV`, `source` = `campagnes/R/R10_plateau/a/a1_results.json`, md5 3d98bbd2…) ; `defects/alignment.py` (`atom_sphere_shifts`, `true_min_image_dist`) ; `local_tmatrix.py` (`defect_mwr`) ; `wannier_interpolation.py` (`ws_images`, `ws_phase`, `Mwr_to_Mwk(ws=)`) ; ligne « [align] » des journaux `results/M2_plateau/logs/` ; R10 A.1 |
| δρ | δρ(ε) = (1/π) Im Tr[t(ε) dg₀/dε] par défaut et par spin ; contrôle Lloyd δρ = −(1/π) d/dε Im ln det[1 − g₀V] | `resonance_metrics.py:62` ; `resonance_criteria.py:57–66` |
| T̄ en K | paire π (bandes 3, 4 de 5), trace/2 de ⟨n K\| t \|n′ K⟩ | `resonance_metrics.py:86–90` |
| Test d'or | η = 0,10 eV, k_int = k_out = grille grossière, R_local = boîte MP-duale complète (exact) ; seuil rel < 1e-8 et Γ ≥ −1e-8 | `scripts/validation/test_local_tmatrix_real.py:20, 53, 59, 61` |
| Coût (16 cœurs) | final (R10 C) : cartes niveau 1 `specwd_<S>` 1 h 17 / 2 h 01 / 2 h 47 / 2 h 38 / **1 h 01** / 1 h 11 (5, 6, 7, 8, 9, 12 ; R_cut 0–4) ; `post_res9` (métriques + critères 9×9) 39 min 45 ; `post_res6` 6 h 25 ; `post_res12` 2 h 20 ; `nkint` 11 min 33 ; `resigma` 2 min 47 + 8 min 14 ; C14 3 min 36 ; test d'or 5×5 dense non rejoué (D12 : recopié de R6, 1 h 52) [v1 : carte (36 combinaisons) 1 h 15 ; 7 min 39 ; 1 h 34 ; test ne_per_eta 3 min 35] | R10_rapport.md, table « Rapport — C » (jobs 22058839–44, 22058850–52, 22058845–47, 22058854) [v1 : 20294202, 20421370, 20414406, 20684497] |

## 2. Valeurs de production (9×9, R_cut 3, 240², η 0,02, N_k^int 300²) — final (R10, `results/M2_plateau/`)

| quantité | valeur | source |
|---|---|---|
| médiane de \|Γ\|·N_cells sur les 41 266 états (carte niveau 1) | **3 189,01 meV** (Γ par défaut, ×81 = intensif) [v1 : 2 473,55] | `results/M2_plateau/level1_summary.csv` l. « 9x9,3,240,0.02 » (3189.0136) ; `specwd_9x9_prod.npz` (job 22058843, 2026-09-29) |
| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 3 771,46 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 1,7944e+05 meV ; Born/T médian 45,846 (min 1,050, max 202,851 sur ±3 eV, lus sur les courbes `Gamma_Born`/`Gamma_T` du npz) [v1 : 2 343,59 ; 7 405,83 ; 3,302 (0,677, 16,45)] | `campagnes/R/R10_plateau/c/table_v2_plateau.md` l.101–103 (`resonance_9x9.npz`) ; min/max : `campagnes/M/ch4/table_v1_final.md`, compléments 5.4 (porte : relecture v1 et non alignée = valeurs publiées) |
| Re Σ médian, R_cut 3 | **+605,08 meV** ; \|Re Σ\|/Γ médian 0,2173 [v1 : 722,3 ; 0,227] (valeur non alignée : §8) | `results/M2_plateau/m_rcut_resigma.csv` (med_ReSigma_meV 605.0832…, med_absReSigma_over_Gamma 0.21732…) |
| position de résonance E_res − E_D (argmax des états à ±1,5 eV) — **retiré du ch. 4 (R10 B.1 : argmax discret sur une couronne de la grille de sortie ; gardé ici pour la traçabilité)** | −0,175 eV (R_cut 2, 3, 4 à η 0,01/0,02 ; R_cut 0 à η 0,01/0,02 ; −0,181 à R_cut 0 et 4, η 0,05) ; −0,227 (R_cut 1) ; contre la grille de sortie (N_k^int 900) : −0,134 / −0,175 / −0,175 / −0,172 eV à 120² / 240² / 480² / 960² (non aligné : −0,227 / −0,202 / −0,200 / −0,191) [v1 : −1,238 ; −1,183 (η 0,01) ; −1,292 (R_cut 2 et 4)] | `results/M2_plateau/level1_summary.csv` (argmax_E_minus_ED_eV) ; `campagnes/R/R10_plateau/b/b_results.json` B1 |
| pics des courbes (ε − E_D) | Γ_T : −0,180 eV ; Γ_Born : +1,695 eV ; Γ_T/ρ₀ : −0,170 eV ; δρ : −0,787 eV ; ρ_dis : +1,710 eV ; \|T̄\| : −0,170 ; −Im T̄ : −0,172 ; min \|Re T̄\| : −2,140 ; Re T̄(E_D) = 12,624 eV, Im T̄(E_D) = −2,119 eV ; zéros de Re T̄ : −2,1425, −2,130, −2,095, −2,070, −2,0475, −2,005, −2,000, −0,210 [v1 : −1,24 ; +1,695 ; −0,015 ; −2,53 ; −2,53 ; −0,905 ; −1,29 ; −2,145 ; 2,521 / −0,091 ; aucun zéro] | `results/M2_plateau/resonance_9x9.npz` (`peak_*`, `ReTbar_at_ED`, `ImTbar_at_ED`, `Tbar_zero_crossings`) |
| critère det / valeur propre | matrice complète (dim 145) : min de \|det[1 − Vg₀]\|/max = 1,261e-4 à −0,785 eV ; min_i \|λ_i\| = 0,0019 à −0,787 eV. Bloc σ (dim 87) : 4,379e-4 à −0,7875 eV, \|λ\| = 0,00186. Bloc π (dim 58) : aucun zéro, min \|det\|/max 1,706e-2 à −0,170 eV, min \|λ\| = 0,3440 à −0,127 eV [v1 : 2,09e-4 à −2,530 eV ; 0,0108 (λ = +0,0003 + 0,0108 i) ; secondaires −2,19 … −1,97 eV] | `campagnes/R/R10_plateau/c/table_v2_plateau.md` l.115–117 ; `results/M2_plateau/resonance_criteria_9x9.npz` (`sigma_flag_at`, `sigma_flag_minlam`, `sigma_flag_det_rel`, `dim_pi`, `dim_sigma`) |
| règle de somme de Friedel | ∫δρ sur toute la bande = **−1,0005** états (Tr[t g₀′]) et −1,0005 (Lloyd) ; bloc π −0,9981, bloc σ −0,0024 ; dans ±3 eV : +0,702 (π −0,256, σ +0,958) ; écart ponctuel max et cumuls aux bords non publiés par R10 [v1 : −0,0569 ; +1,782 ; 0,614 ; −2,190 / −0,409] | `results/M2_plateau/resonance_criteria_9x9.npz` (`sumrule`, `sumrule_lloyd`, `sumrule_window`, `*_pi`, `*_sigma`) |
| Γ_T à c = 0,1 % sur ±1 eV | min 2,48 meV (+1,00 eV), max 40,58 meV (−0,18 eV), 8,47 meV à E_D ; ħ/Γ à ∓0,3 eV : 26 / 212 fs (lus sur la courbe `x_c`, `Gamma_c` du npz, c_compare 0,001, mêmes définitions que R6) [v1 : 0,63 (+0,24) ; 5,55 (−1,00) ; 1,27 ; 419 / 1 025 fs] | `results/M2_plateau/resonance_criteria_9x9.npz` ; `campagnes/M/ch4/table_v1_final.md`, compléments 5.4 (porte : relecture v1 et non alignée = valeurs publiées) |
| localité de M_W (9×9 dense) | ‖M_W(0, 0)‖ = 41,927 eV ; p_z–p_z sur le site de la lacune 31,546 eV, p_z de l'autre sous-réseau 0,640 eV ; ‖M_W(R, 0)‖ hors site inchangés par l'alignement (seule la diagonale sur site change : \|R\| = a : 1,963 / 1,388 / 1,501 ; √3 a : 0,137 / 0,639 / 0,122 …) ; abscisse maximale 15,588 a (images de Wigner-Seitz, contre 22,517 a avec les étiquettes brutes) [v1 : 9,365 ; 6,617 ; 0,044 ; a : 0,658 / 0,460 / 0,437 ; √3 a : 0,165 / 0,079 / 0,035 / 0,029 ; 2a : 0,082 / 0,081 — le texte du 2026-09-18 rangeait 0,460 dans la 2ᵉ couronne] | `results/M2_plateau/mwr_locality.npz` (`9x9_dense_onsite_norm`, `9x9_dense_onsite_pzvac`, `9x9_dense_onsite_pzB`, `9x9_dense_dist`, `9x9_dense_w`) ; table R10 l.72–74 |
| recentrage | R_d = [4, 4, 0] sur la boîte 27×27 (cohérent avec s_red = 13/27 = 4·3 + 1, p = 3) ; inchangé ; C_N soustrait sur les 81 mailles de la boîte (étiquettes brutes de `Mwk_to_Mwr`) puis recentrage (`defect_mwr`) | lignes « [align] » et « [recenter] » des journaux `results/M2_plateau/logs/` ; R10 C.1 |

## 3. Contrôles réellement implémentés (à confronter aux dix du mémoire)

Type : **P** = porte bloquante (le code refuse, `raise`), **T** = test PASS/FAIL, **C** = convergence (paramètre varié),
**K** = cohérence physique (calculée et rapportée, sans seuil).

| n° | contrôle | type | implémentation | seuil | verdict (date, preuve) |
|---|---|---|---|---|---|
| C1 | Provenance de jauge : tb / u / u_dis d'un même run Wannier90, sha256 recalculés | P | `wannier_provenance.py:42–68`, appelé `compute_spectral_wannier.py:52`, `resonance_metrics.py:21`, `resonance_criteria.py:20` | tout écart ⇒ `ValueError` | « [gauge] provenance OK » dans tous les logs de production (manifeste 27×27) |
| C2 | Normalisation de Bloch et unités de M : sidecar `bloch_norm = unit_cell` (potentiel intensif) et `units = hartree`, conversion eV unique ; depuis R6, `M_normalization = v2` (M^L et M^NL en norme unit_cell) | P | `matrix_io.py:50–78` ; contrat commenté `compute_spectral_wannier.py:60–66` ; R6 : `load_M_checked(…, require_normalization=M_NORM_V2)` et `matrix_io.check_manifest` avant les lectures mmap | sidecar absent / mauvais tag / M v1 ⇒ `ValueError` | en vigueur depuis f876605 (2026-09-07) ; porte v2 depuis R6 (2026-09-26) : `results/M2/M_dense_9x9.json` conforme, `results/M/M_dense_9x9.npy` refusé |
| C3 | Cohérence des grilles k : U/U_dis réordonnés sur la grille de M, point sans correspondant ⇒ erreur ; grille MP complète non décalée | P | `wannier_interpolation.py:142–183` (`_match_kpoint_order`, `_infer_mp_grid`) | tol 1e-5 ; N₁N₂N₃ ≠ nk ⇒ `ValueError` | implicite dans chaque run |
| C4 | Garde-fou a : défaut centré (‖M_W(R,0)‖ maximal en R = 0) | P | `local_tmatrix.py:36–51` ; `compute_spectral_wannier.py:83` ; `resonance_metrics.py:35` | argmax ≠ R₀ ⇒ `AssertionError` | passé (valeurs §2) |
| C5 | Garde-fou b : hermiticité de V_loc | P (symétrise + journalise) | `local_tmatrix.py:54–76` | résidu > 1e-10 ⇒ message | v2 : 1,3e-14 (golden 5×5 dense M2, 2026-09-26) [v1 : 2,3e-15 – 4,0e-15, logs golden 2026-09-05/07] |
| C6 | **Test d'or** : matrice t locale = matrice T dense (`compute_T`) sur la même grille grossière, même sous-espace 5 WF | T, bloquant avant production | `scripts/validation/test_local_tmatrix_real.py` ; `scripts/validation/test_local_tmatrix.py` (synthétique) ; `scripts/validation/test_local_rcut.py` (support tronqué ⇒ écart) | rel < 1e-8 et Γ ≥ −1e-8 | v2 : PASS 5×5 dense 1,80e-13 (M2, job 21852238, 2026-09-26, `results/M2/logs/r6golden_21852238.out`) ; 6×6 et 12×12 non rejoués avec M2 [v1 : 5×5 dense 2,3e-14 (2026-09-05, après le correctif d'unités), 6×6 6,06e-9, 12×12 3,68e-9 (2026-09-07) ; logs `golden_dense_20294198/20444392/20435998.out`] |
| C7 | Positivité Γ_nk ≥ 0 | P | `local_tmatrix.py:269–271` ; `rcut_resigma.py:55` | min Γ < −1e-8 ⇒ `AssertionError` | v2 : min Γ_loc = +0,2646 eV (golden 5×5 dense M2) ; jamais déclenché [v1 : +0,26 eV] |
| C8 | g₀ batché = g₀ de référence (restructuration exacte) et évaluation « point d'énergie le plus proche » vs état par état | T | `scripts/validation/test_local_green_batch.py:15` | rel < 1e-12 | ≤ 1,21e-14 pour R_cut 0–3 (rejoué 2026-09-18, grille 90², 7 énergies ; ne dépend pas de M) ; nearest-grid vs exact : 5,4e-3 avec ne_per_eta 8 (cas synthétique) |
| C9 | Pas de la grille d'énergie ne_per_eta | C | `rcut_resigma.py --npe` (9e0d1a9, 2026-09-09) | — | v1 seulement (non rejoué avec M2 ni avec l'alignement) : ne_per_eta 2 → 32 : médiane Γ 2 476,19 → 2 474,33 meV (0,08 %), E_res inchangé (E_res retiré, R10 B.1) ; log `npe_test_20684497.out` |
| C10 | R_cut (support de V_loc) | C | carte niveau 1 R_cut 0–3 (0–4 depuis R6) ; `rcut_resigma.py` R_cut 0–4 | plateau ≤ 5 % (énoncé `compute_spectral_wannier.py:118`) | final, 9×9 : 3 013,96 / 3 273,25 / 3 222,08 / **3 189,01** / 3 200,00 meV (R_cut 0…4) ; écart médian par état à R_cut 4 : 9,79 / 4,52 / 1,98 / **0,66 %** ; Re Σ médian 29,5 / 761,8 / 737,2 / 605,1 / 696,7 meV (`results/M2_plateau/m_rcut_resigma.csv`, `rel_med_dGamma`, `med_ReSigma_meV`) [v1 : 2 596,9 / 2 488,9 / 2 466,8 / 2 473,5 / 2 468,0 ; 10,0 / 4,6 / 1,5 / 0,63 % ; Re Σ 515 / 694 / 723 / 722 / 720] |
| C11 | Grille de sortie × η (plateau conjoint) | C | carte niveau 1 (grilles 60/120/240, η 0,05/0,02/0,01) | ≤ 5 % quand η/2 et grille ×2 | final, 9×9, R_cut 3 (`results/M2_plateau/level1_summary.csv`) : 120² → 240² : 3 200,71 → 3 189,01 (η 0,02 ; 0,37 %, arithmétique), 3 208,42 → 3 202,75 (η 0,01 ; 0,18 %) ; η 0,02 → 0,01 à 240² : 3 189,01 → 3 202,75 (0,43 %) ; 0,05 → 0,02 : 3 232,61 → 3 189,01 (1,37 %) [v1 : 0,07 % ; 0,3 % ; 0,001 % ; 0,74 %] |
| C12 | Taille de super-cellule N (niveau 2, familles N mod 3) | C | `scripts/t/level2_families.py`, `results/M2_plateau/level2_summary.csv`, `level2_families.csv` | N ≥ 7 (config l. 6) | final, 5/6/7/8/9/12 : 2 997,14 / 3 149,15 / 2 942,85 / 2 974,18 / 3 189,01 / 3 264,51 meV ; (max − min)/moyenne, arithmétique (= R9 clôture R.4) : famille 3m (6, 9, 12) 3,60 %, non-3m (5, 7, 8) 1,83 %, 7–9 : 8,11 % [v1 : 2 524,3 / 2 509,2 / 2 487,5 / 2 478,5 / 2 473,5 / 2 458,7 ; 0,56 % ; 2,0 % ; 1,84 %] |
| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | final : Born/T médian 45,846 sur ±3 eV (min 1,050, max 202,851 ; compléments 5.4) [v1 : 3,30] |
| C14 | Sensibilité à l'alignement du potentiel : **C = ±rms du plateau de Kumagai–Oba (±9,05 meV à 9×9) ajouté uniformément au V_loc aligné** (R10 D10) | K | `resonance_metrics.py --shift-L-meV 9.05,-9.05` (V_loc + C·1 sur la boîte) | — | final : médiane des états 3 188,35 → 3 237,29 (+9,05 meV) / 3 157,07 (−9,05 meV) ; pic de la courbe Γ_T −0,180 eV inchangé ; écart relatif de la courbe Γ_T : max 6,36e-2 / 6,10e-2, médian 1,37e-2 / 1,24e-2 (`results/M2_plateau/resonance_9x9_shiftL.npz` ; R10 C.1). Retirés : la variante « M − ⟨M^L⟩·1 » (décalage 5 427,3 meV en v2) et C = ±25 meV autour du V_loc non aligné (R6 3.5) [v1 : décalage 67,0 meV ⇒ 6,4e-4 (médian 1,8e-4)] |
| C15 | Règle de somme de Friedel, deux formules (Tr[t g₀′] vs Lloyd) | K | `resonance_criteria.py:57–75` ; R6 : par bloc (`--blocks full,pi,sigma`) | — | final : −1,0005 / −1,0005 états sur toute la bande (π −0,9981, σ −0,0024) ; +0,702 dans ±3 eV [v1 : −0,0569 / −0,0569 ; +1,78] |
| C16 | Critère de résonance (det, valeur propre minimale) et position du pic | K | `resonance_criteria.py:41–55` ; R6 : par bloc, `--flag-eV` | — | final : minimum global à −0,785 eV (det) / −0,787 eV (\|λ\| = 0,0019), porté par le bloc σ ; bloc π sans zéro (min \|λ\| 0,3440 à −0,127 eV) ; pic de la courbe Γ_T −0,180 eV ; E_res(argmax Γ) retiré (R10 B.1 ; valeur −0,175 eV) [v1 : minimum unique à −2,530 eV ; E_res −1,24 eV] |
| C17 | K sur la grille de sortie et dégénérescence π/π* en K | K | `resonance_metrics.py:86` ; `config/production.json` clé `K_red` | — | k_out[38480] = (2/3, 1/3), E(π) = E(π*) = E_D |
| C18 | **Porte A.2 (R6)** : ΔV^L (grille) et ΔV^NL (projecteur KB de l'atome retiré) appliqués directement aux états de Bloch purs repliés, comparés à M2/N_cells | P | `scripts/m/gate_M_normalization.py` (`defects/deltav_pw.py`, `wavefunctions/sc_projection.py`) | max \|écart\| ≤ 1e-6 eV, sinon refus | 14 M2 OK (8 grossiers ≤ 3,3e-14 eV, 6 denses ≤ 1,2e-8 eV ; job 21820491, 2026-09-25) ; M v1 refusé (M^L = N_cells × direct, rapport 81,000000) |
| — | **Grille interne N_k^int** | C (depuis P13) | `rcut_resigma.py --nk-int` ; `submit_nkint_check.sh` ; `nkint_check_post.py` | 5 % | final, 300 → 600 : médiane −0,09 %, Re Σ médian −0,03 %, Γ_T(E_D) +4,82 % (arithmétique sur `results/M2_plateau/nkint_check_9x9.csv`) [v1 : −0,14 % ; +5,36 %] ; §6e |

## 4. Points ouverts constatés le 2026-09-18

(État v1 du 2026-09-18, inchangé ; valeurs finales aux §2, §3 et §6e, points R6 au §7, R10 au §8. Le point 4 (E_res) est clos par le retrait d'E_res du ch. 4, R10 B.1.)

1. **N_k^int n'a jamais été varié.** Toute la production (cartes niveau 1, niveau 2, résonance, critères, Re Σ, test
ne_per_eta) utilise 300² = 90 000 points ; c'est aussi le défaut codé en dur de `scattering_rate` /
`scattering_rate_fast`. Aucun log ni npz ne porte une autre valeur. Le plateau grille × η (C11) porte sur la grille de
**sortie** et ne teste pas la somme interne de g₀. Près de E_D, la grille interne ne contient que 4 états dans ±η (les
deux K), 52 dans ±0,1 eV ; le pic de Γ_T/ρ₀ est justement à −0,015 eV. Pour chiffrer : `scripts/t/rcut_resigma.py`
accepte désormais `--nk-int` (même mécanisme que `--npe`, valeur inscrite dans le npz) et
`scripts/slurm/submit_nkint_check.sh 9x9 "150 300 450 600"` enchaîne les quatre valeurs à R_cut 3 (≈ 4 × (1 → 4) × 1 min de g₀,
3 h demandées, 64 G). **Fait le 2026-09-18 (P13, job 21337627) : voir §6.**
2. **R_cut est une norme sur les indices réduits, pas une distance cartésienne.** Le disque i² + j² ≤ 9 (29 mailles) va
jusqu'à 3a le long de a₁ ou a₂ mais jusqu'à 2√3 a ≈ 3,46a le long de a₁ + a₂, et exclut (3, −3) dont |R| = 3a. Si le
mémoire écrit « R_cut = 3 mailles » comme un rayon, préciser la norme (ou l'énoncer en nombre de mailles : 1 / 5 / 13 / 29 / 49
pour R_cut 0…4). Le tracé `fig_locality_final` utilise, lui, la distance cartésienne |R| en unités de a
(`mwr_locality.npz`, distances 1, √3, 2, …).
3. **Deux « médianes de Γ » coexistent** : la médiane sur les états (2 473,5 meV, carte niveau 1, `level1_summary.csv`,
fig_convergence / fig_level2) et la médiane de la courbe lorentzienne Γ_T(ε) (2 343,6 meV, log de `resonance_metrics`).
Ne citer que la première comme « Γ N_cells » ; la seconde n'est qu'un diagnostic Born/T.
4. **E_res est discret** : argmax sur des états, fenêtre ±1,5 eV codée en dur (`compute_spectral_wannier.py:112`,
`rcut_resigma.py:56`) ; il saute entre −1,18, −1,24 et −1,29 eV selon (R_cut, η). Le critère det (C16) donne −2,53 eV,
ce n'est pas la même quantité (pôle de [1 − Vg₀]⁻¹ vs maximum de Γ sur couche).
5. **Test d'or 6×6 à 6,1e-9** pour un seuil de 1e-8 (5×5 : 2e-14, 12×12 : 4e-9) : PASS, mais pas « ~1e-10 » comme
l'annonce le message du script.
6. **Règle de somme** : les intégrales Tr[t g₀′] et Lloyd coïncident (−0,0569) alors que les intégrandes diffèrent
ponctuellement jusqu'à 0,61 états/eV ; la compensation de l'excès +1,78 de la fenêtre ±3 eV se fait surtout dans
[−3,4, −2,9] eV (−0,79) et près du bas de bande ([−21, −18,6] eV : +0,94 / −0,49).
7. **Contrôle « k coïncidents » (M dense vs M grossier)** : 6×6 a rendu « CHECK (nscf gauge/convergence tolerance
exceeded) », écart de valeurs singulières 6,4e-2 (log `golden_dense_20444392.out`). C'est un contrôle de M (chapitre 4,
§ matrice), pas de t/Γ, mais il figure dans la même chaîne de jobs.

## 5. Fichiers touchés ce jour (non commis)

Test d'or 5×5 grossier rejoué le 2026-09-18 (nœud de connexion, 1,3 s) après correction du message (P14 ; seuil réel 1e-8 au
lieu de « ~1e-10 ») — sortie :

```
[5x5] REAL GOLDEN: max|local - dense*N_cells| rel = 5.10e-14 (seuil 1e-8) : PASS
[5x5] positivity min Gamma(local) = 6.596e-01
RESULT: PASS
```
(`scripts/validation/test_local_tmatrix_real.py` ; `test_local_rcut.py` corrigé de même, « must be ~1e-8 » → « seuil 1e-8 » ; `test_local_tmatrix.py` n'avait pas de message trompeur.)

- `scripts/t/rcut_resigma.py` : option `--nk-int` (2 lignes, même patron que `--npe` ; défaut = config figée, valeur
  inscrite dans le npz et le log).
- `scripts/slurm/submit_nkint_check.sh` : job de balayage N_k^int (lancé : job 21337627, §6).
- `scripts/t/nkint_check_post.py` : post-traitement du balayage (tableaux de §6, `results/M/nkint_check_9x9.csv`).
- `NOTES_TGAMMA.md` : ce document.

## 6. Balayage de la grille interne N_k^int (P13, 2026-09-18)

Job Slurm 21337627 (`scripts/slurm/submit_nkint_check.sh 9x9 "150 300 450 600"`, rc32615, 16 cœurs, 64 G demandés, MaxRSS 14906672K,
départ 08:27:19, fin 09:04:59 EDT, 37 min 40 ; 8 à 10 min par valeur, dominées par le chargement de M et la rotation, pas par g₀).
Chaîne identique à la production : `scripts/t/rcut_resigma.py --nk-int N` (option ajoutée ce jour, l. 20 et 22 ; `nk_int` inscrit
dans le npz l. 41 et dans le log l. 34), R_cut 3, grille de sortie 240², η 0,02, fenêtre ±3 eV, ne_per_eta 8, portes C1–C7
actives (« [gauge] provenance OK », recentrage R_d = [4, 4, 0], positivité l. 55). `config/production.json` inchangé (nk_int 300).
Dry run préalable (nœud de connexion, nk_int 60, grille 30, 17 s) : valeur présente dans le log et le npz.
Sources : `results/M/resigma_9x9_rc3_nk{150,300,450,600}.npz` (clés `Sigma_rc3`, `E_out`, `E_D`, `nk_int` ; non commis),
log `results/M/logs/nkint_21337627.out`, post-traitement `scripts/t/nkint_check_post.py` → `results/M/nkint_check_9x9.csv` (commis).
Définitions : médiane **sur les états** de |Γ| (Γ = −2 Im Σ_nk, V_loc intensif ⇒ Γ N_cells ; §4.3), Re Σ médian sur les mêmes
états, E_res = argmax de Γ sur les états à |ε − E_D| ≤ 1,5 eV (`rcut_resigma.py:56`), Γ_T(E_D) = moyenne des 4 états dont
ε = E_D exactement (les deux K × π, π* ; il n'y a pas d'autre état dans ±η sur 240²). Écarts relatifs à nk_int = 600.

### 6a. États de la fenêtre ±3 eV (41 266 états)

| nk_int | N_k^int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) | Γ_T(E_D) (meV, 4 états à K) | écart |
|---|---|---|---|---|---|---|---|---|
| 150 | 22 500 | 2 398,96 | −3,15 % | 725,37 | +0,92 % | −1,116 | 248,81 | +44,06 % |
| 300 (production) | 90 000 | 2 473,55 | −0,14 % | 722,28 | +0,49 % | −1,238 | 181,96 | +5,36 % |
| 450 | 202 500 | 2 478,71 | +0,07 % | 716,87 | −0,27 % | −1,238 | 174,15 | +0,83 % |
| 600 | 360 000 | 2 476,95 | 0 | 718,79 | 0 | −1,245 | 172,71 | 0 |

Le run nk 300 reproduit la production à l'identique (2 473,55 meV, −1,238 eV ; `m_rcut_resigma.csv`, `level1_summary.csv`).
Détail des 4 états à K (Γ en meV, nk 150 / 300 / 450 / 600) : 246,8 – 250,8 / 180,5 – 183,4 / 172,7 – 175,6 / 171,3 – 174,1 ;
Re Σ à K : 2 498 – 2 539 / 2 504 – 2 545 / 2 505 – 2 546 / 2 505 – 2 546 meV.

### 6b. États à |ε − E_D| ≤ 0,3 eV (280 états, zone à faible échantillonnage interne)

| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) | min – max de Γ (meV) |
|---|---|---|---|---|---|---|
| 150 | 439,82 | −1,13 % | 2 329,26 | −0,49 % | −0,284 | 170,6 – 1 522,0 |
| 300 (production) | 433,99 | −2,44 % | 2 347,61 | +0,29 % | −0,284 | 180,5 – 1 228,9 |
| 450 | 443,93 | −0,21 % | 2 341,34 | +0,02 % | −0,284 | 172,7 – 1 212,4 |
| 600 | 444,85 | 0 | 2 340,77 | 0 | −0,284 | 171,3 – 1 212,5 |

### 6c. Constats (seuil de C10/C11 : 5 %)

- Médiane Γ N_cells sur la fenêtre : 300 → 600 = −0,14 % ; Re Σ médian : +0,49 % ; zone ±0,3 eV : −2,44 % (Γ) et +0,29 % (Re Σ).
  Tous sous 5 %. 150 → 600 : −3,15 % sur la fenêtre.
- **Γ_T(E_D) (états à K) : 300 → 600 = +5,36 %, au-dessus du seuil de 5 %** ; 450 → 600 = +0,83 % ; 150 → 600 = +44 %.
  La valeur à E_D décroît de façon monotone avec N_k^int (248,8 → 182,0 → 174,2 → 172,7 meV).
- **E_res change de point** : −1,116 eV (150), −1,238 eV (300 et 450), −1,245 eV (600) — écart de 7 meV entre 450 et 600, soit
  un état voisin sur la grille 240² (argmax discret, §4.4). Dans la zone ±0,3 eV, E_res = −0,284 eV pour les quatre valeurs.
- Le maximum de Γ dans ±0,3 eV passe de 1 522 (150) à 1 229 (300) puis 1 212 meV (450, 600).

### 6d. Même balayage avec M2 (R6, job 21857278, 2026-09-26, 7 min)

`scripts/slurm/submit_nkint_check.sh 9x9`, mêmes options ; sources `results/M2/resigma_9x9_rc3_nk{150,300,450,600}.npz`, post-traitement
`nkint_check_post.py` (job post_fig 21872955) → `results/M2/nkint_check_9x9.csv`. Écarts relatifs à nk_int = 600.

| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) | Γ_T(E_D) (meV, 4 états à K) | écart |
|---|---|---|---|---|---|---|---|
| 150 | 3 155,57 | +0,64 % | −23,25 | −7,33 % | −0,238 | 4 801,86 | +38,63 % |
| 300 (production) | 3 132,60 | −0,10 % | −26,65 | +6,23 % | −0,175 | 3 630,65 | +4,81 % |
| 450 | 3 135,87 | +0,01 % | −26,79 | +6,79 % | −0,202 | 3 489,52 | +0,74 % |
| 600 | 3 135,65 | 0 | −25,09 | 0 | −0,202 | 3 463,87 | 0 |

États à |ε − E_D| ≤ 0,3 eV (280 états) :

| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV) | min – max de Γ (meV) |
|---|---|---|---|---|---|---|
| 150 | 3 801,49 | +41,18 % | 5 421,24 | +2,40 % | −0,238 | 1 996,3 – 56 388,2 |
| 300 (production) | 2 749,49 | +2,11 % | 5 269,21 | −0,47 % | −0,175 | 2 508,3 – 49 568,5 |
| 450 | 2 717,85 | +0,93 % | 5 292,56 | −0,03 % | −0,202 | 2 589,5 – 47 666,9 |
| 600 | 2 692,71 | 0 | 5 294,35 | 0 | −0,202 | 2 596,6 – 47 693,8 |

Détail des 4 états à K (nk 150 / 300 / 450 / 600) : Γ 4 759,1 – 4 844,6 / 3 598,2 – 3 663,1 / 3 458,3 – 3 520,7 / 3 432,9 – 3 494,8 meV ;
Re Σ 10 799 – 10 986 / 11 085 – 11 277 / 11 117 – 11 310 / 11 124 – 11 317 meV (mêmes sélections qu'en 6a–6b ; appliquées aux npz v1,
elles redonnent les valeurs v1 ci-dessus).

Constats (seuil 5 %) : médiane de la fenêtre 300 → 600 : −0,10 % ; Γ_T(E_D) 300 → 600 : +4,81 % (v1 : +5,36 %) ; Re Σ médian 300 → 600 :
+6,23 % d'une valeur de −25 meV (écart absolu 1,56 meV) ; zone ±0,3 eV : +2,11 % (Γ), −0,47 % (Re Σ) ; E_res : −0,238 / −0,175 / −0,202 /
−0,202 eV.

### 6e. Même balayage, final (R10 C, job 22058845 `nkint`, 2026-09-29, 11 min 33)

`scripts/slurm/submit_nkint_check.sh 9x9`, mêmes options, V_loc aligné (Kumagai–Oba, C_9 = −25,1437 meV) ; sources `results/M2_plateau/resigma_9x9_rc3_nk{150,300,450,600}.npz`,
post-traitement `nkint_check_post.py` (C1f, job 22058855) → `results/M2_plateau/nkint_check_9x9.csv`. Écarts relatifs à nk_int = 600 (arithmétique sur le csv).

| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV, retiré) | Γ_T(E_D) (meV, 4 états à K) | écart |
|---|---|---|---|---|---|---|---|
| 150 | 3 210,02 | +0,57 % | 607,78 | +0,42 % | −0,238 | 5 662,16 | +38,43 % |
| 300 (production) | 3 189,01 | −0,09 % | 605,08 | −0,03 % | −0,175 | 4 287,11 | +4,82 % |
| 450 | 3 196,11 | +0,13 % | 604,19 | −0,18 % | −0,175 | 4 120,38 | +0,74 % |
| 600 | 3 191,93 | 0 | 605,25 | 0 | −0,175 | 4 090,11 | 0 |

États à |ε − E_D| ≤ 0,3 eV (280 états ; colonnes `z_medG`, `z_medRe`, `z_E_res` du csv) :

| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV, retiré) |
|---|---|---|---|---|---|
| 150 | 4 288,84 | +41,12 % | 6 056,58 | +2,03 % | −0,238 |
| 300 (production) | 3 045,28 | +0,20 % | 5 858,01 | −1,31 % | −0,175 |
| 450 | 3 079,75 | +1,34 % | 5 942,08 | +0,10 % | −0,175 |
| 600 | 3 039,06 | 0 | 5 935,94 | 0 | −0,175 |

Constats (seuil 5 %) : médiane de la fenêtre 300 → 600 : −0,09 % ; Γ_T(E_D) 300 → 600 : +4,82 % (v2 : +4,81 % ; v1 : +5,36 %) ; Re Σ médian 300 → 600 : −0,03 % ;
zone ±0,3 eV : +0,20 % (Γ), −1,31 % (Re Σ). Le run nk 300 reproduit la production finale (3 189,01 meV, `level1_summary.csv`, `m_rcut_resigma.csv`).

## 7. R6 (2026-09-25 → 2026-09-27) : passage à M2

- **Cause** : dans tous les M de production, M^L venait d'états de Bloch normalisés sur la super-cellule et M^NL sur la maille (R5-A.2) :
  M^L était N_cells = N² fois trop petit. Noyau corrigé par Greg (`local_R.py`, commits 99d64da, a223687) ; M^L grossiers 5, 6, 7, 8, 9
  recalculés ; les autres M réassemblés, M2 = N_cells·M^L + M^NL (`scripts/m/assemble_M2.py`, sidecar `assembled_from`, `N_cells`,
  `M_normalization = v2`). 11×11 exclu (`.save` de maille écrasé).
- **Portes** : C18 (porte A.2) sur chaque M2 ; C6 (test d'or 5×5 dense avec M2) ; porte de niveau 1 (R6 §3.3 : même critère que le gel du
  2026-09-05 sur les six tailles) : OK, `config/production.json` inchangé hors clés `M_normalization` et `results_dir`.
- **Code** : `config.results_dir(cfg)` pour tous les chemins ; `matrix_io` : porte v2 et `check_manifest` ; `resonance_criteria.py --blocks
  full,pi,sigma --flag-eV` ; `resonance_metrics.py --shift-L-meV` ; `epw_ed_vs_ep.py` écrit dans `results_dir`.
- **Constats ouverts (non jugés)** :
  1. Test « convention intensive » (tab:tests_M) : 7,7e-2 pour un seuil de 5e-2 (« À VOIR ») [v1 : 2,7e-2] ; max\|M\| des bandes 1–16 :
     23,85 (5×5), 23,82 (7×7), 24,31 (8×8), 25,35 eV (9×9).
  2. Contrôle « k coïncidents » (5×5 dense contre grossier, 16 bandes, job du test d'or) : écart de valeurs singulières 0,348, « CHECK »
     (v1 6×6 : 0,064) ; la version de tab:tests_M (bandes 1–8 / 1–15, 9×9) : 5,0e-4 / 1,3e-3, OK.
  3. C6 : 6×6 et 12×12 non rejoués avec M2 ; C9 non rejoué.
  4. C10 : la porte de R6 §3.3 a comparé les médianes (R_cut 3 → 4 : 0,47 %) ; l'écart médian par état (définition de C10) vaut 1,15 %.
  5. `results/M2/m_rcut_resigma.csv` reconstruit par `R6_production_corrigee/r6_m_rcut_resigma.py` (le script de ff39c7f n'est pas dans
     `scripts/` ; appliqué aux npz v1, il redonne le csv v1 à 0,0 près) ; `rel_med_dReSigma` a un dénominateur signé (définition v1),
     négatif en v2 (Re Σ médian à R_cut 4 : −410,3 meV).
  6. Points du §4 à relire avec les valeurs v2 : §4.3 (médiane des états 3 132,6 meV, de la courbe 3 740,2 meV), §4.4 (E_res −0,175 /
     −0,202 / −0,227 eV ; critère det −0,812 eV), §4.6 (Friedel −1,0007 ; +0,669 dans ±3 eV).


## 8. R10 (2026-09-29 → 2026-09-30) : alignement de Kumagai–Oba, production finale `results/M2_plateau/`

- **Décision** (Greg, R10 étape G) : le potentiel de défaut est aligné par l'alignement de Kumagai–Oba (§1, ligne « Alignement du potentiel » ;
  région, sphères de 1,0 Å, rms = incertitude ; 13 tailles, `config/production.json` bloc `alignment`) ; l'alignement à site unique est abandonné
  (les six tailles 5…12 étaient « sans plateau » au critère de R9 A.1, max|écart| 14–234 meV ; R10 D14 : aucun critère d'arrêt). E_res est retiré
  du ch. 4 (R10 B.1 : argmax discret, saute entre couronnes de la grille de sortie, −0,227/−0,202/−0,200/−0,191 eV non aligné et
  −0,134/−0,175/−0,175/−0,172 eV final à 120²…960², `campagnes/R/R10_plateau/b/b_results.json`).
- **Code** (commits c61c118, ea91d61, c2bc733 ; ce dernier était `23ee3bb` avant la réécriture de l'historique du 2026-09-29, nom gardé dans les archives R10) : `wannier_interpolation.ws_images`, `ws_phase`, `Mwr_to_Mwk(_pairs)(ws=)` (étiquettes hors grille par
  images de Wigner-Seitz) ; `local_tmatrix.defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N)` (M_W(R, R) − C_N·𝕀₅ sur les N² mailles de la boîte, étiquettes
  brutes, puis recentrage) ; `config.matrices_dir` / `alignment_C` ; `results_dir = results/M2_plateau`, `matrices_dir = results/M2`, `results_dir_frozen`.
- **Portes** : C.0 (C_N = 0 redonne R9 à 0,0 ; V_loc aligné = R9 au bit ; `Mwr_to_Mwk(ws)` sur la grille MP = sans ws à 4e-16) ; D6 d'`analyze_M.py`
  2,9e-15 / 1,3e-15 / 1,4e-10 ; P-b2, P-c2 de l'audit de l'image minimale (seuil porté à 1e-10 eV par Greg, 16/16 profils OK) ; test d'or 5×5 dense non
  rejoué (D12, recopié de R6, 1,80e-13). `results/M2_plateau/MD5SUMS_2026-09-30.txt` (33 fichiers).
- **Valeur non alignée → finale** (les deux seules grandeurs pour lesquelles ce document garde la valeur non alignée) : Re Σ médian, R_cut 3 :
  −26,65 → **+605,08 meV** (R_cut 0…4 : 17,67 / 682,54 / 479,23 / −26,65 / −410,29 → 29,47 / 761,79 / 737,20 / 605,08 / 696,74 ;
  `m_rcut_resigma.csv` des deux répertoires ; colonnes Δ et C_9 × Δn_L dans `campagnes/M/ch4/table_v1_final.md`, section j). Grandeur de Kaasbjerg (Fig. 3,
  R10 C.2, `campagnes/R/R10_plateau/c/c2_results.json`), eV Å², moyenne sur le disque de 0,1471 Å⁻¹, images de Wigner-Seitz : **K′ (insensible à
  l'alignement)** valence 81,177 → 81,168, conduction 82,305 → 82,297 ; **K** valence 76,956 → **81,605**, conduction 78,168 → **82,821** ; D.2 (bloc
  2×2 π/π* à (K, K)) : ½ Tr 75,87 → **86,59** (partie locale 59,53 → 70,25 ; non locale 16,34, inchangée). Médiane des états 3 132,60 → 3 189,01 meV
  (table complète : `campagnes/R/R10_plateau/c/table_v2_plateau.md`).
- **Fig. 13–14 de Kaasbjerg (R8, `campagnes/R/R8_kaasbjerg/out/{dos,spec,7a}/*.json`)** : position du maximum de ρ − ρ₀ (600², η_G 50 meV) à
  c_i = 0,1 % / 1 % : non aligné −0,1625 / −0,1650 eV, final −0,1475 / −0,1500 eV ; gap à K (A_k, 1 %, η_G 25 meV) : non aligné 0,1176 eV, final
  0,1182 eV ; valeurs lues sur l'article : maximum à −0,20 / −0,21 eV, points blancs à K +0,010 / +0,110 (gap 0,100 eV). (R8 a une troisième variante,
  η_t = η_G, signalée dans `campagnes/M/ch4/table_v1_final.md`, non utilisée.)
- **Constats ouverts (non jugés)** :
  1. C.6 : ⟨ΔV⟩_3D des 18×18, 24×24, 27×27 anormal (−8,48 / −62,74 / −99,91 meV contre +1…+29 meV pour 5…15 et 21) ; l'anomalie est dans le vide
     (|z| > 5 Å : −30,4 / −167,5 / −252,9 meV en moyenne), la part du feuillet (|z| ≤ 3 Å) est régulière (+2,905 / +1,605 / +1,261 meV) ; SCF tous
     convergés ; cause inconnue (`campagnes/R/R10_plateau/c/c6_results.json`).
  2. nscf dense 27×27 de production (`defect_uc_dense_27`, source de M_dense_9x9 et de la wannierisation) : 7 avertissements « c_bands: 1 eigenvalues
     not converged » (R5 phase 0, `campagnes/R/R5_base_vs_M/R5_rapport.md` l.51 et l.202).
  3. Le modèle replié « aligné » de R10 A.3 utilise C_9 à site unique (−24,65 meV), pas C_9 de Kumagai–Oba (−25,14 meV) ; aucun modèle replié n'a été
     recalculé avec l'alignement final.
  4. Les points 1–2 et 5 de « Constats ouverts » de R6 (§7) restent ouverts ; le point 1 (convention intensive six tailles, 7,7e-2) est remplacé par
     deux lignes par famille (3m 2,3e-2, non-3m 2,1e-2, OK ; R10 C.3).
