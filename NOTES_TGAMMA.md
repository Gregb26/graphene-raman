# Évaluation de t / Γ (matrice t locale de Wannier) — paramètres réels et contrôles implémentés

Transcription pour le mémoire (chapitre 4, §4.1.5 et section spectrale), faite le 2026-09-18 à partir des fichiers réels :
provenances `fichier:ligne` re-grep ce jour ; valeurs de production lues dans `results/M/*.npz`, `results/M/*.csv` et
`results/M/logs/*.out` (jamais recalculées ici, sauf le test `test_local_green_batch.py`, 10 s, et le comptage d'états
de la grille interne, indiqué comme tel). Chemins relatifs à `ab-initio-defects/`. Même rôle que `NOTES_EPW.md` §1f pour
§5.4 : une ligne par paramètre, une ligne par contrôle, pour confronter la liste des dix contrôles du mémoire à ce que le
code fait réellement. **Les jugements sont laissés au lecteur ; ce document constate.**

## 0. Ce que le script fait, dans l'ordre (chaîne de production, taille de référence 9×9)

Scripts : `scripts/compute_spectral_wannier.py` (carte niveau 1 : R_cut × grille × η), `scripts/resonance_metrics.py`
(Γ(ε), Born, δρ, T̄ en K), `scripts/resonance_criteria.py` (det, valeurs propres, règle de somme, Γ à c = 0,1 %),
`scripts/rcut_resigma.py` (Re/Im Σ par R_cut). Tous chargent `config/production.json` (`config.py:10–19`, clés
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

## 1. Paramètres réels (taille de référence 9×9)

| paramètre | valeur | fichier:ligne / preuve |
|---|---|---|
| Super-cellule de référence, M | 9×9 (N_cells = 81), lacune sur le sous-réseau A, s_red = (13/27, 13/27) ; M dense zero-padded p = 3, D = 27, forme (20, 729, 20, 729), Hartree, norme `unit_cell` | `config/production.json:7` ; `results/M/M_dense_9x9.json` (sidecar) |
| Cellule primitive dense (nscf) | grille 27×27×1 = 729 k, `nosym` + `noinv` (grille complète), nbnd = 20, `diago_full_acc`, ecutwfc 50 (unités Ha du XML, soit 100 Ry), ecutrho 200, smearing m-v degauss 5,0e-3 (Ha du XML), conv_thr 5e-15, `assume_isolated = 2D` | `scratch/qe_tmp/defect_uc_dense_27/defect_uc_dense_27.save/data-file-schema.xml` (`<nk>`, `<nbnd>`, `<nosym>`, `<noinv>`, `<ecutwfc>`, `<smearing>`, `<assume_isolated>`) ; `config/production.json:16` (nbnd_dense 20) |
| Wannierisation 27×27 | num_wann = 5 (C1 : sp², p_z ; C2 : p_z), 20 bandes, fenêtre externe [−25, 15] eV, fenêtre gelée [−25, −1,74] eV (= E_D + 2,5 eV), Ω_I = 3,0246, Ω_D = 0,0136, Ω_OD = 0,7285, Ω_tot = 3,7668 Å² ; étalements 0,611 (×3, sp²) et 0,967 Å² (×2, p_z) ; run_id `baa17b88b1b69e51` | `wannier/27x27/wannier.wout` l. 125 (grille), 863 (num_wann), 1012–1013 (fenêtres), 2553–2555 (Ω), « Final State » ; `config/production.json` clé `wannier` ; `wannier/27x27/wannier_manifest.json` (sha256 de tb / u / u_dis) |
| Maille | a₁ = (2,135490, −1,232926, 0), a₂ = (2,135490, 1,232926, 0) Å → a = 2,4659 Å, |b| = 4π/(√3 a) = 2,9423 Å⁻¹ ; c = 15,875 Å | `wannier/27x27/wannier.wout` l. 89–92 |
| Unités | M stocké en Ha, multiplié par 27,211386245988 **une fois** au chargement ; H(R) Wannier et toutes les énergies de la matrice t en eV ; Γ rapporté en eV (×10³ = meV) | `config.py:5`, `matrix_io.py:24, 77–78` ; `config/production.json` clé `units` |
| R_cut | 3, sur la **norme euclidienne des indices réduits** (i² + j² ≤ 9, `np.linalg.norm(R_mwr)`), soit n_L = 29 mailles et dim V_loc = 29 × 5 = 145 | `config/production.json:3` ; `compute_spectral_wannier.py:101` ; `resonance_metrics.py:36` ; log `resonance_9x9_20421370.out` « [cluster] R_cut=3: nL=29 sites, dim=145 » |
| Grille de sortie N_k^out | 240×240 = 57 600 k (Γ-centrée, non décalée, K = (2/3, 1/3) sur la grille : indice 38480) × 5 bandes = 288 000 états ; **41 266 états** dans la fenêtre ±3 eV | `config/production.json:4` ; `local_tmatrix.py:24–28` ; `resonance_metrics.py:49, 51` ; log « 41266 on-shell states in +-3.0 eV on 240x240 » ; « k_out[38480]=[0.667 0.333 0.] » |
| **Grille interne N_k^int (pour g₀)** | **300×300 = 90 000 k**, Δk = |b|/300 = 0,00981 Å⁻¹ ; découplée de la grille de sortie ; valeur figée ET valeur par défaut codée en dur (`k_int=None` → 300) ; **jamais variée** : les 53 occurrences dans `results/M/logs/*.out`, `scripts/*.py`, `scripts/*.sh` valent 300 (grep 2026-09-18) | `config/production.json:13` ; `compute_spectral_wannier.py:49, 89` ; `resonance_metrics.py:22, 46` ; `resonance_criteria.py:23, 32` ; `rcut_resigma.py:22, 36` ; défaut : `local_tmatrix.py:171–172, 240–241` |
| Échantillonnage de la grille interne près de E_D (compté ce jour sur les bandes Wannier 27×27) | états à ‖ε − E_D‖ ≤ η = 0,02 eV : **4** (les deux points K × π, π*) ; ≤ 0,1 eV : 52 ; ≤ 0,5 eV : 1 300 ; ≤ 3 eV : 64 528. Pour comparaison : 240² → 4 / 28 / 832 / 41 266 ; 600² → 4 / 220 / 5 140 / 258 196 | calcul direct (`lt.Hwr_to_Hwk` sur `mp_grid(n, n, 1)`, `wannier/27x27/wannier_tb.dat`), non versionné |
| η | 0,02 eV (lorentzienne ε + iη dans g₀ ; même η pour la moyenne lorentzienne des courbes et pour ρ₀) | `config/production.json:5` ; `local_tmatrix.py:150` ; `resonance_metrics.py:75` |
| Fenêtre d'énergie | ±3,0 eV autour de E_D (états sur couche retenus) ; grille d'énergie de g₀ : [E_D − 3 − η, E_D + 3 + η] | `config/production.json:14` ; `resonance_metrics.py:51, 55` ; `local_tmatrix.py:258` |
| Pas de la grille d'énergie | η/ne_per_eta = 0,02/8 = **2,5 meV**, 2 417 énergies (g₀ et dg₀/dε) ; courbes énergétiques sur `egrid[::2]` = 1 209 points (5 meV) | `config/production.json:15` ; `resonance_metrics.py:55, 74` ; log « [g0] 2417 energies, spacing 2.50 meV » |
| E_D (point de Dirac, Wannier) | **−4,238896 eV** (gap minimal 0,0 meV sur 90×90 : K est sur la grille) ; égal à `efermi_read` d'EPW (chapitre 5, −4,2389 eV) | `resonance_metrics.py:47–48` ; log « [dirac] E_D = -4.2389 eV » ; npz `E_D` |
| Grille de ρ₀ | 900×900 (`--rho0-grid`), ρ₀ par cellule et par spin ; `rho0_240` (sur la grille de sortie) aussi stocké | `resonance_metrics.py:18, 80–82` |
| Concentration | c = 1 % par cellule pour ρ_dis = ρ₀ + c δρ ; c = 0,1 % pour la comparaison Kaasbjerg | `config/production.json` clé `defect_concentration_for_dos` ; `resonance_metrics.py:83` ; `resonance_criteria.py:19` |
| Born | t_B = V + V g₀ V (deux premiers termes) → Γ_Born = −2 Im ⟨V g₀ V⟩ (le terme du 1er ordre est réel) | `resonance_metrics.py:60` |
| Alignement du potentiel | décalage G ≈ 0 = moyenne diagonale de M^L = **67,0 meV** (norme unit_cell) ; Γ_T recalculé avec M − 67 meV·1 | `resonance_metrics.py:25–26, 42, 61` |
| δρ | δρ(ε) = (1/π) Im Tr[t(ε) dg₀/dε] par défaut et par spin ; contrôle Lloyd δρ = −(1/π) d/dε Im ln det[1 − g₀V] | `resonance_metrics.py:62` ; `resonance_criteria.py:57–66` |
| T̄ en K | paire π (bandes 3, 4 de 5), trace/2 de ⟨n K| t |n′ K⟩ | `resonance_metrics.py:86–90` |
| Test d'or | η = 0,10 eV, k_int = k_out = grille grossière, R_local = boîte MP-duale complète (exact) ; seuil rel < 1e-8 et Γ ≥ −1e-8 | `scripts/test_local_tmatrix_real.py:20, 53, 59, 61` |
| Coût (16 cœurs) | carte niveau 1 9×9 (36 combinaisons) 1 h 15 ; `resonance_metrics` 7 min 39 ; `resonance_criteria` 1 h 34 ; test ne_per_eta 3 min 35 | `sacct` jobs 20294202, 20421370, 20414406, 20684497 |

## 2. Valeurs de production (9×9, R_cut 3, 240², η 0,02, N_k^int 300²)

| quantité | valeur | source |
|---|---|---|
| médiane de \|Γ\|·N_cells sur les 41 266 états (carte niveau 1) | **2 473,55 meV** (Γ par défaut, ×81 = intensif) | `results/M/level1_summary.csv` l. « 9x9,3,240,0.02 » ; `specwd_9x9_prod.npz` (job 20294202, 2026-09-05) |
| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 2 343,59 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 7 405,83 meV ; Born/T médian 3,302 (min 0,677, max 16,45) | log `resonance_9x9_20421370.out` |
| Re Σ médian, R_cut 3 | 722,3 meV ; \|Re Σ\|/Γ médian 0,227 | `results/M/m_rcut_resigma.csv` |
| position de résonance E_res − E_D (argmax des états à ±1,5 eV) | −1,238 eV (R_cut 3, 240², 0,02) ; −1,183 (η 0,01) ; −1,292 (R_cut 2 et 4) : valeurs discrètes | `level1_summary.csv` ; `m_rcut_resigma.csv` |
| pics des courbes (ε − E_D) | Γ_T : −1,24 eV ; Γ_Born : +1,695 eV ; Γ_T/ρ₀ : **−0,015 eV** (≈ E_D, zone à 4 états internes dans ±η) ; δρ et ρ_dis : −2,53 eV ; \|T̄\| : −0,905 ; −Im T̄ : −1,29 ; min \|Re T̄\| : −2,145 ; Re T̄(E_D) = 2,521 eV, Im T̄(E_D) = −0,091 eV ; aucun zéro de Re T̄ | `resonance_9x9.npz` (`peak_*`, `Tbar_zero_crossings = []`) |
| critère det / valeur propre | min de \|det[1 − Vg₀]\|/max = 2,09e-4 à −2,530 eV ; min_i \|λ_i\| = 0,0108 au même point (λ = +0,0003 + 0,0108 i) ; minima secondaires à −2,19…−1,97 eV (\|det\|/max 2,6e-2…4,1e-2) | log `resonance_criteria_9x9_20414406.out` |
| règle de somme de Friedel | ∫δρ sur [−25,28, 8,66] eV (6 790 points, pas 5 meV) = **−0,0569** états (Tr[t g₀′]) et −0,0569 (Lloyd) ; dans ±3 eV : +1,782 ; écart ponctuel max entre les deux formules 0,614 états/eV ; cumul aux bords : −2,190 (−3 eV), −0,409 (+3 eV) | idem |
| Γ_T à c = 0,1 % sur ±1 eV | min 0,63 meV (+0,24 eV), max 5,55 meV (−1,00 eV), 1,27 meV à E_D ; ħ/Γ à ∓0,3 eV : 419 / 1 025 fs | idem |
| localité de M_W (9×9 dense) | ‖M_W(0, 0)‖ = 9,365 eV ; p_z–p_z sur le site de la lacune 6,617 eV, p_z de l'autre sous-réseau 0,044 eV ; 1ʳᵉ couronne 0,658 / 0,437 eV ; 2ᵉ 0,460 / 0,165 / 0,035 eV ; 3ᵉ ≈ 0,08 eV | `results/M/mwr_locality.npz` (`9x9_dense_*`) ; log `specwd_20294202.out` |
| recentrage | R_d = [4, 4, 0] sur la boîte 27×27 (cohérent avec s_red = 13/27 = 4·3 + 1, p = 3) | log `specwd_20294202.out` « [recenter] » |

## 3. Contrôles réellement implémentés (à confronter aux dix du mémoire)

Type : **P** = porte bloquante (le code refuse, `raise`), **T** = test PASS/FAIL, **C** = convergence (paramètre varié),
**K** = cohérence physique (calculée et rapportée, sans seuil).

| n° | contrôle | type | implémentation | seuil | verdict (date, preuve) |
|---|---|---|---|---|---|
| C1 | Provenance de jauge : tb / u / u_dis d'un même run Wannier90, sha256 recalculés | P | `wannier_provenance.py:42–68`, appelé `compute_spectral_wannier.py:52`, `resonance_metrics.py:21`, `resonance_criteria.py:20` | tout écart ⇒ `ValueError` | « [gauge] provenance OK » dans tous les logs de production (manifeste 27×27) |
| C2 | Normalisation de Bloch et unités de M : sidecar `bloch_norm = unit_cell` (potentiel intensif) et `units = hartree`, conversion eV unique | P | `matrix_io.py:50–78` ; contrat commenté `compute_spectral_wannier.py:60–66` | sidecar absent / mauvais tag ⇒ `ValueError` | en vigueur depuis d8b26fe (2026-09-07) ; M_dense_9x9.json conforme |
| C3 | Cohérence des grilles k : U/U_dis réordonnés sur la grille de M, point sans correspondant ⇒ erreur ; grille MP complète non décalée | P | `wannier_interpolation.py:142–183` (`_match_kpoint_order`, `_infer_mp_grid`) | tol 1e-5 ; N₁N₂N₃ ≠ nk ⇒ `ValueError` | implicite dans chaque run |
| C4 | Garde-fou a : défaut centré (‖M_W(R,0)‖ maximal en R = 0) | P | `local_tmatrix.py:36–51` ; `compute_spectral_wannier.py:83` ; `resonance_metrics.py:35` | argmax ≠ R₀ ⇒ `AssertionError` | passé (valeurs §2) |
| C5 | Garde-fou b : hermiticité de V_loc | P (symétrise + journalise) | `local_tmatrix.py:54–76` | résidu > 1e-10 ⇒ message | résidu 2,3e-15 – 4,0e-15 (logs golden 2026-09-05/07) |
| C6 | **Test d'or** : matrice t locale = matrice T dense (`compute_T`) sur la même grille grossière, même sous-espace 5 WF | T, bloquant avant production | `scripts/test_local_tmatrix_real.py` ; `scripts/test_local_tmatrix.py` (synthétique) ; `scripts/test_local_rcut.py` (support tronqué ⇒ écart) | rel < 1e-8 et Γ ≥ −1e-8 | PASS : 5×5 dense 2,3e-14 (2026-09-05, après le correctif d'unités), 6×6 6,06e-9, 12×12 3,68e-9 (2026-09-07) ; logs `golden_dense_20294198/20444392/20435998.out` |
| C7 | Positivité Γ_nk ≥ 0 | P | `local_tmatrix.py:269–271` ; `rcut_resigma.py:55` | min Γ < −1e-8 ⇒ `AssertionError` | min Γ = +0,26 eV (golden) ; jamais déclenché |
| C8 | g₀ batché = g₀ de référence (restructuration exacte) et évaluation « point d'énergie le plus proche » vs état par état | T | `scripts/test_local_green_batch.py:15` | rel < 1e-12 | ≤ 1,21e-14 pour R_cut 0–3 (rejoué 2026-09-18, grille 90², 7 énergies) ; nearest-grid vs exact : 5,4e-3 avec ne_per_eta 8 (cas synthétique) |
| C9 | Pas de la grille d'énergie ne_per_eta | C | `rcut_resigma.py --npe` (d2d0bac, 2026-09-09) | — | ne_per_eta 2 → 32 : médiane Γ 2 476,19 → 2 474,33 meV (0,08 %), E_res inchangé ; log `npe_test_20684497.out` |
| C10 | R_cut (support de V_loc) | C | carte niveau 1 R_cut 0–3 ; `rcut_resigma.py` R_cut 0–4 | plateau ≤ 5 % (énoncé `compute_spectral_wannier.py:118`) | 9×9 : 2 596,9 / 2 488,9 / 2 466,8 / 2 473,5 / 2 468,0 meV (R_cut 0…4) ; écart médian par état à R_cut 4 : 10,0 / 4,6 / 1,5 / **0,63 %** ; Re Σ 515 / 694 / 723 / 722 / 720 meV |
| C11 | Grille de sortie × η (plateau conjoint) | C | carte niveau 1 (grilles 60/120/240, η 0,05/0,02/0,01) | ≤ 5 % quand η/2 et grille ×2 | 9×9, R_cut 3 : 120² → 240² : 0,07 % (η 0,02), 0,3 % (η 0,01) ; η 0,02 → 0,01 à 240² : 0,001 % ; 0,05 → 0,02 : 0,74 % ; `level1_summary.csv` |
| C12 | Taille de super-cellule N (niveau 2, familles N mod 3) | C | `scripts/level2_families.py`, `results/M/level2_summary.csv`, `level2_families.csv` | N ≥ 7 (config l. 6) | 5/6/7/8/9/12 : 2 524,3 / 2 509,2 / 2 487,5 / 2 478,5 / 2 473,5 / 2 458,7 meV ; 7–9 : étendue 0,56 % ; famille 3m (6, 9, 12) : 2,0 % |
| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | Born/T médian 3,30 sur ±3 eV |
| C14 | Indépendance à l'alignement du potentiel (composante G ≈ 0) | K | `resonance_metrics.py:25–26, 42, 61` | — | décalage 67,0 meV ⇒ écart relatif max 6,4e-4 (médian 1,8e-4) sur Γ_T |
| C15 | Règle de somme de Friedel, deux formules (Tr[t g₀′] vs Lloyd) | K | `resonance_criteria.py:57–75` | — | −0,0569 / −0,0569 états sur toute la bande ; +1,78 dans ±3 eV |
| C16 | Critère de résonance (det, valeur propre minimale) et position du pic | K | `resonance_criteria.py:41–55` ; `compute_spectral_wannier.py:112` | — | minimum unique à −2,530 eV ; E_res(argmax Γ) = −1,24 eV |
| C17 | K sur la grille de sortie et dégénérescence π/π* en K | K | `resonance_metrics.py:86` ; `config/production.json` clé `K_red` | — | k_out[38480] = (2/3, 1/3), E(π) = E(π*) = E_D |
| — | **Grille interne N_k^int** | **aucun** | valeur 300 partout, jamais balayée | — | voir §4 |

## 4. Points ouverts constatés le 2026-09-18

1. **N_k^int n'a jamais été varié.** Toute la production (cartes niveau 1, niveau 2, résonance, critères, Re Σ, test
ne_per_eta) utilise 300² = 90 000 points ; c'est aussi le défaut codé en dur de `scattering_rate` /
`scattering_rate_fast`. Aucun log ni npz ne porte une autre valeur. Le plateau grille × η (C11) porte sur la grille de
**sortie** et ne teste pas la somme interne de g₀. Près de E_D, la grille interne ne contient que 4 états dans ±η (les
deux K), 52 dans ±0,1 eV ; le pic de Γ_T/ρ₀ est justement à −0,015 eV. Pour chiffrer : `scripts/rcut_resigma.py`
accepte désormais `--nk-int` (même mécanisme que `--npe`, valeur inscrite dans le npz) et
`scripts/submit_nkint_check.sh 9x9 "150 300 450 600"` enchaîne les quatre valeurs à R_cut 3 (≈ 4 × (1 → 4) × 1 min de g₀,
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
(`scripts/test_local_tmatrix_real.py` ; `test_local_rcut.py` corrigé de même, « must be ~1e-8 » → « seuil 1e-8 » ; `test_local_tmatrix.py` n'avait pas de message trompeur.)

- `scripts/rcut_resigma.py` : option `--nk-int` (2 lignes, même patron que `--npe` ; défaut = config figée, valeur
  inscrite dans le npz et le log).
- `scripts/submit_nkint_check.sh` : job de balayage N_k^int (lancé : job 21337627, §6).
- `scripts/nkint_check_post.py` : post-traitement du balayage (tableaux de §6, `results/M/nkint_check_9x9.csv`).
- `NOTES_TGAMMA.md` : ce document.

## 6. Balayage de la grille interne N_k^int (P13, 2026-09-18)

Job Slurm 21337627 (`scripts/submit_nkint_check.sh 9x9 "150 300 450 600"`, rc32615, 16 cœurs, 64 G demandés, MaxRSS 14906672K,
départ 08:27:19, fin 09:04:59 EDT, 37 min 40 ; 8 à 10 min par valeur, dominées par le chargement de M et la rotation, pas par g₀).
Chaîne identique à la production : `scripts/rcut_resigma.py --nk-int N` (option ajoutée ce jour, l. 20 et 22 ; `nk_int` inscrit
dans le npz l. 41 et dans le log l. 34), R_cut 3, grille de sortie 240², η 0,02, fenêtre ±3 eV, ne_per_eta 8, portes C1–C7
actives (« [gauge] provenance OK », recentrage R_d = [4, 4, 0], positivité l. 55). `config/production.json` inchangé (nk_int 300).
Dry run préalable (nœud de connexion, nk_int 60, grille 30, 17 s) : valeur présente dans le log et le npz.
Sources : `results/M/resigma_9x9_rc3_nk{150,300,450,600}.npz` (clés `Sigma_rc3`, `E_out`, `E_D`, `nk_int` ; non commis),
log `results/M/logs/nkint_21337627.out`, post-traitement `scripts/nkint_check_post.py` → `results/M/nkint_check_9x9.csv` (commis).
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
