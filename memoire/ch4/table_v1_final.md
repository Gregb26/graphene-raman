# Chiffres du chapitre 4 — table v1 → final

Date : 2026-09-30 ; HEAD `90d2022` ; md5 `article/R6_production_corrigee/etape3/table_v1_v2.md` f8972c4de041dbd79e70b2899aa4f941 ; `article/R10_plateau/c/table_v2_plateau.md` 807bb352fcb989f01d94481955fe4c10 ; `results/M2_plateau/MD5SUMS_2026-09-30.txt` 9c29daedb9c417ba795bb39c0ed71bee. Généré par `memoire/ch4/ch4_chiffres.py table` (ne pas éditer à la main).

## Nomenclature

- **non aligné** (= « tel quel » des archives) : ΔV brut, v2 (M2 = N_cells·M^L + M^NL), `results/M2` ;
- **alignement à site unique** (= « Lu » des archives) : potentiel au site le plus éloigné de la lacune (Lu et al. 2019) ;
- **alignement de Kumagai–Oba** (= « plateau (i) » des archives) : moyenne des potentiels de site (sphères de 1 Å) sur les atomes à ≥ 0,75 r_max de la lacune (Kumagai et Oba, PRB 89, 195205, 2014 ; région élargie, voir `alignement_regions.md`) ;
- symbole du mémoire ΔV_PA^(N) = C_N = `alignment.C_N_eV` de `config/production.json` ;
- **final** = v2 + alignement de Kumagai–Oba, `results/M2_plateau` (R10). Les archives (rapports R4–R10, tables sources) gardent leurs étiquettes.

## Règles

Appariement : « grandeur » normalisée (NFKC, minuscules, × → x, espaces réduits) après les réécritures d'étiquettes suivantes (R6 → R10) : «  : brut ; aligné (D6) » → «  » ; «  : brut ; aligné » → «  » ; « médiane γ (mev) / re σ (mev) » → « médiane γ / re σ (mev) » ; « médiane de la courbe γ_t sur ±3 ev » → « médiane courbe γ_t ±3 ev » ; «  (resonance_metrics) » → «  » ; « (mev, resonance_metrics) » → « (mev) » ; « (à ε − e_d, ev) » → « (à ε − e_d) » ; « (ratio = 1, ev) » → « (ratio = 1) » ; « c14 redéfini (v_loc + c·1) : » → « c14 m2 (±25 mev, v_loc tel quel) : » ; « (matrice complète) » → « (complet) » ; les 9 lignes « critère/Friedel 9x9 » de R6 (det, λ_min, Friedel × complet/π/σ) s'apparient aux 3 lignes combinées de R10 (cellule découpée sur « ; ») ; les lignes v1 seulement 48, 51, 77, 80 prennent leur valeur finale dans un npz de `results/M2_plateau` (D1). Statut mécanique : tous les nombres d'une cellule sont extraits dans l'ordre ; « inchangé » si même compte et chaque |final − v1| ≤ 1e-6 relatif ; « remplacé » sinon ; « nouveau » sans v1 ; « sans équivalent final » sans final ; « retiré (décision) » pour la liste du prompt seulement (E_res, C14 67/±25 meV, zéro de det −2,530 eV, convention intensive six tailles), avec la référence, le statut mécanique étant conservé après « ; ». « final − v1 » : élément par élément, au nombre de décimales le plus grand des deux cellules ; « — » si les comptes diffèrent. Valeurs copiées telles quelles (format des tables sources ; npz : repr).

Comptes : R6 126 lignes, R10 192 lignes ; appariées 117 ; v1 seulement 9 ; final seulement 81 ; total 207. Statuts : inchangé 6, nouveau 87, remplacé 85, retiré (décision) 24, sans équivalent final 5.

## Table principale (ordre : lignes de R6, puis lignes propres à R10)

| grandeur | v1 (mémoire actuel) | non aligné (v2) | final | final − v1 | statut | source finale (fichier : clé ou ligne) | lieu |
|---|---|---|---|---|---|---|---|
| médiane \|Γ\|·N_cells, 5x5, R_cut 3, 240², η 0.02 (meV) | 2524.30 | 3288.90 | 2997.14 | +472.84 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.18 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 5x5 (eV) | -1.127 | -0.284 | -0.227 | +0.900 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.19 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 6x6, R_cut 3, 240², η 0.02 (meV) | 2509.20 | 3155.98 | 3149.15 | +639.95 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.20 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 6x6 (eV) | -1.127 | -0.227 | -0.175 | +0.952 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.21 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 7x7, R_cut 3, 240², η 0.02 (meV) | 2487.46 | 3052.03 | 2942.85 | +455.39 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.22 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 7x7 (eV) | -1.127 | -0.269 | -0.269 | +0.858 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.23 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 8x8, R_cut 3, 240², η 0.02 (meV) | 2478.47 | 3020.42 | 2974.18 | +495.71 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.24 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 8x8 (eV) | -1.238 | -0.270 | -0.227 | +1.011 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.25 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 3, 240², η 0.02 (meV) | 2473.55 | 3132.60 | 3189.01 | +715.46 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.26 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 9x9 (eV) | -1.238 | -0.175 | -0.175 | +1.063 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.27 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 12x12, R_cut 3, 240², η 0.02 (meV) | 2458.74 | 3162.78 | 3264.51 | +805.77 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.28 | fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA |
| E_res − E_D (argmax états ±1,5 eV), 12x12 (eV) | -1.292 | -0.175 | -0.175 | +1.117 | retiré (décision) — E_res (R10 B.1) ; mécanique : remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.29 | fig_plateau, §2 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 3, 240², η 0.01 (meV) | 2473.51 | 3129.38 | 3202.75 | +729.24 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.30 | fig_plateau, fig_rcut, C10/C11 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 3, 120², η 0.02 (meV) | 2475.21 | 3148.05 | 3200.71 | +725.50 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.31 | fig_plateau, fig_rcut, C10/C11 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 2, 240², η 0.02 (meV) | 2466.84 | 3170.21 | 3222.08 | +755.24 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.32 | fig_plateau, fig_rcut, C10/C11 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 4, 240², η 0.02 (meV) | — | 3147.47 | 3200.00 | — | nouveau | `level1_summary.csv` ; `table_v2_plateau.md` l.33 | fig_plateau, fig_rcut, C10/C11 |
| médiane \|Γ\|·N_cells, 9x9, R_cut 0, 240², η 0.02 (meV) | 2596.92 | 3009.73 | 3013.96 | +417.04 | remplacé | `level1_summary.csv` ; `table_v2_plateau.md` l.34 | fig_plateau, fig_rcut, C10/C11 |
| p_z–p_z sur le site de la lacune, 6x6 dense (eV) | 7.007 | 31.502 | 31.552 | +24.545 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.35 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 6x6 (eV) | +0.3154 / +3.1030 | +11.3561 / +3.1030 | +11.3561 / +3.1030 | +11.0407 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.36 | fig_M_map, §4.1.5 |
| p_z–p_z sur le site de la lacune, 9x9 dense (eV) | 6.617 | 31.521 | 31.546 | +24.929 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.38 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 9x9 (eV) | +0.1396 / +3.1030 | +11.3051 / +3.1030 | +11.3051 / +3.1030 | +11.1655 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.39 | fig_M_map, §4.1.5 |
| p_z–p_z sur le site de la lacune, 12x12 dense (eV) | 6.483 | 31.528 | 31.547 | +25.064 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.41 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 12x12 (eV) | +0.0779 / +3.1030 | +11.2220 / +3.1030 | +11.2220 / +3.1030 | +11.1441 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.42 | fig_M_map, §4.1.5 |
| p_z–p_z sur le site de la lacune, 5x5 dense (eV) | 7.296 | 31.002 | 31.060 | +23.764 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.44 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 5x5 (eV) | +0.4214 / +3.1004 | +10.5343 / +3.1004 | +10.5343 / +3.1004 | +10.1129 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.45 | fig_M_map, §4.1.5 |
| p_z–p_z sur le site de la lacune, 7x7 dense (eV) | 6.807 | 30.936 | 30.963 | +24.156 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.47 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 7x7 (eV) | +0.2140 / +3.1009 | +10.4857 / +3.1009 | +10.4857 / +3.1009 | +10.2717 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.48 | fig_M_map, §4.1.5 |
| p_z–p_z sur le site de la lacune, 8x8 dense (eV) | 6.691 | 31.092 | 31.107 | +24.416 | remplacé | `level2_families.csv (mwr_locality.npz)` ; `table_v2_plateau.md` l.50 | fig_locality, tab familles |
| Re M^L / Re M^NL à K, paire π, 8x8 (eV) | +0.1677 / +3.1009 | +10.7343 / +3.1009 | +10.7343 / +3.1009 | +10.5666 / +0.0000 | remplacé | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.51 | fig_M_map, §4.1.5 |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 5x5 dense | 6.445130123841091 | 0.2578052049536437 | 0.2578052049536437 ; 0.2572153016716182 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.53 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 5x5 (eV) | 0.9774808691110594 | 24.43702172777648 | 24.43702172777648 ; 24.493066135813088 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.54 | tab:L_NL |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 6x6 dense | 9.087989279935528 | 0.25244414666487575 | 0.25244414666487575 ; 0.2519315893483562 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.55 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 6x6 (eV) | 0.6897427076471349 | 24.830737475296853 | 24.830737475296853 ; 24.881255856895844 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.56 | tab:L_NL |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 7x7 dense | 12.656662128876818 | 0.25829922711993497 | 0.25829922711993497 ; 0.25802354135325317 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.57 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 7x7 (eV) | 0.49705210093257957 | 24.355552945696406 | 24.355552945696406 ; 24.381575684751834 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.58 | tab:L_NL |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 8x8 dense | 16.426469175893885 | 0.25666358087334196 | 0.25666358087334196 ; 0.2565179816585835 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.59 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 8x8 (eV) | 0.38306212443794735 | 24.51597596402863 | 24.51597596402863 ; 24.529891194556797 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.60 | tab:L_NL |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 9x9 dense | 20.44280135341195 | 0.25238026362236965 | 0.25238026362236965 ; 0.2521218195055461 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.61 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 9x9 (eV) | 0.3087352524875175 | 25.00755545148893 | 25.00755545148893 ; 25.03319010538462 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.62 | tab:L_NL |
| ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, 12x12 dense | 36.309206345434944 | 0.2521472662877427 | 0.2521472662877427 ; 0.25195367625385306 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.63 | tab:L_NL |
| ⟨‖M^L‖_F⟩, 12x12 (eV) | 0.17263870417258428 | 24.859973400852134 | 24.859973400852134 ; 24.879074702189435 | — | remplacé | `lnl_frobenius.csv` ; `table_v2_plateau.md` l.64 | tab:L_NL |
| localité M_W 9x9 dense : ‖M_W(0,0)‖ (eV) | 9.365 | 41.879 | 41.927 | +32.562 | remplacé | `mwr_locality.npz` ; `table_v2_plateau.md` l.72 | fig_locality, §2 |
| localité M_W 9x9 dense : p_z–p_z lacune (eV) | 6.617 | 31.521 | 31.546 | +24.929 | remplacé | `mwr_locality.npz` ; `table_v2_plateau.md` l.73 | fig_locality, §2 |
| localité M_W 9x9 : ‖M_W(R,0)‖ première couronne (eV, 4 premières) | 0.658, 0.437, 0.437, 0.658 | 1.963, 1.388, 1.388, 1.963 | 1.963, 1.388, 1.388, 1.963 | +1.305 / +0.951 / +0.951 / +1.305 | remplacé | `results/M2_plateau/mwr_locality.npz` : 9x9_dense_w (D1, hors table R10) | fig_locality |
| localité M_W 12x12 dense : ‖M_W(0,0)‖ (eV) | 6.543 | 35.328 | 35.354 | +28.811 | remplacé | `mwr_locality.npz` ; `table_v2_plateau.md` l.75 | fig_locality, §2 |
| localité M_W 12x12 dense : p_z–p_z lacune (eV) | 6.483 | 31.528 | 31.547 | +25.064 | remplacé | `mwr_locality.npz` ; `table_v2_plateau.md` l.76 | fig_locality, §2 |
| localité M_W 12x12 : ‖M_W(R,0)‖ première couronne (eV, 4 premières) | 0.096, 0.096, 2.749, 2.749 | 0.756, 0.756, 2.328, 2.328 | 0.756, 0.756, 2.328, 2.328 | +0.660 / +0.660 / -0.421 / -0.421 | remplacé | `results/M2_plateau/mwr_locality.npz` : 12x12_dense_w (D1, hors table R10) | fig_locality |
| tab:rcut_M 9x9 R_cut 0 : max\|ΔM\|/max\|M\| (60²) | 4.0851e-01 | 3.1569e-01 | 3.1708e-01 (plateau, étiquettes actuelles : 3.1708e-01) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.78 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 1 : max\|ΔM\|/max\|M\| (60²) | 1.7367e-01 | 2.1293e-01 | 1.6301e-01 (plateau, étiquettes actuelles : 1.6301e-01) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.79 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 2 : max\|ΔM\|/max\|M\| (60²) | 4.9711e-02 | 1.1395e-01 | 5.5681e-02 (plateau, étiquettes actuelles : 5.5681e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.80 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 3 : max\|ΔM\|/max\|M\| (60²) | 2.2704e-02 | 6.6789e-02 | 2.5879e-02 (plateau, étiquettes actuelles : 2.4704e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.81 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 0 : max\|ΔM\|/max\|M\| (60²) | 7.3098e-01 | 8.6729e-01 | 8.6782e-01 (plateau, étiquettes actuelles : 8.6781e-01) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.85 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 1 : max\|ΔM\|/max\|M\| (60²) | 1.8606e-01 | 2.1719e-01 | 1.8925e-01 (plateau, étiquettes actuelles : 1.8923e-01) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.86 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 2 : max\|ΔM\|/max\|M\| (60²) | 1.1630e-01 | 1.2490e-01 | 1.2506e-01 (plateau, étiquettes actuelles : 1.2505e-01) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.87 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 3 : max\|ΔM\|/max\|M\| (60²) | 2.5817e-02 | 9.0206e-02 | 3.0979e-02 (plateau, étiquettes actuelles : 3.0966e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.88 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 4 : max\|ΔM\|/max\|M\| (60²) | 1.5740e-02 | 4.4189e-02 | 1.8194e-02 (plateau, étiquettes actuelles : 1.7179e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.82 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 5 : max\|ΔM\|/max\|M\| (60²) | 1.2464e-02 | 1.7488e-02 | 1.4950e-02 (plateau, étiquettes actuelles : 1.3578e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.83 | tab:rcut_M |
| tab:rcut_M 9x9 R_cut 6 : max\|ΔM\|/max\|M\| (60²) | 1.0763e-02 | 1.2068e-02 | 1.3619e-02 (plateau, étiquettes actuelles : 1.1992e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.84 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 4 : max\|ΔM\|/max\|M\| (60²) | 1.5742e-02 | 7.5696e-02 | 1.8573e-02 (plateau, étiquettes actuelles : 1.8398e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.89 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 5 : max\|ΔM\|/max\|M\| (60²) | 1.2370e-02 | 5.3559e-02 | 1.2739e-02 (plateau, étiquettes actuelles : 1.2858e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.90 | tab:rcut_M |
| tab:rcut_M 12x12 R_cut 6 : max\|ΔM\|/max\|M\| (60²) | 1.1527e-02 | 3.2881e-02 | 1.1781e-02 (plateau, étiquettes actuelles : 1.2725e-02) | — | remplacé | `m_rcut_convergence.csv ; c/c2_results.json` ; `table_v2_plateau.md` l.91 | tab:rcut_M |
| N_k^int 150 : médiane Γ (meV) / Re Σ (meV) / E_res / Γ(E_D) | 2398.96 / 725.37 / -1.116 / 248.81 | 3155.57 / -23.25 / -0.238 / 4801.86 | 3210.02 / 607.78 / -0.238 / 5662.16 | +811.060 / -117.590 / +0.878 / +5413.350 | remplacé | `nkint_check_9x9.csv` ; `table_v2_plateau.md` l.97 | §6 NOTES_TGAMMA |
| N_k^int 300 : médiane Γ (meV) / Re Σ (meV) / E_res / Γ(E_D) | 2473.55 / 722.28 / -1.238 / 181.96 | 3132.60 / -26.65 / -0.175 / 3630.65 | 3189.01 / 605.08 / -0.175 / 4287.11 | +715.460 / -117.200 / +1.063 / +4105.150 | remplacé | `nkint_check_9x9.csv` ; `table_v2_plateau.md` l.98 | §6 NOTES_TGAMMA |
| N_k^int 450 : médiane Γ (meV) / Re Σ (meV) / E_res / Γ(E_D) | 2478.71 / 716.87 / -1.238 / 174.15 | 3135.87 / -26.79 / -0.202 / 3489.52 | 3196.11 / 604.19 / -0.175 / 4120.38 | +717.400 / -112.680 / +1.063 / +3946.230 | remplacé | `nkint_check_9x9.csv` ; `table_v2_plateau.md` l.99 | §6 NOTES_TGAMMA |
| N_k^int 600 : médiane Γ (meV) / Re Σ (meV) / E_res / Γ(E_D) | 2476.95 / 718.79 / -1.245 / 172.71 | 3135.65 / -25.09 / -0.202 / 3463.87 | 3191.93 / 605.25 / -0.175 / 4090.11 | +714.980 / -113.540 / +1.070 / +3917.400 | remplacé | `nkint_check_9x9.csv` ; `table_v2_plateau.md` l.100 | §6 NOTES_TGAMMA |
| résonance 9x9 : médiane de la courbe Γ_T sur ±3 eV (meV) | 2343.59 | 3740.17 | 3771.46 | +1427.87 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.101 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : médiane Γ_Born (meV) | 7405.83 | 1.7720e+05 | 1.7944e+05 | +172034.1700 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.102 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : Born/T médian | 3.302 | 46.489 | 45.846 | +42.544 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.103 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : Born/T min | 0.677 | 1.106 | — | — | sans équivalent final | — (v1 : `table_v1_v2.md` l.73, `resonance_9x9.npz`) | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : Born/T max | 16.453 | 211.720 | — | — | sans équivalent final | — (v1 : `table_v1_v2.md` l.74, `resonance_9x9.npz`) | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : moyenne diagonale de M^L (meV, C14 ancien) | 67.00 | 5427.32 | — | — | retiré (décision) — C14 à 67 meV (R10 C.1 : C14 = ±9,05 meV) | — (v1 : `table_v1_v2.md` l.75, `resonance_9x9.npz`) | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic Γ_T (eV) | -1.240 | -0.180 | -0.180 | +1.060 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.104 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic Γ_Born | 1.695 | 1.695 | 1.695 | +0.000 | inchangé | `results/M2_plateau/resonance_9x9.npz` : peak_GB (D1, hors table R10) | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic Γ_T/ρ₀ | -0.015 | -0.175 | -0.170 | -0.155 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.105 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic δρ | -2.530 | -0.815 | -0.787 | +1.743 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.106 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic ρ_dis | -2.530 | 1.710 | 1.710 | +4.240 | remplacé | `results/M2_plateau/resonance_9x9.npz` : peak_rho_dis (D1, hors table R10) | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic \|T̄\| | -0.905 | -0.172 | -0.170 | +0.735 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.107 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : pic −Im T̄ | -1.292 | -0.177 | -0.172 | +1.120 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.108 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : min \|Re T̄\| | -2.145 | -2.447 | -2.140 | +0.005 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.109 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : Re T̄(E_D) (eV) | 2.521 | 11.122 | 12.624 | +10.103 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.110 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : Im T̄(E_D) (eV) | -0.091 | -1.796 | -2.119 | -2.028 | remplacé | `resonance_9x9.npz` ; `table_v2_plateau.md` l.111 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : E_res états ±1,5 eV (resonance_metrics) | — | -0.175 | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_9x9.npz` ; `table_v2_plateau.md` l.112 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : médiane \|Γ_T\| états (meV, resonance_metrics) | — | 3132.872 [R10 : 3132.87] | 3188.35 | — | nouveau | `resonance_9x9.npz` ; `table_v2_plateau.md` l.113 | fig_spectral, fig_epw_vs_ed (ch. 5), §2 |
| résonance 9x9 : zéros de Re T̄ | aucun | -2.447, -0.222, 1.638 [R10 : -2.447, -0.222, +1.638] | -2.142, -2.130, -2.095, -2.070, -2.047, -2.005, -2.000, -0.210 | — | nouveau | `resonance_9x9.npz` ; `table_v2_plateau.md` l.114 | §2 |
| critère 9x9 (matrice complète) : min \|det\|/max (ε − E_D) | 2.091e-04 (-2.530) | 1.318e-04 (-0.812) | 1.261e-04 (-0.785) | -0.000 / +1.745 | retiré (décision) — zéro de det à −2,530 eV = doublet σ (R4 D3) ; critère par bloc (R10 C.1) | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.115 (cellule combinée, élément 1) | §2, C16 |
| critère 9x9 (matrice complète) : λ_min (ε − E_D) | +0.0003+0.0108i (-2.530) | +0.0000+0.0019i (-0.812) [R10 : \|λ\| 0.0019 (-0.812)] | \|λ\| 0.0019 (-0.787) | — | retiré (décision) — idem (R4 D3 ; R10 C.1) | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.115 (cellule combinée, élément 2) | §2, C16 |
| Friedel 9x9 (matrice complète) : ∫δρ bande / Lloyd / ±3 eV (états) | -0.0569 / -0.0569 / +1.782 | -1.0007 / -1.0007 / +0.669 [R10 : -1.0007] | -1.0005 / -1.0005 / +0.702 | -0.9436 / -0.9436 / -1.0800 | remplacé | `results/M2_plateau/resonance_criteria_9x9.npz` : sumrule, sumrule_lloyd, sumrule_window (D1 ; table R10 l.115 : -1.0005) | §2, C15 |
| critère 9x9 (bloc π) : min \|det\|/max (ε − E_D) | — | 2.330e-02 (-0.172) | 1.706e-02 (-0.170) | — | nouveau | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.116 (cellule combinée, élément 1) | §2, C16 |
| critère 9x9 (bloc π) : λ_min (ε − E_D) | — | +0.2186+0.3310i (-0.170) [R10 : \|λ\| 0.3967 (-0.170)] | \|λ\| 0.3440 (-0.127) | — | nouveau | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.116 (cellule combinée, élément 2) | §2, C16 |
| Friedel 9x9 (bloc π) : ∫δρ bande / Lloyd / ±3 eV (états) | — | -0.9980 / -0.9980 / -0.285 [R10 : -0.9980] | -0.9981 / -0.9981 / -0.256 | — | nouveau | `results/M2_plateau/resonance_criteria_9x9.npz` : sumrule_pi, sumrule_lloyd_pi, sumrule_window_pi (D1 ; table R10 l.116 : -0.9981) | §2, C15 |
| critère 9x9 (bloc σ) : min \|det\|/max (ε − E_D) | — | 3.480e-04 (-0.812) | 4.379e-04 (-0.787) | — | nouveau | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.117 (cellule combinée, élément 1) | §2, C16 |
| critère 9x9 (bloc σ) : λ_min (ε − E_D) | — | +0.0000+0.0019i (-0.812) [R10 : \|λ\| 0.0019 (-0.812)] | \|λ\| 0.0019 (-0.787) | — | nouveau | `resonance_criteria_9x9.npz` ; `table_v2_plateau.md` l.117 (cellule combinée, élément 2) | §2, C16 |
| Friedel 9x9 (bloc σ) : ∫δρ bande / Lloyd / ±3 eV (états) | — | -0.0028 / -0.0028 / +0.955 [R10 : -0.0028] | -0.0024 / -0.0024 / +0.958 | — | nouveau | `results/M2_plateau/resonance_criteria_9x9.npz` : sumrule_sigma, sumrule_lloyd_sigma, sumrule_window_sigma (D1 ; table R10 l.117 : -0.0024) | §2, C15 |
| Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV) | 0.63 (+0.24) / 5.55 (-1.00) / 1.27 | 2.25 (+1.00) / 37.84 (-0.18) / 7.79 | — | — | sans équivalent final | — (v1 : `table_v1_v2.md` l.98, `resonance_criteria_9x9.npz`) | §2 (Kaasbjerg) |
| ħ/Γ à ∓0,3 eV, c = 0,1 % (fs) | 419 / 1025 | 24 / 230 | — | — | sans équivalent final | — (v1 : `table_v1_v2.md` l.99, `resonance_criteria_9x9.npz`) | §2 |
| tab:tests_M : hermiticité M dense 5x5 | 1.38e-14 (OK) | 4.09e-15 (OK) | 4.09e-15 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.168 | tab:tests_M |
| tab:tests_M : hermiticité M dense 6x6 | 2.05e-14 (OK) | 5.79e-15 (OK) | 5.79e-15 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.169 | tab:tests_M |
| tab:tests_M : hermiticité M dense 7x7 | 2.83e-14 (OK) | 8.37e-15 (OK) | 8.37e-15 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.170 | tab:tests_M |
| tab:tests_M : hermiticité M dense 8x8 | 4.08e-14 (OK) | 1.19e-14 (OK) | 1.19e-14 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.171 | tab:tests_M |
| tab:tests_M : hermiticité M dense 9x9 | 4.06e-14 (OK) | 1.15e-14 (OK) | 1.15e-14 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.172 | tab:tests_M |
| tab:tests_M : padding : k coïncidents, bandes 1–8, 9x9 | 5.5e-04 (OK) | 5.0e-04 (OK) | 5.0e-04 (OK) | -0.0 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.173 | tab:tests_M |
| tab:tests_M : padding : k coïncidents, bandes 1–15, 9x9 | 1.2e-03 (OK) | 1.3e-03 (OK) | 1.3e-03 (OK) | +0.0 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.174 | tab:tests_M |
| tab:tests_M : non-régression nbnd 16 → 20, 9x9 | 1.2e-07 (OK) | non rejouable (fichier absent) | non rejouable (fichier absent) | — | sans équivalent final | `M_tests_summary.csv` ; `table_v2_plateau.md` l.175 | tab:tests_M |
| tab:tests_M : padding : noyau dense à p=1 vs noyau N×N, 9x9 | 1.3e-15 (OK) | 1.3e-15 (OK) | 1.3e-15 (OK) | +0.0 | inchangé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.176 | tab:tests_M |
| tab:tests_M : fermeture Fourier Wannier (k → R → k), 9x9 | 2.4e-14 (OK) | 2.4e-14 (OK) | 2.4e-14 (OK) | +0.0 | inchangé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.177 | tab:tests_M |
| tab:tests_M : fermeture Bloch → Wannier → Bloch (5 bandes), 9x9 | 1.3e-15 (OK) | 1.3e-15 (OK) | 1.3e-15 (OK) | +0.0 | inchangé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.178 | tab:tests_M |
| tab:tests_M : hermiticité M dense 12x12 | 8.13e-14 (OK) | 2.29e-14 (OK) | 2.29e-14 (OK) | -0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.179 | tab:tests_M |
| tab:tests_M : convention intensive (cellule unitaire) | 2.7e-02 (OK) | 7.7e-02 (À VOIR) | — (ligne retirée) | — | retiré (décision) — test « convention intensive » six tailles (R10 C.3 : lignes par famille) | `M_tests_summary.csv` ; `table_v2_plateau.md` l.180 | tab:tests_M |
| tab:tests_M : test d'or (local vs compute_T), 5×5 dense | 2.3e-14 (OK) | 1.80e-13 (OK) | 1.80e-13 (OK) | +0.00 | remplacé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.181 | tab:tests_M |
| tab:tests_M : g0 par lots vs référence, R_cut 0–3 | 1.2e-14 (OK) | 1.2e-14 (OK) | 1.2e-14 (OK) | +0.0 | inchangé | `M_tests_summary.csv` ; `table_v2_plateau.md` l.182 | tab:tests_M |
| tab:tests_M : porte A.2 : ΔV appliqué aux états de Bloch purs contre M2/N_cells, 8 fichiers grossiers (N = 5x5, 6x6, 7x7, 8x8, 9x9, 9x9 (128 b), 10x10, 12x12) | — | L 3.3e-14 ; NL 8.1e-15 (OK) | L 3.3e-14 ; NL 8.1e-15 (OK) | — | nouveau | `M_tests_summary.csv` ; `table_v2_plateau.md` l.183 | tab:tests_M (nouvelle ligne) |
| tab:tests_M : porte A.2 : ΔV appliqué aux états de Bloch purs contre M2/N_cells, 6 fichiers denses (N = 5x5, 6x6, 7x7, 8x8, 9x9, 12x12) | — | L 1.2e-08 ; NL 4.0e-09 (OK) | L 1.2e-08 ; NL 4.0e-09 (OK) | — | nouveau | `M_tests_summary.csv` ; `table_v2_plateau.md` l.184 | tab:tests_M (nouvelle ligne) |
| C14 redéfini (V_loc + C·1) : E_res_states | — | -0.175 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.152 | 3.5 |
| C14 redéfini (V_loc + C·1) : E_res_states_shiftm25 | — | -0.227 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.153 | 3.5 |
| C14 redéfini (V_loc + C·1) : E_res_states_shiftp25 | — | -0.175 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.154 | 3.5 |
| C14 redéfini (V_loc + C·1) : median_GT_states_meV | — | 3132.872 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.155 | 3.5 |
| C14 redéfini (V_loc + C·1) : median_GT_states_shiftm25_meV | — | 3158.539 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.156 | 3.5 |
| C14 redéfini (V_loc + C·1) : median_GT_states_shiftp25_meV | — | 3187.964 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.157 | 3.5 |
| C14 redéfini (V_loc + C·1) : peak_GT_shiftm25 | — | -0.235 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.158 | 3.5 |
| C14 redéfini (V_loc + C·1) : peak_GT_shiftp25 | — | -0.180 | — | — | retiré (décision) — C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV) | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.159 | 3.5 |
| Γ^ed/Γ^ep (c = 1 %, 300 K, ±3 eV) : médiane | 0.361 | 0.537 | 0.554 | +0.193 | remplacé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.191 | fig_epw_vs_ed, NOTES_EPW (ch. 5) |
| Γ^ed/Γ^ep : min (à ε − E_D, eV) | 0.088 (+1.875) | 0.123 (+1.875) | 0.128 (+1.875) | +0.040 / +0.000 | remplacé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.192 | NOTES_EPW (ch. 5) |
| Γ^ed/Γ^ep : max (à ε − E_D, eV) | 2.011 (-0.755) | 32.442 (-0.180) | 34.797 (-0.180) | +32.786 / +0.575 | remplacé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.193 | NOTES_EPW (ch. 5) |
| Γ^ed/Γ^ep : croisements (ratio = 1, eV) | -1.619, -0.215 | -1.504, +0.690 | -1.486, +0.737 | +0.133 / +0.952 | remplacé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.194 | NOTES_EPW (ch. 5) |
| médiane Γ^ed, c = 1 % (meV) : ±3 eV / ±1,2 eV | 23.497 / 10.737 | 37.462 / 59.298 | 37.730 / 57.508 | +14.233 / +46.771 | remplacé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.195 | fig_epw_vs_ed (ch. 5) |
| médiane Γ^ep, 300 K (meV) : ±3 eV / ±1,2 eV (indépendant de M) | 56.956 / 21.653 | 56.956 / 21.653 | 56.956 / 21.653 | +0.000 / +0.000 | inchangé | `ed_vs_ep_24k24q_mv0.02.npz` ; `table_v2_plateau.md` l.196 | fig_epw_vs_ed (ch. 5) |
| C_N plateau (i), 5x5 (meV) | — | 0 (tel quel) | -57.5626 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.5 | D14, R10 A.1 |
| C_N plateau (i), 6x6 (meV) | — | 0 (tel quel) | -50.5091 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.6 | D14, R10 A.1 |
| C_N plateau (i), 7x7 (meV) | — | 0 (tel quel) | -26.9143 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.7 | D14, R10 A.1 |
| C_N plateau (i), 8x8 (meV) | — | 0 (tel quel) | -14.4103 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.8 | D14, R10 A.1 |
| C_N plateau (i), 9x9 (meV) | — | 0 (tel quel) | -25.1437 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.9 | D14, R10 A.1 |
| C_N plateau (i), 10x10 (meV) | — | 0 (tel quel) | -9.3417 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.10 | D14, R10 A.1 |
| C_N plateau (i), 11x11 (meV) | — | 0 (tel quel) | -9.5507 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.11 | D14, R10 A.1 |
| C_N plateau (i), 12x12 (meV) | — | 0 (tel quel) | -18.6896 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.12 | D14, R10 A.1 |
| C_N plateau (i), 15x15 (meV) | — | 0 (tel quel) | -16.2185 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.13 | D14, R10 A.1 |
| C_N plateau (i), 18x18 (meV) | — | 0 (tel quel) | -15.0621 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.14 | D14, R10 A.1 |
| C_N plateau (i), 21x21 (meV) | — | 0 (tel quel) | -14.2343 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.15 | D14, R10 A.1 |
| C_N plateau (i), 24x24 (meV) | — | 0 (tel quel) | -13.4815 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.16 | D14, R10 A.1 |
| C_N plateau (i), 27x27 (meV) | — | 0 (tel quel) | -13.0058 | — | nouveau | `config/production.json (alignment)` ; `table_v2_plateau.md` l.17 | D14, R10 A.1 |
| Re M^L à K aligné (D6), 6x6 (eV) | — | — | +13.1744 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.37 | tab familles (D6) |
| Re M^L à K aligné (D6), 9x9 (eV) | — | — | +13.3417 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.40 | tab familles (D6) |
| Re M^L à K aligné (D6), 12x12 (eV) | — | — | +13.9133 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.43 | tab familles (D6) |
| Re M^L à K aligné (D6), 5x5 (eV) | — | — | +11.9733 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.46 | tab familles (D6) |
| Re M^L à K aligné (D6), 7x7 (eV) | — | — | +11.8045 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.49 | tab familles (D6) |
| Re M^L à K aligné (D6), 8x8 (eV) | — | — | +11.6565 | — | nouveau | `level2_families.csv (M_analysis.npz)` ; `table_v2_plateau.md` l.52 | tab familles (D6) |
| carte fig_M_map, Ṽ_π à K (eV Å²) : brut ; aligné | — | 107.22 | 107.22 ; 114.89 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.65 | fig_M_map (D6) |
| carte fig_M_map, Ṽ_π* à K (eV Å²) : brut ; aligné | — | 109.79 | 109.79 ; 117.62 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.66 | fig_M_map (D6) |
| carte fig_M_map, M^L_∥ à K (eV Å²) : brut ; aligné | — | 84.39 | 84.39 ; 92.09 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.67 | fig_M_map (D6) |
| porte D6 align_frozen_projector | — | — | 1.4e-10 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.68 | D6 |
| porte D6 align_gate_closed_vs_wannier | — | — | 2.9e-15 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.69 | D6 |
| porte D6 align_gate_coarse_wannier | — | — | 1.3e-15 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.70 | D6 |
| porte D6 align_gate_frozen_bloch | — | — | 1.4e-10 | — | nouveau | `M_analysis.npz` ; `table_v2_plateau.md` l.71 | D6 |
| localité 9x9 dense : abscisse max (a) | — | 22.517 | 15.588 | — | nouveau | `mwr_locality.npz` ; `table_v2_plateau.md` l.74 | fig_locality (P-b2) |
| localité 12x12 dense : abscisse max (a) | — | 20.785 | 13.856 | — | nouveau | `mwr_locality.npz` ; `table_v2_plateau.md` l.77 | fig_locality (P-b2) |
| Re Σ / Γ, 9x9, R_cut 0 : méd. Γ ; méd. Re Σ (meV) ; E_res | — | 3009.73 ; 17.67 ; -0.175 | 3013.96 ; 29.47 ; -0.175 | — | nouveau | `m_rcut_resigma.csv` ; `table_v2_plateau.md` l.92 | C10, C11 |
| Re Σ / Γ, 9x9, R_cut 1 : méd. Γ ; méd. Re Σ (meV) ; E_res | — | 3223.74 ; 682.54 ; -0.227 | 3273.25 ; 761.79 ; -0.227 | — | nouveau | `m_rcut_resigma.csv` ; `table_v2_plateau.md` l.93 | C10, C11 |
| Re Σ / Γ, 9x9, R_cut 2 : méd. Γ ; méd. Re Σ (meV) ; E_res | — | 3170.21 ; 479.23 ; -0.175 | 3222.08 ; 737.20 ; -0.175 | — | nouveau | `m_rcut_resigma.csv` ; `table_v2_plateau.md` l.94 | C10, C11 |
| Re Σ / Γ, 9x9, R_cut 3 : méd. Γ ; méd. Re Σ (meV) ; E_res | — | 3132.60 ; -26.65 ; -0.175 | 3189.01 ; 605.08 ; -0.175 | — | nouveau | `m_rcut_resigma.csv` ; `table_v2_plateau.md` l.95 | C10, C11 |
| Re Σ / Γ, 9x9, R_cut 4 : méd. Γ ; méd. Re Σ (meV) ; E_res | — | 3147.47 ; -410.29 ; -0.175 | 3200.00 ; 696.74 ; -0.175 | — | nouveau | `m_rcut_resigma.csv` ; `table_v2_plateau.md` l.96 | C10, C11 |
| résonance 6x6 : médiane courbe Γ_T ±3 eV (meV) | — | 3780.12 | 3713.28 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.118 | R6 3.4 (D9) |
| résonance 6x6 : médiane Γ_Born (meV) | — | 1.7433e+05 | 1.7855e+05 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.119 | R6 3.4 (D9) |
| résonance 6x6 : Born/T médian | — | 46.512 | 46.286 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.120 | R6 3.4 (D9) |
| résonance 6x6 : pic Γ_T (eV) | — | -0.240 | -0.180 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.121 | R6 3.4 (D9) |
| résonance 6x6 : pic Γ_T/ρ₀ | — | -0.175 | -0.175 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.122 | R6 3.4 (D9) |
| résonance 6x6 : pic δρ | — | -0.847 | -0.797 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.123 | R6 3.4 (D9) |
| résonance 6x6 : pic \|T̄\| | — | -0.177 | -0.172 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.124 | R6 3.4 (D9) |
| résonance 6x6 : pic −Im T̄ | — | -0.220 | -0.177 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.125 | R6 3.4 (D9) |
| résonance 6x6 : min \|Re T̄\| | — | 1.435 | 2.875 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.126 | R6 3.4 (D9) |
| résonance 6x6 : Re T̄(E_D) (eV) | — | 9.218 | 11.884 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.127 | R6 3.4 (D9) |
| résonance 6x6 : Im T̄(E_D) (eV) | — | -1.394 | -1.892 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.128 | R6 3.4 (D9) |
| résonance 6x6 : E_res états ±1,5 eV | — | -0.227 | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.129 | R6 3.4 (D9) |
| résonance 6x6 : médiane \|Γ_T\| états (meV) | — | 3156.45 | 3151.65 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.130 | R6 3.4 (D9) |
| résonance 6x6 : zéros de Re T̄ | — | -0.262, +1.433 | -2.202, -2.172, -2.157, -2.115, -2.110, -0.220, +2.693, +2.703, +2.755, +2.778, +2.818, +2.850, +2.875, +2.928, +2.943 | — | nouveau | `resonance_6x6.npz` ; `table_v2_plateau.md` l.131 | R6 3.4 (D9) |
| critère 6x6 (complet) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 1.432e-04 (-0.847) ; \|λ\| 0.0019 (-0.847) ; Friedel -1.0009 | 1.412e-04 (-0.795) ; \|λ\| 0.0019 (-0.797) ; Friedel -1.0004 | — | nouveau | `resonance_criteria_6x6.npz` ; `table_v2_plateau.md` l.132 | §2, C15, C16 |
| critère 6x6 (bloc π) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 3.299e-02 (-0.175) ; \|λ\| 0.4962 (-0.175) ; Friedel -0.9978 | 1.876e-02 (-0.172) ; \|λ\| 0.3832 (-0.127) ; Friedel -0.9981 | — | nouveau | `resonance_criteria_6x6.npz` ; `table_v2_plateau.md` l.133 | §2, C15, C16 |
| critère 6x6 (bloc σ) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 3.025e-04 (-0.847) ; \|λ\| 0.0019 (-0.847) ; Friedel -0.0030 | 4.579e-04 (-0.797) ; \|λ\| 0.0019 (-0.797) ; Friedel -0.0023 | — | nouveau | `resonance_criteria_6x6.npz` ; `table_v2_plateau.md` l.134 | §2, C15, C16 |
| résonance 12x12 : médiane courbe Γ_T ±3 eV (meV) | — | 3749.54 | 3823.41 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.135 | R6 3.4 (D9) |
| résonance 12x12 : médiane Γ_Born (meV) | — | 1.7906e+05 | 1.8067e+05 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.136 | R6 3.4 (D9) |
| résonance 12x12 : Born/T médian | — | 45.886 | 45.393 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.137 | R6 3.4 (D9) |
| résonance 12x12 : pic Γ_T (eV) | — | -0.180 | -0.175 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.138 | R6 3.4 (D9) |
| résonance 12x12 : pic Γ_T/ρ₀ | — | -0.170 | -0.130 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.139 | R6 3.4 (D9) |
| résonance 12x12 : pic δρ | — | -0.820 | -0.797 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.140 | R6 3.4 (D9) |
| résonance 12x12 : pic \|T̄\| | — | -0.170 | -0.127 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.141 | R6 3.4 (D9) |
| résonance 12x12 : pic −Im T̄ | — | -0.175 | -0.170 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.142 | R6 3.4 (D9) |
| résonance 12x12 : min \|Re T̄\| | — | -2.207 | -1.762 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.143 | R6 3.4 (D9) |
| résonance 12x12 : Re T̄(E_D) (eV) | — | 12.429 | 13.619 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.144 | R6 3.4 (D9) |
| résonance 12x12 : Im T̄(E_D) (eV) | — | -2.103 | -2.390 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.145 | R6 3.4 (D9) |
| résonance 12x12 : E_res états ±1,5 eV | — | -0.175 | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.146 | R6 3.4 (D9) |
| résonance 12x12 : médiane \|Γ_T\| états (meV) | — | 3166.57 | 3263.13 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.147 | R6 3.4 (D9) |
| résonance 12x12 : zéros de Re T̄ | — | -2.207, -2.170, -2.160, -0.210, +2.690, +2.705, +2.753, +2.778, +2.815, +2.850, +2.875, +2.928, +2.943 | -1.812, -1.785, -1.762, -0.197 | — | nouveau | `resonance_12x12.npz` ; `table_v2_plateau.md` l.148 | R6 3.4 (D9) |
| critère 12x12 (complet) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 1.443e-04 (-0.817) ; \|λ\| 0.0019 (-0.817) ; Friedel -1.0005 | 1.437e-04 (-0.797) ; \|λ\| 0.0019 (-0.797) ; Friedel -1.0003 | — | nouveau | `resonance_criteria_12x12.npz` ; `table_v2_plateau.md` l.149 | §2, C15, C16 |
| critère 12x12 (bloc π) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 1.826e-02 (-0.170) ; \|λ\| 0.3440 (-0.127) ; Friedel -0.9981 | 1.408e-02 (-0.127) ; \|λ\| 0.3108 (-0.127) ; Friedel -0.9982 | — | nouveau | `resonance_criteria_12x12.npz` ; `table_v2_plateau.md` l.150 | §2, C15, C16 |
| critère 12x12 (bloc σ) : min\|det\|/max ; min\|λ\| ; ∫δρ | — | 4.442e-04 (-0.817) ; \|λ\| 0.0019 (-0.817) ; Friedel -0.0024 | 5.307e-04 (-0.797) ; \|λ\| 0.0019 (-0.797) ; Friedel -0.0021 | — | nouveau | `resonance_criteria_12x12.npz` ; `table_v2_plateau.md` l.151 | §2, C15, C16 |
| C14 plateau (±9,05 meV, V_loc aligné) : E_res_states | — | — | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.160 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : E_res_states_shiftm9_05 | — | — | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.161 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : E_res_states_shiftp9_05 | — | — | -0.175 | — | retiré (décision) — E_res (R10 B.1) ; mécanique : nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.162 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : median_GT_states_meV | — | — | 3188.350 | — | nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.163 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : median_GT_states_shiftm9_05_meV | — | — | 3157.069 | — | nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.164 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : median_GT_states_shiftp9_05_meV | — | — | 3237.293 | — | nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.165 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : peak_GT_shiftm9_05 | — | — | -0.180 | — | nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.166 | C14 (D10) |
| C14 plateau (±9,05 meV, V_loc aligné) : peak_GT_shiftp9_05 | — | — | -0.180 | — | nouveau | `resonance_9x9_shiftL.npz` ; `table_v2_plateau.md` l.167 | C14 (D10) |
| tab:tests_M : convention intensive (cellule unitaire), famille 3m (N = 6, 9, 12) | — | — | 2.3e-02 (OK) | — | nouveau | `M_tests_summary.csv` ; `table_v2_plateau.md` l.185 | tab:tests_M (nouvelle ligne) |
| tab:tests_M : convention intensive (cellule unitaire), famille non-3m (N = 5, 7, 8) | — | — | 2.1e-02 (OK) | — | nouveau | `M_tests_summary.csv` ; `table_v2_plateau.md` l.186 | tab:tests_M (nouvelle ligne) |
| D.1 Kaasbjerg 9×9, valence \| K (eV Å²) : R9 brut → tel quel WS ; plateau WS | — | 76.955 | 76.956 ; 81.605 | — | nouveau | `c/c2_results.json` ; `table_v2_plateau.md` l.187 | carte de Kaasbjerg |
| D.1 Kaasbjerg 9×9, conduction \| K (eV Å²) : R9 brut → tel quel WS ; plateau WS | — | 78.159 | 78.168 ; 82.821 | — | nouveau | `c/c2_results.json` ; `table_v2_plateau.md` l.188 | carte de Kaasbjerg |
| D.1 Kaasbjerg 9×9, valence \| K' (eV Å²) : R9 brut → tel quel WS ; plateau WS | — | 81.176 | 81.177 ; 81.168 | — | nouveau | `c/c2_results.json` ; `table_v2_plateau.md` l.189 | carte de Kaasbjerg |
| D.1 Kaasbjerg 9×9, conduction \| K' (eV Å²) : R9 brut → tel quel WS ; plateau WS | — | 82.297 | 82.305 ; 82.297 | — | nouveau | `c/c2_results.json` ; `table_v2_plateau.md` l.190 | carte de Kaasbjerg |

## Écarts entre les deux tables sources (rapportés, non corrigés)

- R6 l.5 / R10 l.18 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.7 / R10 l.20 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.9 / R10 l.22 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.11 / R10 l.24 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.13 / R10 l.26 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.15 / R10 l.28 : fichier source différent (`level1_summary.csv / level2_summary.csv` vs `level1_summary.csv`)
- R6 l.52 / R10 l.78 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.53 / R10 l.79 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.54 / R10 l.80 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.55 / R10 l.81 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.56 / R10 l.85 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.57 / R10 l.86 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.58 / R10 l.87 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.59 / R10 l.88 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.60 / R10 l.82 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.61 / R10 l.83 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.62 / R10 l.84 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.63 / R10 l.89 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.64 / R10 l.90 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.65 / R10 l.91 : fichier source différent (`m_rcut_convergence.csv` vs `m_rcut_convergence.csv ; c/c2_results.json`)
- R6 l.87 v2 « 3132.872 » ≠ R10 v2 « 3132.87 »
- R6 l.88 v2 « -2.447, -0.222, 1.638 » ≠ R10 v2 « -2.447, -0.222, +1.638 »
- R6 l.90 : v1/v2 = λ_min complexe, final = |λ| (forme de la table R10)
- R6 l.90 v2 « +0.0000+0.0019i (-0.812) » ≠ R10 v2 « |λ| 0.0019 (-0.812) »
- R6 l.91 v2 « -1.0007 / -1.0007 / +0.669 » ≠ R10 v2 « -1.0007 »
- R6 l.93 : v1/v2 = λ_min complexe, final = |λ| (forme de la table R10)
- R6 l.93 v2 « +0.2186+0.3310i (-0.170) » ≠ R10 v2 « |λ| 0.3967 (-0.170) »
- R6 l.94 v2 « -0.9980 / -0.9980 / -0.285 » ≠ R10 v2 « -0.9980 »
- R6 l.96 : v1/v2 = λ_min complexe, final = |λ| (forme de la table R10)
- R6 l.96 v2 « +0.0000+0.0019i (-0.812) » ≠ R10 v2 « |λ| 0.0019 (-0.812) »
- R6 l.97 v2 « -0.0028 / -0.0028 / +0.955 » ≠ R10 v2 « -0.0028 »

## Chiffres nouveaux

Chaque section nomme la source (fichier : clé ou ligne), l'alignement et l'E_D utilisés.

### a. Porte A.2 (ΔV appliqué directement aux états de Bloch purs contre M/N_cells)

Alignement : aucun (M bruts, non alignés) ; E_D sans objet (éléments de matrice). Sources : R5 A.2 = table de `article/R5_base_vs_M/R5_rapport.md` (l.374–377 ; détail par paire dans `a/a_results.json` : A2.pure["16"/"20"]) ; R6 1.3 = `article/R6_production_corrigee/etape1/gate/gate_table.md` (json `gate_<S>_<niveau>.json`) ; R6 1.2 = `etape1/assemble_summary.jsonl`.

R5 A.2 (9×9, convention de production v1 = M^L en norme super-cellule) :

| ligne | partie | base 16 : max\|Δ\| | base 20 : max\|Δ\| | rapport direct / (M/81) |
|---|---|---|---|---|
| R5_rapport.md l.376 | non locale (M^NL) | 3,4e-15 eV | 1,5e-9 eV (bruit 1,3e-7 des k du XML dense) | 1,0000 |
| R5_rapport.md l.377 | locale (M^L) | 1,7e-1 eV | 1,5e-1 eV | **81,0000** sur toutes les paires non nulles (ex. (3, K, 3, K) : direct +0,155642, M^L/81 = +0,001922, M^L du fichier = +0,155642) |

R6 1.3, porte A.2 sur chaque M2 (seuil 1e-6 eV ; `gate_table.md`) :

| ligne | N | niveau | fichier M^L | k testés | paires | max\|Δ_L\| (eV) | max\|Δ_NL\| (eV) | direct/(M^L/N_cells) médian | direct/(M^NL/N_cells) médian | seuil | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gate_table.md l.3 | 5x5 | coarse | `M_L_5x5.npy` | 25 | 8 | 4.16e-15 | 2.91e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.4 | 5x5 | dense | `M_L_dense_5x5.npy` | 25 | 8 | 3.55e-15 | 9.44e-16 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.5 | 6x6 | coarse | `M_L_6x6.npy` | 36 | 8 | 3.30e-14 | 6.40e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.6 | 6x6 | dense | `M_L_dense_6x6.npy` | 36 | 8 | 3.45e-10 | 5.20e-11 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.7 | 7x7 | coarse | `M_L_7x7.npy` | 49 | 8 | 1.51e-14 | 8.08e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.8 | 7x7 | dense | `M_L_dense_7x7.npy` | 49 | 8 | 1.21e-08 | 4.04e-09 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.9 | 8x8 | coarse | `M_L_8x8.npy` | 64 | 8 | 1.32e-14 | 6.60e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.10 | 8x8 | dense | `M_L_dense_8x8.npy` | 64 | 8 | 4.66e-15 | 1.56e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.11 | 9x9 | coarse | `M_L_9x9_nb128.npy` | 81 | 8 | 3.52e-15 | 5.76e-16 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.12 | 9x9 | coarse | `M_L_9x9.npy` | 81 | 8 | 6.80e-15 | 3.41e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.13 | 9x9 | dense | `M_L_dense_9x9.npy` | 81 | 8 | 3.54e-09 | 1.54e-09 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.14 | 10x10 | coarse | `M_L_10x10.npy` | 100 | 8 | 4.34e-15 | 5.55e-16 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.15 | 12x12 | coarse | `M_L_12x12.npy` | 144 | 8 | 6.63e-15 | 3.10e-15 | 1.000000 | 1.000000 | 1e-06 | OK |
| gate_table.md l.16 | 12x12 | dense | `M_L_dense_12x12.npy` | 144 | 8 | 1.10e-10 | 3.19e-11 | 1.000000 | 1.000000 | 1e-06 | OK |

R6 1.2, réassemblage M2 = N_cells·M^L + M^NL (`assemble_summary.jsonl`, une ligne par fichier ; max|·| en Ha) :

| ligne | fichier | N_cells | max\|N_cells·M^L\| | max\|M^NL\| | max\|M2\| | herm. rel. | ok |
|---|---|---|---|---|---|---|---|
| l.1 | M_ed_7x7.npy | 49 | 0.6399768457059315 | 0.2698373101523188 | 0.8544400764646501 | 4.2898783938949e-15 | True |
| l.2 | M_ed_9x9.npy | 81 | 0.6922845683408538 | 0.2698317093078657 | 0.9114336258387719 | 5.927560996665594e-15 | True |
| l.3 | M_ed_10x10.npy | 100 | 0.6459919139709387 | 0.269836889966939 | 0.8556773923361867 | 9.298551955921207e-15 | True |
| l.4 | M_ed_12x12.npy | 144 | 0.7116008459012562 | 0.26983353337533345 | 0.9354881441758727 | 1.5190844270636046e-14 | True |
| l.5 | M_ed_5x5.npy | 25 | 0.6710135769598589 | 0.2697723824585996 | 0.92599739850567 | 2.4165575942945836e-16 | True |
| l.6 | M_ed_6x6.npy | 36 | 0.8025294325333105 | 0.26982037166914097 | 1.0034140447483177 | 2.23011234232543e-16 | True |
| l.7 | M_ed_8x8.npy | 64 | 0.651879043998663 | 0.26982879844776875 | 0.8934049274945324 | 3.728062126561142e-16 | True |
| l.8 | M_dense_5x5.npy | 25 | 0.687931838115321 | 0.28933195892156505 | 0.9632214812571466 | 3.720753972771275e-15 | True |
| l.9 | M_dense_6x6.npy | 36 | 0.7409449759274279 | 0.29128022763032035 | 1.0235396031267303 | 5.224558195398481e-15 | True |
| l.10 | M_dense_7x7.npy | 49 | 0.6952352524888566 | 0.29128038199925754 | 0.9733055371935604 | 7.52498223605507e-15 | True |
| l.11 | M_dense_8x8.npy | 64 | 0.7051909598590628 | 0.291279653698183 | 0.9825689097168625 | 1.0847367029740555e-14 | True |
| l.12 | M_dense_9x9.npy | 81 | 0.7590675752447206 | 0.2896060247190917 | 1.0425396184975904 | 1.0251565063180901e-14 | True |
| l.13 | M_dense_12x12.npy | 144 | 0.7780732024010126 | 0.2912802276303239 | 1.0612419906589767 | 2.035575461327524e-14 | True |
| l.14 | M_L_dense_9x9_coarsecheck.npy | 81 | 0.6922845683408542 | 0.0 | 0.0 | 0.0 | True |
| l.15 | M_ed_9x9_nb128.npy | 81 | 0.8315135771581877 | 0.37838048802797886 | 1.2098940651861665 | 7.915004667653365e-15 | True |

### b. Escalier D4 9×9 (états localisés, ε − E_D en eV ; (w₂) ; R6 2.1)

E_D par marche = clé `E_D` de chaque variante (`d4_results.json` : E_D.{SC_P, uc16_at_K, uc20dense_at_K, wannier_at_K}) ; « décalage rigide » = `rigid_shift_fit.shift_eV` (269 états < E_D − 4 eV, QE aligné site unique contre modèle sur son E_D). Portes de R6 2.1 : gate1 4.884981308350689e-14 eV, gate2 (k = m/27) 2.2595258997171186e-12 eV.

| marche | E_D (eV) | alignement | décalage rigide (meV) | σ (pairs) localisés | π (impairs) localisés | source |
|---|---|---|---|---|---|---|
| QE (D1, R6 2.1 ; E_D = E_D.SC_P, quadruplet de la parfaite) | -4.238471201574294 | site unique, Lu = -0.02465140127665677 eV (R5 C) | — | 0.1014 (0.713); 0.1014 (0.713) | -1.7595 (0.114); -1.7595 (0.114); -0.7373 (0.246); 0.2693 (0.102) | `d4_results.json` : variants.QE_D1 |
| (a1) M2 16 bandes, total | -4.238470462286799 | aucun (non aligné) | +24.3 (269 états) | aucun | -1.7751 (0.116); -1.7751 (0.116); -0.7274 (0.236); 0.2626 (0.107) | `d4_results.json` : variants.a1_tot |
| (a1) M2^L seul | -4.238470462286799 | aucun | — | -2.2406 (0.618) | -1.7751 (0.116); -1.7751 (0.116); -0.7419 (0.239); 0.2551 (0.105) | `d4_results.json` : variants.a1_L |
| (a1) M^NL seul | -4.238470462286799 | aucun | — | -2.8467 (0.255); -2.8467 (0.255) | -1.3701 (0.222) | `d4_results.json` : variants.a1_NL |
| (a2) M2 dense ⊂ 81 k, 16 b. | -4.23889512404274 | aucun | +24.4 (269 états) | aucun | -1.7751 (0.116); -1.7751 (0.116); -0.7274 (0.236); 0.2626 (0.107) | `d4_results.json` : variants.a2_16 |
| (a3) idem 20 b. | -4.23889512404274 | aucun | +24.4 (269 états) | aucun | -1.7756 (0.116); -1.7756 (0.116); -0.7344 (0.238); 0.2589 (0.106) | `d4_results.json` : variants.a3_20 |
| (b) 5 WF, V†εV ⊕ M2_W/81 | -4.238895124081603 | aucun | +24.0 (269 états) | -0.8120 (0.821) | -1.7736 (0.116); -1.7736 (0.116); -0.6886 (0.225); 0.2864 (0.115) | `d4_results.json` : variants.b_5wf |
| (c-all) H(R) + M2_W replié | -4.2388962954558025 | aucun | +24.0 (269 états) | -0.8120 (0.821) | -1.7736 (0.116); -1.7736 (0.116); -0.6886 (0.225); 0.2864 (0.115) | `d4_results.json` : variants.c_all |
| (c-3) idem, R_cut 3 | -4.2388962954558025 | aucun | +6.8 (269 états) | -0.8132 (0.821) | -1.7652 (0.115); -1.7627 (0.115); -0.6774 (0.227); 0.3024 (0.113) | `d4_results.json` : variants.c_3 |
| QE, alignement de Kumagai–Oba (R10 A.2) | -4.238471201574294 | C_9 = -0.02514371002267805 eV | — | 0.1019 (0.713); 0.1019 (0.713) | -0.7368 (0.246) | `R10_plateau/a/a2_results.json` : rows.9x9."R5 C".*.x_plateau |
| (c-3) « aligné » de R10 A.3 = modèle de R9 C.1 avec C_9 à site unique (-0.02465140127665677 eV), **pas Kumagai–Oba** ; aucune valeur K–O | -4.2388962954558025 | site unique (modèles de R9 C.1 non recalculés (aligné : C_9 = Lu de R9) ; QE réaligné avec C_N de R10 A.1) | — | -0.7886 | -0.6607 | `R10_plateau/a/a3_results.json` : rows.9x9.aligne ; `R9_controles/c/c_results.json` |
| (c-3) « brut » de R10 A.3 (= R6 (c-3), non aligné) | idem | aucun | — | -0.8132 | -0.6774 | `a3_results.json` : rows.9x9.brut |

### c. Convergence en bandes, 9×9, (a1) à n bandes (R6 2.2 ; R5 B.3)

Non aligné (M2 nb128 réassemblé, `M_ed_9x9_nb128`) ; E_D = -4.238470783590611 eV (`b3_results_M2.json` : E_D_128) ; premier n avec une paire σ à < 0,3 eV de +0,101 : 48 (v2) / None (v1, convention de production) / 48 (v1 diagnostique M^L × 81). QE : π −0,737 (0,246), σ +0,101 ×2 (0,713) (R6 2.2).

| n | dim | π v2 (ε − E_D ; w₂) | σ v2 | π v1 (production) | σ v1 | π v1 M^L × 81 |
|---|---|---|---|---|---|---|
| 16 | 1296 | -0.7274 (0.236) | aucun | -1.3496 (0.226) | -2.8117 (0.287); -2.8117 (0.287) | -0.7274 |
| 24 | 1944 | -0.7406 (0.240) | 0.8429 (0.674); 0.8429 (0.674) | -1.3676 (0.225) | -2.8671 (0.236); -2.8671 (0.236) | -0.7406 |
| 32 | 2592 | -0.7469 (0.242) | 0.4278 (0.697); 0.4278 (0.697) | -1.3786 (0.224) | -2.9014 (0.201); -2.9014 (0.201) | -0.7469 |
| 48 | 3888 | -0.7544 (0.244) | 0.3449 (0.700); 0.3449 (0.700) | -1.3969 (0.222) | -2.9087 (0.193); -2.9087 (0.193) | -0.7544 |
| 64 | 5184 | -0.7573 (0.245) | 0.2107 (0.708); 0.2107 (0.708) | -1.4080 (0.220) | -2.9216 (0.180); -2.9216 (0.180) | -0.7573 |
| 96 | 7776 | -0.7595 (0.245) | 0.1367 (0.713); 0.1367 (0.713) | -1.4226 (0.218) | -2.9331 (0.168); -2.9331 (0.168) | -0.7595 |
| 128 | 10368 | -0.7604 (0.245) | 0.1157 (0.714); 0.1157 (0.714) | -1.4314 (0.216) | -2.9400 (0.160); -2.9400 (0.160) | -0.7604 |

Sources : `R6_production_corrigee/etape2/b3/b3_results_M2.json` (pi_vs_n, per_n[n].blocks.even.localized) ; `R5_base_vs_M/b/b3_results.json` (pi_vs_n, per_n, pi_vs_n_Lx81).

### d. Niveaux QE π et σ par taille, site unique → Kumagai–Oba (R10 A.2)

x = ε − C − E_D (eV) ; x_plateau = x_Lu + Lu − C_N ; w₂ inchangé ; E_D par taille (quadruplet pour 3m, maille pour N ≠ 3m) ; famille lue dans `a/A2_tables.md`.

| N | famille | E_D (eV) | Lu (meV) | C_N (meV) | π bande | π : x site unique → x K–O | w₂ | σ bandes | σ : x site unique → x K–O | w₂ | source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5x5 | non-3m | -4.2313837467852125 | -122.33 | -57.5626 | 98 | -0.3284 → -0.3931 | 0.267 | 99, 100 | -0.1308 → -0.1955 / -0.1308 → -0.1955 | 0.702 / 0.702 | R5 C |
| 6x6 | 3m | -4.237452738557307 | -98.62 | -50.5091 | 140 | -1.0654 → -1.1135 | 0.284 | 144, 145 | 0.1560 → 0.1078 / 0.1560 → 0.1078 | 0.711 / 0.711 | R7 D1 |
| 7x7 | non-3m | -4.239039893859123 | -38.54 | -26.9143 | 194 | -0.4948 → -0.5064 | 0.301 | 195, 196 | -0.2823 → -0.2939 / -0.2823 → -0.2939 | 0.703 / 0.703 | R5 C |
| 8x8 | non-3m | -4.237684021596026 | -22.97 | -14.4103 | 254 | -0.3698 → -0.3784 | 0.223 | 255, 256 | -0.1954 → -0.2039 / -0.1953 → -0.2039 | 0.708 / 0.708 | R5 C |
| 9x9 | 3m | -4.238471201574294 | -24.65 | -25.1437 | 320 | -0.7373 → -0.7368 | 0.246 | 324, 325 | 0.1014 → 0.1019 / 0.1014 → 0.1019 | 0.713 / 0.713 | R7 D1 |
| 10x10 | non-3m | -4.238894813481761 | -8.22 | -9.3417 | 398 | -0.4224 → -0.4212 | 0.241 | 399, 400 | -0.2383 → -0.2372 / -0.2383 → -0.2372 | 0.708 / 0.708 | R5 C |
| 11x11 | non-3m | -4.238510933406554 | -12.16 | -9.5507 | 482 | -0.3319 → -0.3345 | 0.191 | 483, 484 | -0.1641 → -0.1668 / -0.1641 → -0.1667 | 0.707 / 0.708 | R5 C |
| 12x12 | 3m | -4.238688683752139 | -29.24 | -18.6896 | 572 | -0.5508 → -0.5614 | 0.223 | 576, 577 | 0.1140 → 0.1034 / 0.1140 → 0.1034 | 0.713 / 0.713 | R7 D1 |
| 15x15 | 3m | -4.238781478817844 | -15.64 | -16.2185 | 896 | -0.4593 → -0.4588 | 0.207 | 900, 901 | 0.1033 → 0.1039 / 0.1033 → 0.1039 | 0.713 / 0.713 | R7 D1 |
| 18x18 | 3m | -4.238836282712659 | -18.18 | -15.0621 | 1292 | -0.3881 → -0.3912 | 0.194 | 1297, 1298 | 0.1070 → 0.1038 / 0.1070 → 0.1038 | 0.712 / 0.712 | R7 D1 |
| 21x21 | 3m | -4.238864868363458 | -13.77 | -14.2343 | 1760 | -0.3437 → -0.3432 | 0.183 | 1765, 1766 | 0.1030 → 0.1035 / 0.1030 → 0.1035 | 0.711 / 0.710 | R7 D1 |
| 24x24 | 3m | -4.238883603865964 | -13.12 | -13.4815 | 2300 | -0.3074 → -0.3071 | 0.174 | 2305, 2306 | 0.1026 → 0.1030 / 0.1026 → 0.1030 | 0.712 / 0.712 | R7 D1 |
| 27x27 | 3m | -4.23889533619687 | -12.77 | -13.0058 | 2912 | -0.2785 → -0.2783 | 0.166 | 2917, 2918 | 0.1024 → 0.1027 / 0.1024 → 0.1027 | 0.712 / 0.712 | R7 D1 |

Source : `article/R10_plateau/a/a2_results.json` : rows[S][source].{E_D, Lu_eV, C_N_eV, pi, sigma}.

### e. Chaîne repliée − QE, N = 6…27 (R10 A.3)

Modèles de R9 C.1 non recalculés. **La colonne « aligné » utilise C_9 à site unique = -0.02465140127665677 eV (−24,65 meV), et non C_9 de Kumagai–Oba (−25,14 meV)** (`a3_results.json` : C9_model_label = « sans plateau : valeur Lu (1,0 Å) » ; note : « modèles de R9 C.1 non recalculés (aligné : C_9 = Lu de R9) ; QE réaligné avec C_N de R10 A.1 »). QE : E_D par taille (A.2), colonne « plateau » = Kumagai–Oba. Écarts en meV, chaîne − QE.

| N | QE π : site unique → K–O | modèle non aligné π | écart vs QE site unique / K–O | modèle aligné (site unique) π | écart | σ modèle non aligné / aligné | σ QE K–O |
|---|---|---|---|---|---|---|---|
| 6x6 | -1.0654 → -1.1135 | -1.0438 | +21.6 / +69.7 | -1.0202 | +45.2 / +93.3 | -0.8308 / -0.8061 | 0.1078 / 0.1078 |
| 9x9 | -0.7373 → -0.7368 | -0.6774 | +59.8 / +59.3 | -0.6607 | +76.5 / +76.0 | -0.8132 / -0.7886 | 0.1019 / 0.1019 |
| 12x12 | -0.5508 → -0.5614 | -0.5168 | +34.0 / +44.5 | -0.5022 | +48.6 / +59.2 | -0.8125 / -0.7879 | 0.1034 / 0.1034 |
| 15x15 | -0.4593 → -0.4588 | -0.4277 | +31.6 / +31.1 | -0.4140 | +45.4 / +44.8 | -0.8125 / -0.7879 | 0.1039 / 0.1039 |
| 18x18 | -0.3881 → -0.3912 | -0.3711 | +16.9 / +20.0 | -0.3580 | +30.1 / +33.2 | -0.8125 / -0.7879 | 0.1038 / 0.1038 |
| 21x21 | -0.3437 → -0.3432 | -0.3320 | +11.6 / +11.2 | -0.3193 | +24.4 / +23.9 | -0.8125 / -0.7879 | 0.1035 / 0.1035 |
| 24x24 | -0.3074 → -0.3071 | -0.3033 | +4.1 / +3.8 | -0.2909 | +16.5 / +16.2 | -0.8125 / -0.7879 | 0.1030 / 0.1030 |
| 27x27 | -0.2785 → -0.2783 | -0.2812 | -2.7 / -2.9 | -0.2692 | +9.4 / +9.1 | -0.8125 / -0.7879 | 0.1027 / 0.1027 |

Source : `article/R10_plateau/a/a3_results.json` : rows[S].{QE_pi_Lu, QE_pi_plateau, brut, aligne, QE_sigma_plateau}.

### f. E_res et pic de la courbe Γ_T contre la grille de sortie (R10 B.1 ; 9×9, R_cut 3, η 0,02, N_k^int 900, E_D Wannier)

E_res = argmax |Γ| des états à ±1,5 eV (grandeur retirée du ch. 4, R10 B.1 ; ici pour mémoire) ; pics de Γ_T au pas 5 meV (prod) et 2,5 meV (fine).

| variante | grille | états | E_res − E_D (eV) | couronne x (rang) | pic Γ_T 5 meV | pic Γ_T 2,5 meV | médiane états (meV) | C (eV) |
|---|---|---|---|---|---|---|---|---|
| non aligné | 120² | 10336 | -0.22689085781678653 | 3 (2) | -0.22999999999956344 | -0.22999999999956344 | 3148.67 | 0.0 |
| non aligné | 240² | 41266 | -0.20183959702885002 | 9 (5) | -0.17999999999955563 | -0.18249999999955602 | 3136.61 | 0.0 |
| non aligné | 480² | 165234 | -0.19992939997973025 | 37 (16) | -0.19999999999955875 | -0.20249999999955914 | 3138.81 | 0.0 |
| non aligné | 960² | 660782 | -0.19050323410507453 | 133 (47) | -0.19999999999955875 | -0.19749999999955836 | 3137.74 | 0.0 |
| final (Kumagai–Oba) | 120² | 10336 | -0.13441123100174046 | 1 (1) | -0.1349999999995486 | -0.1349999999995486 | 3182.76 | -0.02514371002267805 |
| final (Kumagai–Oba) | 240² | 41266 | -0.17484882850105077 | 7 (4) | -0.17999999999955563 | -0.17999999999955563 | 3191.51 | -0.02514371002267805 |
| final (Kumagai–Oba) | 480² | 165234 | -0.17484882850105077 | 28 (13) | -0.17999999999955563 | -0.17999999999955563 | 3189.40 | -0.02514371002267805 |
| final (Kumagai–Oba) | 960² | 660782 | -0.1721023616137609 | 109 (39) | -0.17999999999955563 | -0.17999999999955563 | 3190.31 | -0.02514371002267805 |

Source : `article/R10_plateau/b/b_results.json` : B1[brut|plateau][grille].{E_res, crown, peak_GT_prod, peak_GT_fine, median_GT_states_meV, C_eV}.

### g. Familles N mod 3 (R10 B.2, C.1, C.3 ; R9 R.2)

B.2 : ε − E_F de la cellule avec lacune (non aligné) ; occupations XML par état (Γ, poids 2) ; gap de la parfaite à Γ.

| N | famille | gap parfaite Γ (eV) | ΔE_F (meV) | π : ε − E_F ; occ. | σ : ε − E_F ; occ. (chacun) | C_N (meV) | max\|M\| dense (eV) |
|---|---|---|---|---|---|---|---|
| 5x5 | non-3m | 2.809257681448693 | -235.8 | -0.0912 ; 0.8778 | 0.1064 ; 0.0611 / 0.1064 ; 0.0611 | -57.56 | 23.849840956591645 |
| 6x6 | 3m | 1.748198181772409e-07 | -20.5 | -1.1625 ; 1.0000 | 0.0589 ; 0.1622 / 0.0589 ; 0.1622 | -50.51 | 25.147561946808214 |
| 7x7 | non-3m | 3.0030460087554256 | -866.1 | -0.0969 ; 0.9016 | 0.1156 ; 0.0492 / 0.1156 ; 0.0492 | -26.91 | 23.817548744360064 |
| 8x8 | non-3m | 1.9760988591478879 | -583.3 | -0.0817 ; 0.8346 | 0.0928 ; 0.0827 / 0.0928 ; 0.0827 | -14.41 | 24.31117565954618 |
| 9x9 | 3m | 7.672752122545035e-08 | -2.4 | -0.7784 ; 1.0000 | 0.0602 ; 0.1584 / 0.0602 ; 0.1584 | -25.14 | 25.348304266624794 |
| 10x10 | non-3m | 2.038864168938389 | -508.7 | -0.0858 ; 0.8534 | 0.0983 ; 0.0733 / 0.0983 ; 0.0733 | -9.34 | — |
| 11x11 | non-3m | 1.5092934299440066 | -209.3 | -0.0788 ; 0.8205 | 0.0890 ; 0.0897 / 0.0890 ; 0.0897 | -9.55 | — |
| 12x12 | 3m | 3.9651267158546943e-08 | 2.5 | -0.6015 ; 1.0000 | 0.0633 ; 0.1497 / 0.0633 ; 0.1497 | -18.69 | 25.72092524023282 |

Source : `article/R10_plateau/b/b_results.json` : B2.rows[S].

Médianes de niveau 1 par taille et par famille (final, 240², η 0,02, R_cut 3, N_k^int 300 ; `results/M2_plateau/level2_families.csv`, lignes dans l'ordre du fichier) :

| ligne | N | famille | sous-réseau | médiane Γ·N_cells (meV) | p_z–p_z lacune (eV) | Re M^L à K | Re M^NL à K | Re M^L à K aligné (D6) |
|---|---|---|---|---|---|---|---|---|
| l.2 | 6x6 | 3m | B | 3149.15 | 31.552 | +11.3561 | +3.1030 | +13.1744 |
| l.3 | 9x9 | 3m | A | 3189.01 | 31.546 | +11.3051 | +3.1030 | +13.3417 |
| l.4 | 12x12 | 3m | B | 3264.51 | 31.547 | +11.2220 | +3.1030 | +13.9133 |
| l.5 | 5x5 | non-3m | A | 2997.14 | 31.060 | +10.5343 | +3.1004 | +11.9733 |
| l.6 | 7x7 | non-3m | A | 2942.85 | 30.963 | +10.4857 | +3.1009 | +11.8045 |
| l.7 | 8x8 | non-3m | A | 2974.18 | 31.107 | +10.7343 | +3.1009 | +11.6565 |

**Test C.3 de R10** (son vrai nom : « convention intensive (cellule unitaire) », (max − min)/moyenne de max|M| des M denses complets, bandes 1–16, par famille ; `c/c3_results.json` : tests_M) :

| famille | tailles | écart relatif (valeur pleine) | csv (`M_tests_summary.csv`) | seuil | verdict |
|---|---|---|---|---|---|
| 3m | 6x6, 9x9, 12x12 | 0.022568384833144092 | 2.3e-02 | 5e-2 | OK |
| non-3m | 5x5, 7x7, 8x8 | 0.02057391305510878 | 2.1e-02 | 5e-2 | OK |

**Test R9 R.2** (son vrai nom : `offdiag_intensive`, max|M| hors k = k′, trois variantes tel quel / site unique / Kumagai–Oba ; table de `article/R9_controles/R9_rapport.md` l.733–737, valeurs de la clôture de R9) :

| ligne | ensemble | taille | tel quel | Lu | plateau | mode |
|---|---|---|---|---|---|---|
| R9_rapport.md l.735 | grossiers | 5×5 / 7×7 / 8×8 / 9×9 | 25,198 / 23,250 / 24,311 / 24,801 | idem | idem | identique par construction (a) |
| R9_rapport.md l.736 | denses | 5×5 | 23,850 | 25,546 | 24,038 | exact : effet mesuré |
| R9_rapport.md l.737 | denses | 6×6 / 7×7 / 8×8 / 12×12 | 25,148 / 23,818 / 24,311 / 25,721 | idem | idem | identique par construction (b) |

(R9 l.739 : « (max − min)/moyenne : grossiers 7,98e-2 (trois variantes) ; denses six tailles 7,71e-2 / 7,61e-2 / 7,69e-2 ».)

### h. Kaasbjerg (R10 C.2 : Fig. 3 ; R8 : Fig. 13–14)

D.1 (9×9, k = K + δx̂, disques de 0,1471 Å⁻¹ sur la carte 240², eV Å², moyenne [min ; max]) ; **K′ en tête** (insensible à l'alignement). E_D Wannier ; variantes = alignement de M_W.

| variante | valence \| K′ | conduction \| K′ | valence \| K | conduction \| K |
|---|---|---|---|---|
| non aligné + Wigner-Seitz (final tel quel) | 81.177 [80.48 ; 82.33] | 82.305 [81.52 ; 83.71] | 76.956 [75.51 ; 78.73] | 78.168 [76.16 ; 80.36] |
| final (Kumagai–Oba) + Wigner-Seitz | 81.168 [80.53 ; 82.31] | 82.297 [81.52 ; 83.69] | 81.605 [77.40 ; 86.36] | 82.821 [77.81 ; 87.47] |
| R9 brut (étiquettes brutes) | 81.176 [80.49 ; 82.31] | 82.297 [81.50 ; 83.70] | 76.955 [75.51 ; 78.67] | 78.159 [76.14 ; 80.37] |
| R9 exact (F_W, C = site unique) | 81.168 [80.53 ; 82.28] | 82.289 [81.50 ; 83.68] | 81.512 [77.36 ; 86.15] | 82.706 [77.79 ; 87.32] |

D.2 (bloc 2×2 π/π* à (K, K), eV Å²) :

| bloc | valeurs propres | ½ Tr | normes de ligne |
|---|---|---|---|
| tot \| tel quel | -1.71, 153.45 | 75.87 | 107.22, 109.79 |
| tot \| plateau | 9.02, 164.17 | 86.59 | 114.89, 117.62 |
| L \| tel quel | -1.71, 120.77 | 59.53 | 84.39, 86.41 |
| L \| plateau | 9.02, 131.49 | 70.25 | 92.10, 94.29 |
| NL \| tel quel | -0.00, 32.68 | 16.34 | 22.83, 23.38 |

Source : `article/R10_plateau/c/c2_results.json` : kaasbjerg.{D1, R9_D1, D2}.

Fig. 13 (DOS, états/eV/maille/spin, ρ − ρ₀ sur [−1, 0] eV, grille 600², η_G 50 meV, N_k^int 900, η_t 20 meV ; R8 étape 2) et Fig. 14 (A_k, η_G 25 meV, R8 étape 3) ; alignement = variante de V_loc, E_D Wannier :

| variante | c_i | position du max de ρ − ρ₀ (eV) | hauteur | largeur à mi-hauteur (eV) | ρ(E_D) | maxima de A_K (ε ; hauteur) | gap à K (eV) |
|---|---|---|---|---|---|---|---|
| non aligné | 0.1 % | -0.1624999999995529 | 0.002828839970526133 | 0.19063160286222433 | 0.009487892679800856 | 0.0091 (133.2) | — (un seul maximum) |
| non aligné | 1 % | -0.16499999999955328 | 0.026762869150290886 | 0.20124621847661056 | 0.013130811303999448 | 0.0076 (82.5); 0.1252 (42.2) | 0.11762346994893341 |
| final (Kumagai–Oba) | 0.1 % | -0.14749999999955055 | 0.0031151430321301986 | 0.1751883201569013 | 0.009540290760435992 | 0.0101 (130.6) | — (un seul maximum) |
| final (Kumagai–Oba) | 1 % | -0.14999999999955094 | 0.029087584195899296 | 0.1870328301810967 | 0.014166264211068568 | 0.0149 (82.4); 0.1331 (40.0) | 0.11824139676128233 |
| η_t = η_G (R8 `eta_unique` ; **signalée, non utilisée**) | 0.1 % | -0.1599999999995525 | 0.002227542964562474 | 0.24102480973834098 | 0.009734780412207883 | 0.0087 (132.6) | — |
| η_t = η_G (R8 `eta_unique` ; **signalée, non utilisée**) | 1 % | -0.1624999999995529 | 0.02132320519000496 | 0.25158032534443464 | 0.015233765435631843 | 0.0076 (82.8); 0.1253 (40.9) | 0.11774597002664208 |
| Kaasbjerg, valeurs lues (R8 7a) | 0.1 % | -0.2 | 0.0012561614461217945 | 0.16269441405413548 | 0.0054607006870102965 | colonne K saturée, non lisible | — |
| Kaasbjerg, valeurs lues (R8 7a) | 1.0 % | -0.21 | 0.030746893024432773 | 0.30635758813655406 | 0.012780548996175641 | 0.0100; 0.1100 (points blancs, colonne K) | 0,100 (R8_rapport.md, 7a) |

Sources : `article/R8_kaasbjerg/out/dos/dos_results.json` : metrics[<variante>_600_c<c>] ; `out/spec/spec_results.json` : runs[<variante>_c<c>].{K_two_highest, K_gap_eV} ; `out/7a/7a_results.json` : fig13_metrics[c], fig14["1.0"].dots (col 0) ; gap 0,100 : `article/R8_kaasbjerg/R8_rapport.md`, section 7a.

### i. C14 à ±9,05 meV (R10 C.1 ; 9×9, V_loc aligné Kumagai–Oba, ±rms du plateau)

E_D Wannier ; `results/M2_plateau/resonance_9x9_shiftL.npz` :

| clé | valeur |
|---|---|
| median_GT_states_meV | 3188.3501861192167 |
| median_GT_states_shiftp9_05_meV | 3237.2932737202404 |
| median_GT_states_shiftm9_05_meV | 3157.0691425783234 |
| E_res_states | -0.17484882850105077 |
| E_res_states_shiftp9_05 | -0.17484882850105077 |
| E_res_states_shiftm9_05 | -0.17484882850105077 |
| peak_GT | -0.17999999999955563 |
| peak_GT_shiftp9_05 | -0.17999999999955563 |
| peak_GT_shiftm9_05 | -0.17999999999955563 |
| shifts_meV | 9.05, -9.05 |
| C_N_eV | -0.02514371002267805 |

Pour mémoire (retiré) : C14 à ±25 meV autour du V_loc non aligné, `results/M2/resonance_9x9_shiftL.npz` :

| clé | valeur |
|---|---|
| median_GT_states_meV | 3132.8722520020106 |
| median_GT_states_shiftp25_meV | 3187.963814142512 |
| median_GT_states_shiftm25_meV | 3158.539299897675 |
| E_res_states | -0.17484882850105077 |
| E_res_states_shiftp25 | -0.17484882850105077 |
| E_res_states_shiftm25 | -0.22689085781678653 |
| peak_GT_shiftp25 | -0.17999999999955563 |
| peak_GT_shiftm25 | -0.23499999999956422 |
| shifts_meV | 25.0, -25.0 |
| ML_diag_mean | 5.427317937284203 |

### j. Re Σ médian contre R_cut (9×9, 240², η 0,02, N_k^int 300, fenêtre ±3 eV ; `m_rcut_resigma.csv`)

Colonnes arithmétiques : Δ = valeur(R_cut) − valeur(R_cut − 1) ; C_9 × Δn = -25.1437 meV × (n_L(R_cut) − n_L(R_cut − 1)), n_L = 1, 5, 13, 29, 49.

| R_cut | n_L | Re Σ méd. non aligné (meV) | Re Σ méd. final (meV) | Δ non aligné | Δ final | Δn_L | C_9 × Δn_L (meV) | Γ méd. non aligné (meV) | Γ méd. final (meV) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 17.67022542112226 | 29.46609782315468 | — | — | — | — | 3009.7275806072057 | 3013.956173409821 |
| 1 | 5 | 682.5387958904906 | 761.7874642951455 | +664.8686 | +732.3214 | 4 | -100.5748 | 3223.7370815896707 | 3273.251518746341 |
| 2 | 13 | 479.2299680161756 | 737.1973345558654 | -203.3088 | -24.5901 | 8 | -201.1497 | 3170.2095343835294 | 3222.0766437880466 |
| 3 | 29 | -26.6485679201531 | 605.0832702757496 | -505.8785 | -132.1141 | 16 | -402.2994 | 3132.596829713962 | 3189.013588325133 |
| 4 | 49 | -410.2885036414499 | 696.7410369288754 | -383.6399 | +91.6578 | 20 | -502.8742 | 3147.4740341233364 | 3199.9967487435815 |

Sources : `results/M2/m_rcut_resigma.csv` et `results/M2_plateau/m_rcut_resigma.csv` : med_ReSigma_meV, med_Gamma_meV, nL ; C_9 : `config/production.json` alignment.C_N_eV.9x9.

### k. ⟨ΔV⟩_3D et profil en z (R10 C.6)

Colonne arithmétique : N² × part |z| ≤ 3 Å. Aucun alignement (potentiels bruts).

| N | ⟨ΔV⟩_3D (meV) | N²⟨ΔV⟩ (meV) | part \|z\| ≤ 3 Å (meV) | N² × part | vide \|z\| > 5 Å : moyenne ± écart-type | bord de cellule (meV) | c (Å) | nz |
|---|---|---|---|---|---|---|---|---|
| 5x5 | 28.618369670298083 | 715.459241757452 | 34.081027052464876 | 852.0257 | -10.876 ± 9.756 | -13.295663946552835 | 15.87531632709 | 192 |
| 6x6 | 20.82580554609759 | 749.7289996595132 | 27.671068917816314 | 996.1585 | -16.183 ± 20.406 | -59.58921190476218 | 15.87531632709 | 192 |
| 7x7 | 14.553554054275692 | 713.1241486595089 | 17.137030623778365 | 839.7145 | -4.941 ± 4.285 | -3.823126428454459 | 15.87531632709 | 192 |
| 8x8 | 11.229852574145752 | 718.7105647453282 | 13.794013711479277 | 882.8169 | -5.344 ± 5.340 | -11.970974921093305 | 15.87531632709 | 192 |
| 9x9 | 9.087407856160715 | 736.0800363490179 | 12.1357332414816 | 982.9944 | -7.007 ± 8.667 | -25.497524846732713 | 15.87531632709 | 192 |
| 10x10 | 7.173000097175685 | 717.3000097175686 | 8.72367231394369 | 872.3672 | -3.170 ± 3.047 | -6.165771094898824 | 15.87531632709 | 192 |
| 11x11 | 5.946608321612843 | 719.5396069151541 | 7.365901794668849 | 891.2741 | -2.996 ± 3.076 | -7.334288429633597 | 15.87531632709 | 192 |
| 12x12 | 5.058902986927519 | 728.4820301175628 | 6.713204165069737 | 966.7014 | -3.732 ± 4.518 | -13.186695636690152 | 15.87531632709 | 192 |
| 15x15 | 1.4157659381566785 | 318.54733608525265 | 4.232656506705147 | 952.3477 | -7.122 ± 8.585 | -25.12348373125766 | 15.87531632709 | 192 |
| 18x18 | -8.481730789349267 | -2748.080775749163 | 2.9047627911496976 | 941.1431 | -30.426 ± 22.533 | -64.98992783897768 | 15.87531632709 | 192 |
| 21x21 | 1.033712816714808 | 455.8673521712303 | 2.1145810034463883 | 932.5302 | -2.577 ± 1.406 | -3.966487122453761 | 15.87531632709 | 192 |
| 24x24 | -62.73899109959698 | -36137.65887336786 | 1.6052656547074151 | 924.6330 | -167.482 ± 36.973 | -211.89472728983557 | 15.87531632709 | 192 |
| 27x27 | -99.91305943910427 | -72836.62033110701 | 1.2607314558298595 | 919.0732 | -252.913 ± 42.283 | -329.66893185965324 | 15.87531632709 | 192 |

Source : `article/R10_plateau/c/c6_results.json` : sizes[S].{mean3d_meV, N2_mean3d_meV, contributions["|z|<=3A"].sum_over_nz_meV, vacuum_mean_meV, vacuum_std_meV, at_cell_edge_meV, c_A, nz}.

### l. Chapitre 5 : Γ^ed/Γ^ep et Γ^ed médians (R10 C.1 ; c = 1 %, 300 K, ±3 eV, eV pour Γ)

| variante | médiane Γ^ed/Γ^ep | min (à ε − E_D) | max (à ε − E_D) | croisements | médiane Γ^ed (eV) | médiane Γ^ep (eV) | c | T |
|---|---|---|---|---|---|---|---|---|
| non aligné | 0.5368254182139114 | 0.12277441935151594 (à 1.8749999999998854) | 32.44212653751012 (à -0.18000000000007077) | -1.5035558608060515, 0.6896599672545027 | 0.037461534592945725 | 0.05695597837361438 | 0.01 | 300.0 |
| final | 0.5538362001212567 | 0.12819644290886822 (à 1.8749999999998854) | 34.79724200632459 (à -0.18000000000007077) | -1.4857070118600786, 0.7367106232088095 | 0.037730394206289 | 0.05695597837361438 | 0.01 | 300.0 |

Médianes sur ±1,2 eV (pas de clé scalaire dans le npz) : table R10 l.195 : « médiane Γ^ed, c = 1 % (meV) : ±3 eV / ±1,2 eV | 37.462 / 59.298 | 37.730 / 57.508 ».

Sources : `results/M2/ed_vs_ep_24k24q_mv0.02.npz`, `results/M2_plateau/ed_vs_ep_24k24q_mv0.02.npz` : median, min, x_min, max, x_max, crossings, median_ed, median_ep.

## Grandeurs retirées (décision), pour la traçabilité

1. **E_res à toutes les tailles** (R10 B.1) — `level1_summary.csv` (R_cut 3, 240², η 0,02), clé argmax_E_minus_ED_eV :

| N | non aligné | final |
|---|---|---|
| 5x5 | -0.2843 | -0.2268 |
| 6x6 | -0.2270 | -0.1749 |
| 7x7 | -0.2694 | -0.2694 |
| 8x8 | -0.2695 | -0.2270 |
| 9x9 | -0.1748 | -0.1748 |
| 12x12 | -0.1749 | -0.1749 |

2. **ε_∞ des extrapolations 1/N et 1/N²** (R9 C) — `article/R9_controles/c/c_results.json` : fits, fits_9_27 :

| série \| grandeur \| loi | ε_∞ (eV) | a | rms |
|---|---|---|---|
| brut \| quasi_bound \| 1/N^1 | -0.04855399082495497 | -5.847698644952903 | 0.014794805189120198 |
| brut \| quasi_bound \| 1/N^2 | -0.28045324756573903 | -28.51044656297081 | 0.029782585734803274 |
| brut \| ldos \| 1/N^1 | -0.44974720292629217 | 1.1325690676895381 | 0.1777997397060839 |
| brut \| ldos \| 1/N^2 | -0.4319246224576739 | 9.135545776873798 | 0.16708433597263403 |
| aligne \| quasi_bound \| 1/N^1 | -0.04020128769211526 | -5.762006102100654 | 0.014131978882239986 |
| aligne \| quasi_bound \| 1/N^2 | -0.2687755653946629 | -28.08287685593415 | 0.029781981976409792 |
| aligne \| ldos \| 1/N^1 | -0.04930015866526845 | -5.606601519619329 | 0.010286438016682655 |
| aligne \| ldos \| 1/N^2 | -0.2722915807227312 | -27.24784261589469 | 0.03246622731590943 |
| QE \| quasi_bound \| 1/N^1 | -0.05231759087455246 | -6.088698372730029 | 0.004541593662950508 |
| QE \| quasi_bound \| 1/N^2 | -0.29632885785734364 | -29.344653186550218 | 0.046454481053102806 |
| brut \| quasi_bound \| 1/N^1 (9…27) | -0.07810551598433148 | -5.3340517911128495 | 0.0048438153276139095 |
| brut \| quasi_bound \| 1/N^2 (9…27) | -0.25139157059363304 | -35.71596811068065 | 0.01411137598838742 |
| brut \| ldos \| 1/N^1 (9…27) | -0.07679109420811825 | -5.349929917738805 | 0.004650796814841861 |
| brut \| ldos \| 1/N^2 (9…27) | -0.250637058725524 | -35.81270163698375 | 0.014399866476982066 |
| aligne \| quasi_bound \| 1/N^1 (9…27) | -0.06834898341715127 | -5.27275976797619 | 0.0047514582950173925 |
| aligne \| quasi_bound \| 1/N^2 (9…27) | -0.23965024921573969 | -35.30417704148523 | 0.013982788266581177 |
| aligne \| ldos \| 1/N^1 (9…27) | -0.06897253520029585 | -5.264668091132053 | 0.004524028769670403 |
| aligne \| ldos \| 1/N^2 (9…27) | -0.24005931789369225 | -35.2394758781117 | 0.014234231561132248 |
| QE \| quasi_bound \| 1/N^1 (9…27) | -0.049894680606526784 | -6.130811943166623 | 0.004695611347209205 |
| QE \| quasi_bound \| 1/N^2 (9…27) | -0.2494673115732318 | -40.96345556444689 | 0.018662607245673354 |

3. **Valeurs d'alignement à site unique** (R10 A) — `a1_results.json` : sizes[S].Lu.{published_atom_1based, published_meV, dist_true_A}, far_true.mean10_meV, C_N_eV :

| N | atome (1-based) | site unique (meV) | source publiée | distance vraie (Å) | site vraiment le plus loin (meV) | C_N K–O (meV) |
|---|---|---|---|---|---|---|
| 5x5 | 49 | -122.32935833955594 | R5 C | 6.205591091537014 | 176.14060112116903 | -57.562606201523394 |
| 6x6 | 71 | -98.61571051991547 | R5 C | 7.3975555757151765 | -93.56163748696034 | -50.50914517629934 |
| 7x7 | 97 | -38.538818440962075 | R5 C | 8.65978712922328 | -37.03801350504676 | -26.914294254341268 |
| 8x8 | 1 | -22.965857127417166 | R5 C | 9.863407434286902 | 44.784226121549864 | -14.410274533323342 |
| 9x9 | 161 | -24.65140127665677 | R5 C | 11.119141883716845 | -45.88049389134241 | -25.14371002267805 |
| 10x10 | 198 | -8.223781758715631 | R5 C | 12.41118218307403 | -14.285155951668571 | -9.341716029942202 |
| 11x11 | 241 | -12.163021491289783 | R5 C | 13.580853070113466 | 8.020268283374321 | -9.550651557225422 |
| 12x12 | 287 | -29.235084249116028 | R5 C | 14.795111151430351 | -30.613966949193383 | -18.689593100831818 |
| 15x15 | 449 | -15.642378585788208 | R7 D1 | 18.507583046238242 | -23.337763030777836 | -16.21854502926427 |
| 18x18 | 1 | -18.18048243260506 | R7 D1 | 22.19266672714553 | -19.315705038259523 | -15.062070057554946 |
| 21x21 | 881 | -13.771481369850846 | R7 D1 | 25.90122779323973 | -16.87260601901386 | -14.234270413329362 |
| 24x24 | 2 | -13.1200424148048 | R7 D1 | 29.624450448890425 | -15.36179017854522 | -13.481540317473474 |
| 27x27 | 1457 | -12.765946493296099 | R7 D1 | 33.296609874798534 | -14.357374051972016 | -13.005771752253901 |

4. **C14 à 67 meV et à ±25 meV** : table principale, R6 l.75 et l.117–124 (statut « retiré ») ; remplacé par C14 à ±9,05 meV (section i).
5. **Zéro de det à −2,530 eV présenté comme état π** : table principale, R6 l.89–90 (v1) ; R4 D3 (`article/R4_quasi_lie/R4_rapport.md` l.669–675) : vecteur propre 100 % bloc σ, doublet E des trois sp² de l'atome retiré ; remplacé par le pôle σ du critère par bloc (final −0,785/−0,787 eV).
6. **Test « convention intensive » six tailles** : table principale, R6 l.112 / R10 l.180 (v1 2,7e-2, v2 7,7e-2 « À VOIR ») ; remplacé par les lignes par famille (section g).

## Compléments (partie 5, 2026-09-30)

Mêmes règles : valeurs copiées ou lues sur les fichiers de `results/` ; colonnes arithmétiques seulement là où c'est dit.

### 5.1 tab:rcut_M du mémoire : π–π* (`sv_mismatch_pi_blocks`) et diag (`diag_max_dM_over_max`), en %

Sources : `results/M/m_rcut_convergence.csv` (v1), `results/M2/m_rcut_convergence.csv` (non aligné), `results/M2_plateau/m_rcut_convergence.csv` (final) ; grille fine 60² ; % = fraction du csv × 100 (décalage de la virgule, chiffres du csv conservés).

| taille | R_cut | cellules | π–π* v1 | π–π* non aligné | π–π* final | diag v1 | diag non aligné | diag final |
|---|---|---|---|---|---|---|---|---|
| 9x9 | 0 | 1 | 40.851 | 31.515 | 31.708 | 40.851 | 25.913 | 31.708 |
| 9x9 | 1 | 5 | 23.589 | 18.945 | 18.508 | 17.113 | 22.987 | 15.005 |
| 9x9 | 2 | 13 | 7.1884 | 13.431 | 6.6075 | 4.8084 | 12.302 | 4.6200 |
| 9x9 | 3 | 29 | 2.5718 | 8.8431 | 2.5737 | 2.0915 | 7.2103 | 1.7074 |
| 9x9 | 4 | 49 | 1.5807 | 5.7559 | 1.9846 | 1.3997 | 4.7705 | 1.4106 |
| 9x9 | 5 | 81 | 1.0844 | 1.3976 | 1.1040 | 1.0836 | 1.5671 | 1.1945 |
| 9x9 | 6 | 113 | 0.99708 | 0.85302 | 0.91728 | 0.96235 | 0.96803 | 1.0052 |
| 12x12 | 0 | 1 | 97.203 | 99.342 | 99.271 | 72.972 | 85.345 | 86.782 |
| 12x12 | 1 | 5 | 18.404 | 19.948 | 18.948 | 18.245 | 24.005 | 18.925 |
| 12x12 | 2 | 13 | 11.813 | 14.949 | 12.521 | 11.480 | 13.805 | 12.506 |
| 12x12 | 3 | 29 | 3.0806 | 12.429 | 3.2908 | 2.3921 | 9.9701 | 3.0979 |
| 12x12 | 4 | 49 | 1.6892 | 10.364 | 2.5296 | 1.4996 | 8.3664 | 1.5291 |
| 12x12 | 5 | 81 | 1.2730 | 6.9336 | 1.6202 | 1.2249 | 5.9197 | 0.85068 |
| 12x12 | 6 | 113 | 1.1745 | 3.8144 | 1.1178 | 1.1327 | 3.6342 | 0.86769 |

Lignes en double : v1 : 14 lignes, couples (taille, R_cut) répétés : aucun ; non aligné : 14 lignes, couples (taille, R_cut) répétés : aucun ; final : 14 lignes, couples (taille, R_cut) répétés : aucun. **Aucune ligne en double dans les trois csv** (le csv v1 est seulement écrit en deux blocs : 9×9 et 12×12 R_cut 0–3, puis 9×9 et 12×12 R_cut 4–6). Colonnes égales dans le csv final (`max_dM_over_maxM` = `diag_max_dM_over_max`) : 12x12 R_cut 0, 12x12 R_cut 1, 12x12 R_cut 2, 12x12 R_cut 3, 9x9 R_cut 0 ; `diag_max_dM_over_max` = `diag_abs_mismatch` sur toutes les lignes des trois csv.

### 5.2 Tableau niveau 1 du mémoire (tab:convergence_gamma, `memoire/défauts.tex` l.664–670) : médiane Γ·N_cells (meV), 240², η 0,02, N_k^int 300

Sources : `level1_summary.csv` de `results/M` (v1), `results/M2` (non aligné), `results/M2_plateau` (final), colonne `median_Gamma_Ncells_meV` ; ordre des tailles = celui du mémoire.

| taille | N mod 3 | R_cut | v1 | non aligné | final | cité dans le mémoire |
|---|---|---|---|---|---|---|
| 6x6 | 0 | 0 | 2611.4029 | 2995.9157 | 3011.1803 | oui |
| 6x6 | 0 | 1 | 2536.6293 | 3155.8267 | 3242.6703 | oui |
| 6x6 | 0 | 2 | 2511.8726 | 3109.0882 | 3186.3363 | oui |
| 6x6 | 0 | 3 | 2509.2036 | 3155.9796 | 3149.1500 | oui |
| 6x6 | 0 | 4 | — | 3206.6405 | 3161.1535 | non (« --- » dans le mémoire) |
| 9x9 | 0 | 0 | 2596.9154 | 3009.7276 | 3013.9562 | oui |
| 9x9 | 0 | 1 | 2488.8991 | 3223.7371 | 3273.2515 | oui |
| 9x9 | 0 | 2 | 2466.8378 | 3170.2095 | 3222.0766 | oui |
| 9x9 | 0 | 3 | 2473.5466 | 3132.5968 | 3189.0136 | oui |
| 9x9 | 0 | 4 | — | 3147.4740 | 3199.9967 | oui |
| 12x12 | 0 | 0 | 2585.1923 | 3012.7732 | 3016.8913 | oui |
| 12x12 | 0 | 1 | 2475.3196 | 3262.2761 | 3299.4530 | oui |
| 12x12 | 0 | 2 | 2453.8468 | 3208.2824 | 3268.5115 | oui |
| 12x12 | 0 | 3 | 2458.7355 | 3162.7806 | 3264.5060 | oui |
| 12x12 | 0 | 4 | — | 3161.4916 | 3283.1683 | non (« --- » dans le mémoire) |
| 5x5 | 2 | 0 | 2627.5712 | 2924.3692 | 2936.1437 | oui |
| 5x5 | 2 | 1 | 2553.6336 | 2950.6992 | 2966.1883 | oui |
| 5x5 | 2 | 2 | 2534.0139 | 3089.9575 | 2967.5634 | oui |
| 5x5 | 2 | 3 | 2524.3006 | 3288.8973 | 2997.1354 | oui |
| 5x5 | 2 | 4 | — | 3302.7391 | 2997.0730 | non (« --- » dans le mémoire) |
| 7x7 | 1 | 0 | 2602.9720 | 2925.5822 | 2929.2014 | oui |
| 7x7 | 1 | 1 | 2508.3428 | 2935.0598 | 2947.9510 | oui |
| 7x7 | 1 | 2 | 2492.0954 | 2982.4421 | 2946.2803 | oui |
| 7x7 | 1 | 3 | 2487.4561 | 3052.0301 | 2942.8454 | oui |
| 7x7 | 1 | 4 | — | 3146.1805 | 2953.2845 | non (« --- » dans le mémoire) |
| 8x8 | 2 | 0 | 2597.6349 | 2939.4219 | 2939.0479 | oui |
| 8x8 | 2 | 1 | 2497.4444 | 2981.0770 | 2991.2716 | oui |
| 8x8 | 2 | 2 | 2474.8065 | 2988.9258 | 2976.2448 | oui |
| 8x8 | 2 | 3 | 2478.4724 | 3020.4216 | 2974.1813 | oui |
| 8x8 | 2 | 4 | — | 3084.4304 | 2991.2441 | non (« --- » dans le mémoire) |

`results/M/level1_summary.csv` n'a aucune ligne R_cut 4 (carte v1 à R_cut 0–3). Le point 9×9, R_cut 4 du mémoire (2468) est dans `results/M/m_rcut_resigma.csv` : med_Gamma_meV = 2468.003090813311.

### 5.3 tab:échantillonnage : ΔE_F et V̄^L_ed(a_CC), 8 tailles (meV)

Sources : `results/M2_plateau/sampling_table.csv` (colonnes `dE_F_meV`, `Ved_radial_1.42A_meV` ; fichier identique dans `results/M` et `results/M2` : md5 égaux) ; C_N : `config/production.json` alignment.C_N_eV (× 10³, 4 décimales). Colonne arithmétique : aligné = non aligné − C_N. ΔE_F ne dépend pas de l'alignement.

| N | N mod 3 | ΔE_F | V̄^L(a_CC) non aligné | C_N | V̄^L(a_CC) aligné = non aligné − C_N |
|---|---|---|---|---|---|
| 6 | 0 | -20.5 | +8.2 | -50.5091 | +58.7091 |
| 9 | 0 | -2.4 | +43.3 | -25.1437 | +68.4437 |
| 12 | 0 | +2.5 | +51.1 | -18.6896 | +69.7896 |
| 5 | 2 | -235.8 | -364.4 | -57.5626 | -306.8374 |
| 7 | 1 | -866.1 | -400.1 | -26.9143 | -373.1857 |
| 8 | 2 | -583.3 | -269.9 | -14.4103 | -255.4897 |
| 10 | 1 | -508.7 | -297.7 | -9.3417 | -288.3583 |
| 11 | 2 | -209.3 | -253.9 | -9.5507 | -244.3493 |

### 5.4 Scalaires lus sur les courbes des npz (définitions de `table_v1_v2.md` l.73, 74, 98, 99)

Lecture : Born/T = Γ_Born/Γ_T de `resonance_9x9.npz` (`eg`, `Gamma_Born`, `Gamma_T`) sur |ε − E_D| ≤ 3 eV, min et max ; Γ_T à c = 0,1 % = `Gamma_c` × 10³ de `resonance_criteria_9x9.npz` (`x_c` sur ±1 eV) : min (position), max (position), interpolation linéaire à 0 ; ħ/Γ = 658,2 meV·fs / Γ interpolé à −0,3 et +0,3 eV — code de `article/R6_production_corrigee/etape3/r6_compare_v1_v2.py` (rstats, cstats).

**Porte** (la même lecture sur les npz v1 et non alignés redonne les valeurs publiées à la dernière décimale) : **PASS**.

| ligne publiée | grandeur | npz | publié | relu | verdict |
|---|---|---|---|---|---|
| R6 l.73 | Born/T min sur ±3 eV | v1 | 0.677 | 0.677 | OK |
| R6 l.73 | Born/T min sur ±3 eV | non aligné | 1.106 | 1.106 | OK |
| R6 l.74 | Born/T max sur ±3 eV | v1 | 16.453 | 16.453 | OK |
| R6 l.74 | Born/T max sur ±3 eV | non aligné | 211.720 | 211.720 | OK |
| R6 l.98 | Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV) | v1 | 0.63 (+0.24) / 5.55 (-1.00) / 1.27 | 0.63 (+0.24) / 5.55 (-1.00) / 1.27 | OK |
| R6 l.98 | Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV) | non aligné | 2.25 (+1.00) / 37.84 (-0.18) / 7.79 | 2.25 (+1.00) / 37.84 (-0.18) / 7.79 | OK |
| R6 l.99 | ħ/Γ à ∓0,3 eV, c = 0,1 % (fs) | v1 | 419 / 1025 | 419 / 1025 | OK |
| R6 l.99 | ħ/Γ à ∓0,3 eV, c = 0,1 % (fs) | non aligné | 24 / 230 | 24 / 230 | OK |

Valeurs finales (même lecture sur les npz de `results/M2_plateau`) :

| grandeur | v1 | non aligné (v2) | final | final − v1 | statut | source finale | lieu |
|---|---|---|---|---|---|---|---|
| Born/T min sur ±3 eV | 0.677 | 1.106 | 1.050 | +0.373 | remplacé | `results/M2_plateau/resonance_9x9.npz` : eg, Gamma_Born, Gamma_T | R6 l.73 (sans équivalent final dans la table principale) |
| Born/T max sur ±3 eV | 16.453 | 211.720 | 202.851 | +186.398 | remplacé | `results/M2_plateau/resonance_9x9.npz` : eg, Gamma_Born, Gamma_T | R6 l.74 (sans équivalent final dans la table principale) |
| Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV) | 0.63 (+0.24) / 5.55 (-1.00) / 1.27 | 2.25 (+1.00) / 37.84 (-0.18) / 7.79 | 2.48 (+1.00) / 40.58 (-0.18) / 8.47 | +1.85 / +0.76 / +35.03 / +0.82 / +7.20 | remplacé | `results/M2_plateau/resonance_criteria_9x9.npz` : x_c, Gamma_c | R6 l.98 (sans équivalent final dans la table principale) |
| ħ/Γ à ∓0,3 eV, c = 0,1 % (fs) | 419 / 1025 | 24 / 230 | 26 / 212 | -393 / -813 | remplacé | `results/M2_plateau/resonance_criteria_9x9.npz` : x_c, Gamma_c | R6 l.99 (sans équivalent final dans la table principale) |
