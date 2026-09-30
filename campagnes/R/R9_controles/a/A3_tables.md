# R9 — A.3 : variantes alignées (tables générées par r9_driver.py a3)

C14 de R6 (`resonance_9x9_shiftL.npz`, lecture) : C = +25 meV : écart relatif maximal de Γ_T 1.920e-01 à -0.115 eV (médian 3.40e-02) ; C = -25 meV : écart relatif maximal de Γ_T 1.633e-01 à -0.120 eV (médian 4.02e-02)

## Convention intensive (max|M|, bandes 1–16, eV)

| ensemble | taille | brut | aligné | mode aligné | C_N (meV) | max en k = k' (brut ; aligné) |
|---|---|---|---|---|---|---|
| grossiers 5, 7, 8, 9 | 5x5 | 25.198 | 26.625 | exact (N_cells C_N 𝕀) | -122.33 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| grossiers 5, 7, 8, 9 | 7x7 | 23.250 | 24.498 | exact (N_cells C_N 𝕀) | -38.54 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| grossiers 5, 7, 8, 9 | 8x8 | 24.311 | 24.444 | exact (N_cells C_N 𝕀) | -22.97 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| grossiers 5, 7, 8, 9 | 9x9 | 24.801 | 25.561 | exact (N_cells C_N 𝕀) | -24.65 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 5x5 | 23.850 | 25.717 | exact (M2 − C_N M^L[1_boîte]) | -122.33 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 6x6 | 25.148 | 27.454 | blocs k = k' seulement | -98.62 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 7x7 | 23.818 | 24.498 | blocs k = k' seulement | -38.54 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 8x8 | 24.311 | 24.444 | blocs k = k' seulement | -22.97 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 9x9 | 25.348 | 25.677 | exact (M2 − C_N M^L[1_boîte]) | -24.65 (sans plateau : valeur Lu (1,0 Å)) | False ; True |
| denses six tailles | 12x12 | 25.721 | 27.939 | blocs k = k' seulement | -29.24 (sans plateau : valeur Lu (1,0 Å)) | False ; True |

| ensemble | variante | (max − min)/moyenne |
|---|---|---|
| grossiers 5, 7, 8, 9 | brut | 7.983e-02 |
| grossiers 5, 7, 8, 9 | aligne | 8.627e-02 |
| denses six tailles | brut | 7.706e-02 |
| denses six tailles | aligne | 1.347e-01 |
| denses 5, 7, 8, 9 | brut | 6.291e-02 |
| denses 5, 7, 8, 9 | aligne | 5.073e-02 |

## tab:rcut_M, 9x9 (max|ΔM|/max|M| de la paire π, grille fine 60² ; C_N -24.65 meV, sans plateau : valeur Lu (1,0 Å))

| R_cut | sites | M2 tel quel | aligné (i) | aligné exact (F_W) |
|---|---|---|---|---|
| 0 | 1 | 3.157e-01 | 3.160e-01 | 3.160e-01 |
| 1 | 5 | 2.129e-01 | 1.632e-01 | 1.624e-01 |
| 2 | 13 | 1.139e-01 | 5.583e-02 | 5.576e-02 |
| 3 | 29 | 6.679e-02 | 2.476e-02 | 2.496e-02 |
| 4 | 49 | 4.419e-02 | 1.713e-02 | 1.722e-02 |
| 5 | 81 | 1.749e-02 | 1.367e-02 | 1.398e-02 |
| 6 | 113 | 1.207e-02 | 1.201e-02 | 1.210e-02 |

## tab:rcut_M, 12x12 (max|ΔM|/max|M| de la paire π, grille fine 60² ; C_N -29.24 meV, sans plateau : valeur Lu (1,0 Å))

| R_cut | sites | M2 tel quel | aligné (i) |
|---|---|---|---|
| 0 | 1 | 8.673e-01 | 8.747e-01 |
| 1 | 5 | 2.172e-01 | 2.318e-01 |
| 2 | 13 | 1.249e-01 | 1.681e-01 |
| 3 | 29 | 9.021e-02 | 7.308e-02 |
| 4 | 49 | 7.570e-02 | 5.198e-02 |
| 5 | 81 | 5.356e-02 | 3.260e-02 |
| 6 | 113 | 3.288e-02 | 2.168e-02 |

## Niveau 1 (240², η 0,02, N_k^int 300) et courbe Γ_T (grille de resonance_metrics)

| taille | R_cut | variante | médiane Γ·N_cells (meV) | E_res | Re Σ médian (meV) | pic Γ_T (5 meV ; 2,5 meV) | pic −Im T̄(K) | Re T̄(E_D) | sites hors boîte | C_N (meV) | régression R6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5x5 | 3 | M2 tel quel | 3288.90 | -0.284 | -1429.90 | -0.305 ; -0.3025 | -0.2700 | +6.519 | 4 | +0.00 (—) | OK |
| 5x5 | 3 | aligné (i) | 3140.34 | -0.175 | 1174.43 | -0.180 ; -0.1800 | -0.1775 | +12.025 | 4 | -122.33 (sans plateau : valeur Lu (1,0 Å)) | — |
| 5x5 | 3 | aligné exact (F_W) | 3130.22 | -0.175 | 1133.42 | -0.180 ; -0.1800 | -0.1775 | +11.948 | 4 | -122.33 (sans plateau : valeur Lu (1,0 Å)) | — |
| 5x5 | 4 | M2 tel quel | 3302.74 | -0.284 | -1437.73 | -0.305 ; -0.3025 | -0.2700 | +6.452 | 24 | +0.00 (—) | OK |
| 5x5 | 4 | aligné (i) | 3143.55 | -0.181 | 1169.89 | -0.180 ; -0.1800 | -0.1775 | +11.921 | 24 | -122.33 (sans plateau : valeur Lu (1,0 Å)) | — |
| 5x5 | 4 | aligné exact (F_W) | 3143.25 | -0.181 | 1172.10 | -0.180 ; -0.1825 | -0.1775 | +11.889 | 24 | -122.33 (sans plateau : valeur Lu (1,0 Å)) | — |
| 6x6 | 3 | M2 tel quel | 3155.98 | -0.227 | -695.92 | -0.240 ; -0.2375 | -0.2200 | +9.218 | 2 | +0.00 (—) | OK |
| 6x6 | 3 | aligné (i) | 3439.97 | -0.134 | 1560.98 | -0.175 ; -0.1775 | -0.1675 | +14.805 | 2 | -98.62 (sans plateau : valeur Lu (1,0 Å)) | — |
| 6x6 | 4 | M2 tel quel | 3206.64 | -0.227 | -996.33 | -0.240 ; -0.2400 | -0.2225 | +8.764 | 14 | +0.00 (—) | OK |
| 6x6 | 4 | aligné (i) | 3556.00 | -0.134 | 1999.29 | -0.175 ; -0.1775 | -0.1675 | +15.260 | 14 | -98.62 (sans plateau : valeur Lu (1,0 Å)) | — |
| 7x7 | 3 | M2 tel quel | 3052.03 | -0.269 | -1033.37 | -0.300 ; -0.3000 | -0.2675 | +7.165 | 0 | +0.00 (—) | OK |
| 7x7 | 3 | aligné (i) | 2927.30 | -0.227 | -62.97 | -0.240 ; -0.2425 | -0.2625 | +8.961 | 0 | -38.54 (sans plateau : valeur Lu (1,0 Å)) | — |
| 7x7 | 4 | M2 tel quel | 3146.18 | -0.269 | -1398.98 | -0.300 ; -0.3025 | -0.2675 | +6.696 | 4 | +0.00 (—) | OK |
| 7x7 | 4 | aligné (i) | 2944.65 | -0.227 | 141.59 | -0.240 ; -0.2425 | -0.2625 | +9.185 | 4 | -38.54 (sans plateau : valeur Lu (1,0 Å)) | — |
| 8x8 | 3 | M2 tel quel | 3020.42 | -0.269 | -796.51 | -0.245 ; -0.2450 | -0.2650 | +7.930 | 0 | +0.00 (—) | OK |
| 8x8 | 3 | aligné (i) | 2958.71 | -0.227 | -216.25 | -0.240 ; -0.2425 | -0.2250 | +9.036 | 0 | -22.97 (sans plateau : valeur Lu (1,0 Å)) | — |
| 8x8 | 4 | M2 tel quel | 3084.43 | -0.269 | -1141.13 | -0.275 ; -0.2775 | -0.2650 | +7.486 | 2 | +0.00 (—) | OK |
| 8x8 | 4 | aligné (i) | 2970.15 | -0.227 | -172.31 | -0.240 ; -0.2425 | -0.2250 | +9.057 | 2 | -22.97 (sans plateau : valeur Lu (1,0 Å)) | — |
| 9x9 | 3 | M2 tel quel | 3132.60 | -0.175 | -26.65 | -0.180 ; -0.1800 | -0.1775 | +11.122 | 0 | +0.00 (—) | OK |
| 9x9 | 3 | aligné (i) | 3185.03 | -0.175 | 592.62 | -0.180 ; -0.1775 | -0.1725 | +12.594 | 0 | -24.65 (sans plateau : valeur Lu (1,0 Å)) | — |
| 9x9 | 3 | aligné exact (F_W) | 3185.02 | -0.175 | 592.55 | -0.180 ; -0.1775 | -0.1725 | +12.594 | 0 | -24.65 (sans plateau : valeur Lu (1,0 Å)) | — |
| 9x9 | 4 | M2 tel quel | 3147.47 | -0.175 | -410.29 | -0.180 ; -0.1800 | -0.1775 | +10.556 | 0 | +0.00 (—) | OK |
| 9x9 | 4 | aligné (i) | 3196.41 | -0.175 | 674.74 | -0.180 ; -0.1775 | -0.1750 | +12.606 | 0 | -24.65 (sans plateau : valeur Lu (1,0 Å)) | — |
| 9x9 | 4 | aligné exact (F_W) | 3196.02 | -0.175 | 671.57 | -0.180 ; -0.1775 | -0.1750 | +12.603 | 0 | -24.65 (sans plateau : valeur Lu (1,0 Å)) | — |
| 12x12 | 3 | M2 tel quel | 3162.78 | -0.175 | 448.29 | -0.180 ; -0.1775 | -0.1750 | +12.429 | 0 | +0.00 (—) | OK |
| 12x12 | 3 | aligné (i) | 3332.74 | -0.134 | 1166.10 | -0.175 ; -0.1775 | -0.1675 | +14.323 | 0 | -29.24 (sans plateau : valeur Lu (1,0 Å)) | — |
| 12x12 | 4 | M2 tel quel | 3161.49 | -0.175 | 249.19 | -0.180 ; -0.1800 | -0.1750 | +12.104 | 0 | +0.00 (—) | OK |
| 12x12 | 4 | aligné (i) | 3395.10 | -0.134 | 1514.67 | -0.175 ; -0.1775 | -0.1675 | +14.722 | 0 | -29.24 (sans plateau : valeur Lu (1,0 Å)) | — |

Familles (médianes Γ·N_cells) : R_cut 3 | brut | 3m : moyenne 3150.5 meV, (max − min)/moyenne 0.96 % ; R_cut 3 | brut | non-3m : moyenne 3120.4 meV, (max − min)/moyenne 8.60 % ; R_cut 3 | aligne | 3m : moyenne 3319.2 meV, (max − min)/moyenne 7.68 % ; R_cut 3 | aligne | non-3m : moyenne 3008.8 meV, (max − min)/moyenne 7.08 % ; R_cut 4 | brut | 3m : moyenne 3171.9 meV, (max − min)/moyenne 1.87 % ; R_cut 4 | brut | non-3m : moyenne 3177.8 meV, (max − min)/moyenne 6.87 % ; R_cut 4 | aligne | 3m : moyenne 3382.5 meV, (max − min)/moyenne 10.63 % ; R_cut 4 | aligne | non-3m : moyenne 3019.4 meV, (max − min)/moyenne 6.59 %

## Critère de pôle, 9×9, 300², fenêtre ±3 eV

| variante | bloc | min\|det\|/max (ε − E_D) | min\|λ\| (ε − E_D) ; λ |
|---|---|---|---|
| M2 tel quel | pi | 2.330e-02 (-0.172) | 0.3967 (-0.170) ; +0.2186+0.3310i |
| M2 tel quel | complet | 1.318e-04 (-0.812) | 0.0019 (-0.812) ; +0.0000+0.0019i |
| aligné (i) | pi | 1.717e-02 (-0.170) | 0.3450 (-0.127) ; +0.2507+0.2370i |
| aligné (i) | complet | 1.260e-04 (-0.787) | 0.0019 (-0.787) ; +0.0000+0.0019i |
| aligné exact (F_W) | pi | 1.717e-02 (-0.170) | 0.3450 (-0.127) ; +0.2507+0.2370i |
| aligné exact (F_W) | complet | 1.260e-04 (-0.787) | 0.0019 (-0.787) ; +0.0000+0.0019i |
R6 (M2 tel quel) : π 2,33e-2 (−0,172), |λ| 0,397 (−0,170) ; complet 1,32e-4 (−0,812), |λ| 0,0019 (−0,812).

