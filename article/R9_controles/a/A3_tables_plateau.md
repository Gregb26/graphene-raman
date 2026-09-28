# R9 — A.3 : variantes alignées (tables générées par r9_driver.py a3)

C14 de R6 (`resonance_9x9_shiftL.npz`, lecture) : C = +25 meV : écart relatif maximal de Γ_T 1.920e-01 à -0.115 eV (médian 3.40e-02) ; C = -25 meV : écart relatif maximal de Γ_T 1.633e-01 à -0.120 eV (médian 4.02e-02)

## Convention intensive (max|M|, bandes 1–16, eV)

| ensemble | taille | brut | aligné | mode aligné | C_N (meV) | max en k = k' (brut ; aligné) |
|---|---|---|---|---|---|---|
| grossiers 5, 7, 8, 9 | 5x5 | 25.198 | 25.198 | exact (N_cells C_N 𝕀) | -57.56 (plateau (i)) | False ; False |
| grossiers 5, 7, 8, 9 | 7x7 | 23.250 | 23.928 | exact (N_cells C_N 𝕀) | -26.91 (plateau (i)) | False ; True |
| grossiers 5, 7, 8, 9 | 8x8 | 24.311 | 24.311 | exact (N_cells C_N 𝕀) | -14.41 (plateau (i)) | False ; False |
| grossiers 5, 7, 8, 9 | 9x9 | 24.801 | 25.601 | exact (N_cells C_N 𝕀) | -25.14 (plateau (i)) | False ; True |
| denses six tailles | 5x5 | 23.850 | 24.097 | exact (M2 − C_N M^L[1_boîte]) | -57.56 (plateau (i)) | False ; True |
| denses six tailles | 6x6 | 25.148 | 25.722 | blocs k = k' seulement | -50.51 (plateau (i)) | False ; True |
| denses six tailles | 7x7 | 23.818 | 23.928 | blocs k = k' seulement | -26.91 (plateau (i)) | False ; True |
| denses six tailles | 8x8 | 24.311 | 24.311 | blocs k = k' seulement | -14.41 (plateau (i)) | False ; False |
| denses six tailles | 9x9 | 25.348 | 25.717 | exact (M2 − C_N M^L[1_boîte]) | -25.14 (plateau (i)) | False ; True |
| denses six tailles | 12x12 | 25.721 | 26.420 | blocs k = k' seulement | -18.69 (plateau (i)) | False ; True |

| ensemble | variante | (max − min)/moyenne |
|---|---|---|
| grossiers 5, 7, 8, 9 | brut | 7.983e-02 |
| grossiers 5, 7, 8, 9 | aligne_plateau | 6.757e-02 |
| denses six tailles | brut | 7.706e-02 |
| denses six tailles | aligne_plateau | 9.956e-02 |
| denses 5, 7, 8, 9 | brut | 6.291e-02 |
| denses 5, 7, 8, 9 | aligne_plateau | 7.299e-02 |

## tab:rcut_M, 9x9 (max|ΔM|/max|M| de la paire π, grille fine 60² ; C_N -25.14 meV, plateau (i))

| R_cut | sites | M2 tel quel | aligné (i), C_N plateau | aligné exact (F_W), C_N plateau |
|---|---|---|---|---|
| 0 | 1 | 3.157e-01 | 3.171e-01 | 3.171e-01 |
| 1 | 5 | 2.129e-01 | 1.630e-01 | 1.621e-01 |
| 2 | 13 | 1.139e-01 | 5.568e-02 | 5.562e-02 |
| 3 | 29 | 6.679e-02 | 2.470e-02 | 2.495e-02 |
| 4 | 49 | 4.419e-02 | 1.718e-02 | 1.726e-02 |
| 5 | 81 | 1.749e-02 | 1.358e-02 | 1.389e-02 |
| 6 | 113 | 1.207e-02 | 1.199e-02 | 1.209e-02 |

## tab:rcut_M, 12x12 (max|ΔM|/max|M| de la paire π, grille fine 60² ; C_N -18.69 meV, plateau (i))

| R_cut | sites | M2 tel quel | aligné (i), C_N plateau |
|---|---|---|---|
| 0 | 1 | 8.673e-01 | 8.678e-01 |
| 1 | 5 | 2.172e-01 | 1.892e-01 |
| 2 | 13 | 1.249e-01 | 1.250e-01 |
| 3 | 29 | 9.021e-02 | 3.097e-02 |
| 4 | 49 | 7.570e-02 | 1.840e-02 |
| 5 | 81 | 5.356e-02 | 1.286e-02 |
| 6 | 113 | 3.288e-02 | 1.272e-02 |

## Niveau 1 (240², η 0,02, N_k^int 300) et courbe Γ_T (grille de resonance_metrics)

| taille | R_cut | variante | médiane Γ·N_cells (meV) | E_res | Re Σ médian (meV) | pic Γ_T (5 meV ; 2,5 meV) | pic −Im T̄(K) | Re T̄(E_D) | sites hors boîte | C_N (meV) | régression R6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5x5 | 3 | aligné (i), C_N plateau | 2997.14 | -0.227 | -226.16 | -0.240 ; -0.2425 | -0.2250 | +8.940 | 4 | -57.56 (plateau (i)) | — |
| 5x5 | 3 | aligné exact (F_W), C_N plateau | 2999.20 | -0.227 | -246.08 | -0.240 ; -0.2425 | -0.2250 | +8.907 | 4 | -57.56 (plateau (i)) | — |
| 5x5 | 4 | aligné (i), C_N plateau | 2997.07 | -0.227 | -233.34 | -0.240 ; -0.2425 | -0.2625 | +8.858 | 24 | -57.56 (plateau (i)) | — |
| 5x5 | 4 | aligné exact (F_W), C_N plateau | 2999.31 | -0.227 | -232.32 | -0.240 ; -0.2425 | -0.2625 | +8.846 | 24 | -57.56 (plateau (i)) | — |
| 6x6 | 3 | aligné (i), C_N plateau | 3149.15 | -0.175 | 474.13 | -0.180 ; -0.1800 | -0.1775 | +11.884 | 2 | -50.51 (plateau (i)) | — |
| 6x6 | 4 | aligné (i), C_N plateau | 3161.15 | -0.175 | 553.63 | -0.180 ; -0.1800 | -0.1775 | +11.872 | 14 | -50.51 (plateau (i)) | — |
| 7x7 | 3 | aligné (i), C_N plateau | 2942.85 | -0.269 | -356.03 | -0.245 ; -0.2450 | -0.2625 | +8.407 | 0 | -26.91 (plateau (i)) | — |
| 7x7 | 4 | aligné (i), C_N plateau | 2953.28 | -0.269 | -328.65 | -0.245 ; -0.2450 | -0.2625 | +8.418 | 4 | -26.91 (plateau (i)) | — |
| 8x8 | 3 | aligné (i), C_N plateau | 2974.18 | -0.227 | -435.27 | -0.245 ; -0.2425 | -0.2625 | +8.618 | 0 | -14.41 (plateau (i)) | — |
| 8x8 | 4 | aligné (i), C_N plateau | 2991.24 | -0.269 | -537.40 | -0.245 ; -0.2425 | -0.2625 | +8.464 | 2 | -14.41 (plateau (i)) | — |
| 9x9 | 3 | aligné (i), C_N plateau | 3189.01 | -0.175 | 605.08 | -0.180 ; -0.1775 | -0.1725 | +12.624 | 0 | -25.14 (plateau (i)) | — |
| 9x9 | 3 | aligné exact (F_W), C_N plateau | 3189.00 | -0.175 | 605.01 | -0.180 ; -0.1775 | -0.1725 | +12.624 | 0 | -25.14 (plateau (i)) | — |
| 9x9 | 4 | aligné (i), C_N plateau | 3200.00 | -0.175 | 696.74 | -0.180 ; -0.1775 | -0.1750 | +12.649 | 0 | -25.14 (plateau (i)) | — |
| 9x9 | 4 | aligné exact (F_W), C_N plateau | 3199.55 | -0.175 | 693.63 | -0.180 ; -0.1775 | -0.1750 | +12.645 | 0 | -25.14 (plateau (i)) | — |
| 12x12 | 3 | aligné (i), C_N plateau | 3264.51 | -0.175 | 913.84 | -0.175 ; -0.1775 | -0.1700 | +13.619 | 0 | -18.69 (plateau (i)) | — |
| 12x12 | 4 | aligné (i), C_N plateau | 3283.17 | -0.175 | 1070.69 | -0.175 ; -0.1775 | -0.1700 | +13.749 | 0 | -18.69 (plateau (i)) | — |

Familles (médianes Γ·N_cells) : R_cut 3 | aligne_plateau | 3m : moyenne 3200.9 meV, (max − min)/moyenne 3.60 % ; R_cut 3 | aligne_plateau | non-3m : moyenne 2971.4 meV, (max − min)/moyenne 1.83 % ; R_cut 4 | aligne_plateau | 3m : moyenne 3214.8 meV, (max − min)/moyenne 3.80 % ; R_cut 4 | aligne_plateau | non-3m : moyenne 2980.5 meV, (max − min)/moyenne 1.47 %

## Critère de pôle, 9×9, 300², fenêtre ±3 eV

| variante | bloc | min\|det\|/max (ε − E_D) | min\|λ\| (ε − E_D) ; λ |
|---|---|---|---|
| aligné (i), C_N plateau | pi | 1.706e-02 (-0.170) | 0.3440 (-0.127) ; +0.2493+0.2370i |
| aligné (i), C_N plateau | complet | 1.261e-04 (-0.785) | 0.0019 (-0.787) ; -0.0000+0.0019i |
| aligné exact (F_W), C_N plateau | pi | 1.706e-02 (-0.170) | 0.3440 (-0.127) ; +0.2493+0.2370i |
| aligné exact (F_W), C_N plateau | complet | 1.261e-04 (-0.785) | 0.0019 (-0.787) ; -0.0000+0.0019i |
R6 (M2 tel quel) : π 2,33e-2 (−0,172), |λ| 0,397 (−0,170) ; complet 1,32e-4 (−0,812), |λ| 0,0019 (−0,812).

