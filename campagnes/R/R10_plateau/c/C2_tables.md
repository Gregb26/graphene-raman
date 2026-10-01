# R10 — C.2 : tab:rcut_M à trois colonnes et grandeur de Kaasbjerg (tables générées par r10_driver.py c2)

| porte | résultat |
|---|---|
| tab:rcut_M 9x9 : colonne Wigner-Seitz = m_rcut_convergence.csv de C.1 | OK |
| tab:rcut_M 12x12 : colonne Wigner-Seitz = m_rcut_convergence.csv de C.1 | OK |

## tab:rcut_M, 9x9 (max|ΔM|/max|M| de la paire π/π*, grille fine 60² ; C_N −25,1437 meV ; 26 étiquettes à égalité)

| R_cut | sites | tel quel | plateau (i), étiquettes actuelles | plateau (i), Wigner-Seitz | SV (WS) | diagonale (WS) |
|---|---|---|---|---|---|---|
| 0 | 1 | 3.1569e-01 | 3.1708e-01 | 3.1708e-01 | 3.1708e-01 | 3.1708e-01 |
| 1 | 5 | 2.1293e-01 | 1.6301e-01 | 1.6301e-01 | 1.8508e-01 | 1.5005e-01 |
| 2 | 13 | 1.1395e-01 | 5.5681e-02 | 5.5681e-02 | 6.6075e-02 | 4.6200e-02 |
| 3 | 29 | 6.6789e-02 | 2.4704e-02 | 2.5879e-02 | 2.5737e-02 | 1.7074e-02 |
| 4 | 49 | 4.4189e-02 | 1.7179e-02 | 1.8194e-02 | 1.9846e-02 | 1.4106e-02 |
| 5 | 81 | 1.7488e-02 | 1.3578e-02 | 1.4950e-02 | 1.1040e-02 | 1.1945e-02 |
| 6 | 113 | 1.2068e-02 | 1.1992e-02 | 1.3619e-02 | 9.1728e-03 | 1.0052e-02 |

max|M_π| (eV) : tel quel 25.4397, plateau 25.6014, plateau WS 25.6014 ; écart relatif max à R9 : tel quel 0.0e+00, plateau (étiquettes actuelles) 0.0e+00 ; colonne WS contre le csv de C.1 : identique.

## tab:rcut_M, 12x12 (max|ΔM|/max|M| de la paire π/π*, grille fine 60² ; C_N −18,6896 meV ; 23 étiquettes à égalité)

| R_cut | sites | tel quel | plateau (i), étiquettes actuelles | plateau (i), Wigner-Seitz | SV (WS) | diagonale (WS) |
|---|---|---|---|---|---|---|
| 0 | 1 | 8.6729e-01 | 8.6781e-01 | 8.6782e-01 | 9.9271e-01 | 8.6782e-01 |
| 1 | 5 | 2.1719e-01 | 1.8923e-01 | 1.8925e-01 | 1.8948e-01 | 1.8925e-01 |
| 2 | 13 | 1.2490e-01 | 1.2505e-01 | 1.2506e-01 | 1.2521e-01 | 1.2506e-01 |
| 3 | 29 | 9.0206e-02 | 3.0966e-02 | 3.0979e-02 | 3.2908e-02 | 3.0979e-02 |
| 4 | 49 | 7.5696e-02 | 1.8398e-02 | 1.8573e-02 | 2.5296e-02 | 1.5291e-02 |
| 5 | 81 | 5.3559e-02 | 1.2858e-02 | 1.2739e-02 | 1.6202e-02 | 8.5068e-03 |
| 6 | 113 | 3.2881e-02 | 1.2725e-02 | 1.1781e-02 | 1.1178e-02 | 8.6769e-03 |

max|M_π| (eV) : tel quel 25.9402, plateau 26.1611, plateau WS 26.1611 ; écart relatif max à R9 : tel quel 0.0e+00, plateau (étiquettes actuelles) 0.0e+00 ; colonne WS contre le csv de C.1 : identique.

## D.1 : Ṽ = A_cell·|M_k'k|, k = K + δx̂ (δ = 0.0294 Å⁻¹), carte 240², disques de rayon 0.1471 Å⁻¹ (2 points exclus) ; eV Å²

| variante | valence K | conduction K | valence K′ | conduction K′ |
|---|---|---|---|---|
| tel quel + Wigner-Seitz | 76.956 [75.51 ; 78.73] | 78.168 [76.16 ; 80.36] | 81.177 [80.48 ; 82.33] | 82.305 [81.52 ; 83.71] |
| plateau + Wigner-Seitz | 81.605 [77.40 ; 86.36] | 82.821 [77.81 ; 87.47] | 81.168 [80.53 ; 82.31] | 82.297 [81.52 ; 83.69] |
| R9 brut (étiquettes brutes) | 76.955 | 78.159 | 81.176 | 82.297 |
| R9 exact (étiquettes brutes) | 81.512 | 82.706 | 81.168 | 82.289 |

## D.2 : bloc 2 × 2 de la paire π/π* à (K, K) (eV Å²)

| partie | variante | valeurs propres | ½ Tr | normes des lignes | ‖·‖_F |
|---|---|---|---|---|---|
| tot | tel quel | -1.71 ; 153.45 | 75.87 | 107.22 ; 109.79 | 153.46 |
| tot | plateau | 9.02 ; 164.17 | 86.59 | 114.89 ; 117.62 | 164.42 |
| L | tel quel | -1.71 ; 120.77 | 59.53 | 84.39 ; 86.41 | 120.78 |
| L | plateau | 9.02 ; 131.49 | 70.25 | 92.10 ; 94.29 | 131.80 |
| NL | tel quel | -0.00 ; 32.68 | 16.34 | 22.83 ; 23.38 | 32.68 |

Kaasbjerg (PRB 101, 045433, Fig. 3) : ~70 eV Å² près de K et K′ ; V₀ ≈ 27 eV ; super-cellule 11×11.
