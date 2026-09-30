# R6 3.3 — porte de niveau 1 sur les paramètres gelés (seuil 5 %, critère du 2026-09-05)

Config : R_cut 3, grille 240, η 0.02 eV, N_k^int 300, référence 9x9. Médianes de |Γ|·N_cells (meV) sur la fenêtre ±3.0 eV.

| taille | médiane (R_cut 3, 240², η 0,02) | E_res (eV) | C11 grille 120 → 240 | C11 η 0,02 → 0,01 | C10 R_cut 3 → 4 | verdict |
|---|---|---|---|---|---|---|
| 5x5 | 3288.90 | -0.284 | 0.34 % | 1.06 % | 0.42 % | OK |
| 6x6 | 3155.98 | -0.227 | 0.01 % | 0.52 % | 1.58 % | OK |
| 7x7 | 3052.03 | -0.269 | 0.09 % | 0.53 % | 2.99 % | OK |
| 8x8 | 3020.42 | -0.269 | 0.10 % | 0.45 % | 2.08 % | OK |
| 9x9 | 3132.60 | -0.175 | 0.49 % | 0.10 % | 0.47 % | OK |
| 12x12 | 3162.78 | -0.175 | 0.41 % | 0.44 % | 0.04 % | OK |

N_k^int (9x9, R_cut 3, grille 240, η 0.02) : médiane de |Γ|·N_cells sur les états de la fenêtre ; référence N_k^int = 600
| N_k^int | n états | médiane (meV) | écart vs max | écart vs 300 | E_res (eV) | Γ_T(E_D) (meV, états à K) | écart Γ(E_D) vs max |
|---|---|---|---|---|---|---|---|
| 150 | 41266 | 3155.57 | +0.64 % | +0.73 % | -0.238 | 4801.86 | +38.63 % |
| 300 | 41266 | 3132.60 | +0.10 % | +0.00 % | -0.175 | 3630.65 | +4.81 % |
| 450 | 41266 | 3135.87 | +0.01 % | +0.10 % | -0.202 | 3489.52 | +0.74 % |
| 600 | 41266 | 3135.65 | +0.00 % | +0.10 % | -0.202 | 3463.87 | +0.00 % |

R_cut par `rcut_resigma.py` (9x9, grille 240, η 0.02, N_k^int 300) :
| R_cut | médiane (meV) | écart vs R_cut max | E_res (eV) |
|---|---|---|---|
| 0 | 3009.73 | +4.38 % | -0.175 |
| 1 | 3223.74 | +2.42 % | -0.227 |
| 2 | 3170.21 | +0.72 % | -0.175 |
| 3 | 3132.60 | +0.47 % | -0.175 |
| 4 | 3147.47 | +0.00 % | -0.175 |

**Verdict** : OK — paramètres gelés conformes au critère du 2026-09-05 avec M2 ; 3.4 peut démarrer
