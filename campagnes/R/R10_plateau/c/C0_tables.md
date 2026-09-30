# R10 — C.0 : portes (tables générées par r10_driver.py c0)

| porte | résultat |
|---|---|
| a niveau 1 9x9 | OK |
| b niveau 1 9x9 | OK |
| b V_loc au bit 9x9 | OK |
| a tab:rcut_M 9x9 | OK |
| b tab:rcut_M 9x9 | OK |
| c M sur la grille MP | OK |
| c D.1 sans ws = R9 D.1 | OK |
| c D.1 ws = audit « vraie » | OK |
| b niveau 1 12x12 | OK |
| b V_loc au bit 12x12 | OK |
| b tab:rcut_M 12x12 | OK |

| taille | variante | médiane Γ·N_cells (meV) ; R9 | E_res (eV) ; R9 | tab:rcut_M R_cut 3 ; R9 (écart relatif) |
|---|---|---|---|---|
| 9x9 | C_N = 0 | 3132,5968 ; 3132,5968 | −0,174849 ; −0,174849 | 6.678859e-02 ; 6.678859e-02 (0.0e+00) |
| 9x9 | plateau (C_N −25,1437 meV) | 3189,0136 ; 3189,0136 | −0,174849 ; −0,174849 | 2.470441e-02 ; 2.470441e-02 (0.0e+00) |
| 12x12 | plateau (C_N −18,6896 meV) | 3264,5060 ; 3264,5060 | −0,174929 ; −0,174929 | 3.096624e-02 ; 3.096624e-02 (0.0e+00) |

(c) Wigner-Seitz sur la grille MP 27² : écart relatif 4.1e-16 (étiquettes brutes), 3.9e-16 (recentrées) ; 26 étiquettes à égalité, 757 images.

| D.1 (eV Å²) | valence K | conduction K | valence K′ | conduction K′ |
|---|---|---|---|---|
| tel quel, sans ws | 76.955098 | 78.159478 | 81.176183 | 82.296735 |
| tel quel, Wigner-Seitz | 76.955553 | 78.168368 | 81.176724 | 82.305361 |
| (i) C = Lu (−24,6514 meV), sans ws | 81.512905 | 82.720185 | 81.167595 | 82.288146 |
| R9 D.1 exact (F_W) | 81.512262 | 82.706262 | 81.167865 | 82.288636 |
| écart (i) − exact | +0.000643 | +0.013923 | -0.000270 | -0.000490 |

Écarts max des moyennes : sans ws − R9 D.1 5.7e-14 ; ws − audit « vraie | brut » 1.4e-14 eV Å².
