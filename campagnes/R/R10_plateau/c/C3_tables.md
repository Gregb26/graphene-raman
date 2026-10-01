# R10 — C.3 : contrôles des sorties touchées par l'audit (tables générées par r10_driver.py c3)

| porte | résultat |
|---|---|
| P-c2 5x5 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 5x5 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 5x5 : appartenance aux anneaux identique | OK |
| P-c2 6x6 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 6x6 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | REFUSÉE |
| P-c2 6x6 : appartenance aux anneaux identique | OK |
| P-c2 7x7 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 7x7 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 7x7 : appartenance aux anneaux identique | OK |
| P-c2 8x8 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 8x8 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 8x8 : appartenance aux anneaux identique | OK |
| P-c2 9x9 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 9x9 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 9x9 : appartenance aux anneaux identique | OK |
| P-c2 10x10 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 10x10 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 10x10 : appartenance aux anneaux identique | OK |
| P-c2 11x11 rad : anneaux sous 0,433 a_sc à 1e-13 eV | REFUSÉE |
| P-c2 11x11 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 11x11 : appartenance aux anneaux identique | OK |
| P-c2 12x12 rad : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 12x12 rad_masked : anneaux sous 0,433 a_sc à 1e-13 eV | OK |
| P-c2 12x12 : appartenance aux anneaux identique | OK |
| P-b2 5x5_dense | OK |
| P-b2 5x5_coarse | OK |
| P-b2 6x6_dense | OK |
| P-b2 7x7_dense | OK |
| P-b2 7x7_coarse | OK |
| P-b2 8x8_dense | OK |
| P-b2 8x8_coarse | OK |
| P-b2 9x9_dense | OK |
| P-b2 12x12_dense | OK |
| tab:tests_M famille 3m | OK |
| tab:tests_M famille non-3m | OK |
| tab:tests_M : ancienne ligne six tailles retirée | OK |

## P-c2 : anneaux d'analyze_Ved (pas 0,05 Å)

| N | 0,433 a_sc (Å) | anneaux entiers dessous | points ; changent d'anneau | écart max profil (eV) à r | écart max masqué (eV) à r | anneaux changés au-delà (audit) | max au-delà (meV) (audit) | écart à l'audit « vraie » (meV) |
|---|---|---|---|---|---|---|---|---|
| 5x5 | 5,339 | 106 | 15079 ; 0 | 0.00e+00 à 0,025 | 0.00e+00 à 0,025 | 18 (18) | 16,46 (16,46) | 0.0e+00 |
| 6x6 | 6,406 | 128 | 21955 ; 0 | 0.00e+00 à 0,025 | 1.21e-13 à 1,775 | 20 (20) | 14,24 (14,24) | 0.0e+00 |
| 7x7 | 7,474 | 149 | 31523 ; 0 | 0.00e+00 à 0,025 | 4.97e-14 à 0,625 | 24 (24) | 6,73 (6,73) | 0.0e+00 |
| 8x8 | 8,542 | 170 | 38797 ; 0 | 4.26e-14 à 0,425 | 0.00e+00 à 0,025 | 28 (28) | 1,84 (1,84) | 0.0e+00 |
| 9x9 | 9,610 | 192 | 49459 ; 0 | 1.17e-14 à 9,525 | 4.26e-14 à 0,175 | 30 (30) | 1,39 (1,39) | 0.0e+00 |
| 10x10 | 10,677 | 213 | 60907 ; 0 | 2.08e-14 à 7,425 | 0.00e+00 à 0,025 | 34 (34) | 0,93 (0,93) | 0.0e+00 |
| 11x11 | 11,745 | 234 | 87464 ; 0 | 1.01e-13 à 9,025 | 7.11e-14 à 0,375 | 38 (37) | 1,15 (1,15) | 0.0e+00 |
| 12x12 | 12,813 | 256 | 87949 ; 0 | 6.06e-14 à 2,875 | 2.17e-14 à 8,875 | 40 (40) | 2,17 (2,17) | 0.0e+00 |

9×9, effectifs en vraie image : anneau 9.625 Å : 534 points ; anneau 11.075 Å : 518 points (audit : 534, 518).

## P-b2 : abscisses de mwr_locality (unités de a)

| taille | D | abscisse max (audit) | points déplacés (audit) | poids hors site − M2 | p_z–p_z sur site − M2 (meV) ; −C_N |
|---|---|---|---|---|---|
| 5x5_dense | 25 | 13.8924 (13.8924) | 96 (96) | 0.0e+00 | +57,5626 ; +57,5626 |
| 5x5_coarse | 5 | 2.6458 (2.6458) | 2 (2) | 0.0e+00 | +57,5626 ; +57,5626 |
| 6x6_dense | 24 | 13.8564 (13.8564) | 89 (89) | 0.0e+00 | +50,5091 ; +50,5091 |
| 7x7_dense | 28 | 15.6205 (15.6205) | 123 (123) | 0.0e+00 | +26,9143 ; +26,9143 |
| 7x7_coarse | 7 | 3.6056 (3.6056) | 6 (6) | 0.0e+00 | +26,9143 ; +26,9143 |
| 8x8_dense | 32 | 18.1934 (18.1934) | 161 (161) | 0.0e+00 | +14,4103 ; +14,4103 |
| 8x8_coarse | 8 | 4.3589 (4.3589) | 9 (9) | 0.0e+00 | +14,4103 ; +14,4103 |
| 9x9_dense | 27 | 15.5885 (15.5885) | 112 (112) | 0.0e+00 | +25,1437 ; +25,1437 |
| 12x12_dense | 24 | 13.8564 (13.8564) | 89 (89) | 0.0e+00 | +18,6896 ; +18,6896 |

## tab:tests_M par famille

- famille 3m (6x6, 9x9, 12x12) : (max − min)/moyenne = 2.257e-02 ; csv « 2.3e-02 » (OK)
- famille non-3m (5x5, 7x7, 8x8) : (max − min)/moyenne = 2.057e-02 ; csv « 2.1e-02 » (OK)
- max|M| (denses et grossiers) = M2 : True ; ligne « six tailles » absente : True

Portes D6 d'`analyze_M.py` (C.1) : align_gate_closed_vs_wannier 2.9e-15 ; align_gate_coarse_wannier 1.3e-15 ; align_gate_frozen_bloch 1.4e-10 ; align_frozen_projector 1.4e-10

Révision du seuil de la porte P-c2 (Greg, 2026-09-30) : 1e-10 eV au lieu de 1e-13 ; réévaluation sans recalcul sur les écarts enregistrés par le job 22058857 (max 1.21e-13 eV, 6×6 masqué) : 16/16 profils OK ; appartenance aux anneaux inchangée (8/8 OK).
