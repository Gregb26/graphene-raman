# R10 — A.1 : C_N = moyenne du plateau (i) pour les 13 tailles (tables générées par r10_driver.py)

Plateau : atomes de la cellule avec lacune à distance vraie (image minimale) ≥ 0,75 r_max de la lacune ; sphères de P1 (`alignment.atom_sphere_shifts`) ; ΔV = ⟨V_d⟩ − ⟨V_p⟩ ; « un atome vrai » : atomes à r_max en vraie image (égalités à 1e-6 Å), moyenne en cas d'égalité ; Lu publié : atome de `far_atom` (réduction axe par axe), R5 C (5…12) et R7 D1 (15…27). meV sauf mention.

Porte A.0 (six tailles de R9) : 5×5 C_i au bit oui, un atome [2] (écart max à l'audit 0.0e+00 meV) ; 6×6 C_i au bit oui, un atome [2, 57] (écart max à l'audit 0.0e+00 meV) ; 7×7 C_i au bit oui, un atome [81] (écart max à l'audit 0.0e+00 meV) ; 8×8 C_i au bit oui, un atome [20] (écart max à l'audit 0.0e+00 meV) ; 9×9 C_i au bit oui, un atome [21, 140] (écart max à l'audit 0.0e+00 meV) ; 12×12 C_i au bit oui, un atome [28, 235] (écart max à l'audit 0.0e+00 meV) → OK

| N | r_max (Å) | plateau (i) 1,0 Å : moyenne ; rms ; max\|écart\| ; n | plateau 0,5 Å : moyenne | un atome vrai : indices ; valeurs 1,0 Å ; moyenne | Lu publié : atome ; valeur | C_N − Lu | ⟨ΔV⟩_3D |
|---|---|---|---|---|---|---|---|
| 5×5 | 7,12 | −57,56 ; 94,14 ; 233,70 ; 13 | −59,98 | 2 ; +176,14 ; +176,14 | 49 ; −122,33 | +64,77 | +28,62 |
| 6×6 | 8,54 | −50,51 ; 24,56 ; 50,04 ; 26 | −52,33 | 2, 57 ; −100,55 ; −86,57 ; −93,56 | 71 ; −98,62 | +48,11 | +20,83 |
| 7×7 | 9,97 | −26,91 ; 18,53 ; 35,58 ; 31 | −27,81 | 81 ; −37,04 ; −37,04 | 97 ; −38,54 | +11,62 | +14,55 |
| 8×8 | 11,39 | −14,41 ; 18,87 ; 59,19 ; 49 | −14,58 | 20 ; +44,78 ; +44,78 | 1 ; −22,97 | +8,56 | +11,23 |
| 9×9 | 12,81 | −25,14 ; 9,05 ; 24,10 ; 53 | −25,29 | 21, 140 ; −42,52 ; −49,24 ; −45,88 | 161 ; −24,65 | −0,49 | +9,09 |
| 10×10 | 14,24 | −9,34 ; 6,13 ; 12,54 ; 70 | −9,46 | 23 ; −14,29 ; −14,29 | 198 ; −8,22 | −1,12 | +7,17 |
| 11×11 | 15,66 | −9,55 ; 5,76 ; 17,57 ; 73 | −9,56 | 26 ; +8,02 ; +8,02 | 241 ; −12,16 | +2,61 | +5,95 |
| 12×12 | 17,08 | −18,69 ; 5,42 ; 13,88 ; 95 | −19,08 | 28, 235 ; −32,57 ; −28,66 ; −30,61 | 287 ; −29,24 | +10,55 | +5,06 |
| 15×15 | 21,35 | −16,22 ; 2,76 ; 8,30 ; 149 | −16,22 | 65, 384 ; −22,16 ; −24,52 ; −23,34 | 449 ; −15,64 | −0,58 | +1,42 |
| 18×18 | 25,63 | −15,06 ; 1,74 ; 4,99 ; 212 | −15,18 | 115, 570 ; −18,58 ; −20,06 ; −19,32 | 1 ; −18,18 | +3,12 | −8,48 |
| 21×21 | 29,90 | −14,23 ; 1,03 ; 3,13 ; 275 | −14,33 | 133, 748 ; −16,38 ; −17,37 ; −16,87 | 881 ; −13,77 | −0,46 | +1,03 |
| 24×24 | 34,17 | −13,48 ; 0,72 ; 2,25 ; 383 | −13,56 | 201, 1000 ; −14,99 ; −15,73 ; −15,36 | 2 ; −13,12 | −0,36 | −62,74 |
| 27×27 | 38,44 | −13,01 ; 0,53 ; 1,62 ; 461 | −13,06 | 225, 1232 ; −14,09 ; −14,63 ; −14,36 | 1457 ; −12,77 | −0,24 | −99,91 |

Contrôles par taille : P1 à l'atome de Lu = `far_atom_alignment` au bit ; écart au Lu publié :

- 5×5 : atome 49 (publié 49, R5 C) à 9,97 Å axe par axe, 6,21 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
- 6×6 : atome 71 (publié 71, R5 C) à 12,81 Å axe par axe, 7,40 Å vrai ; P1 = fonction : oui ; écart au publié -7.1e-12 meV
- 7×7 : atome 97 (publié 97, R5 C) à 14,24 Å axe par axe, 8,66 Å vrai ; P1 = fonction : oui ; écart au publié -7.1e-12 meV
- 8×8 : atome 1 (publié 1, R5 C) à 17,08 Å axe par axe, 9,86 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
- 9×9 : atome 161 (publié 161, R5 C) à 18,51 Å axe par axe, 11,12 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
- 10×10 : atome 198 (publié 198, R5 C) à 19,93 Å axe par axe, 12,41 Å vrai ; P1 = fonction : oui ; écart au publié -1.4e-11 meV
- 11×11 : atome 241 (publié 241, R5 C) à 22,78 Å axe par axe, 13,58 Å vrai ; P1 = fonction : oui ; écart au publié -7.1e-12 meV
- 12×12 : atome 287 (publié 287, R5 C) à 25,63 Å axe par axe, 14,80 Å vrai ; P1 = fonction : oui ; écart au publié -7.1e-12 meV
- 15×15 : atome 449 (publié 449, R7 D1) à 31,32 Å axe par axe, 18,51 Å vrai ; P1 = fonction : oui ; écart au publié -7.1e-12 meV
- 18×18 : atome 1 (publié 1, R7 D1) à 38,44 Å axe par axe, 22,19 Å vrai ; P1 = fonction : oui ; écart au publié 7.1e-12 meV
- 21×21 : atome 881 (publié 881, R7 D1) à 44,13 Å axe par axe, 25,90 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
- 24×24 : atome 2 (publié 2, R7 D1) à 49,83 Å axe par axe, 29,62 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
- 27×27 : atome 1457 (publié 1457, R7 D1) à 56,95 Å axe par axe, 33,30 Å vrai ; P1 = fonction : oui ; écart au publié 0.0e+00 meV
