# results/M2_plateau — produits de production du chapitre 4, base plateau (i) (R10, 2026-09-30)

- Base : ΔV aligné par la moyenne du plateau (i) (config/production.json, bloc `alignment` ; C_N des 13 tailles, source `article/R10_plateau/a/a1_results.json`) ; M_W(R, R) − C_N·𝕀₅ sur les mailles de la boîte N×N (approximation (i), `local_tmatrix.defect_mwr`) ; interpolation hors grille par les images de Wigner-Seitz (`ws_images`, `ws_phase`) ; corrections de l'audit de l'image minimale (P-b2, P-c2, P-c4).
- Matrices : `results/M2` (M2 brutes, lecture seule) ; ce répertoire ne contient que des produits (npz, csv) et les journaux (`logs/`, non versionnés).
- Rejeu : campagne R10 (`graphene/qe/defects/R10_plateau/`, copie `article/R10_plateau/`), lanceur `submit_r10.sh` (tâches c0, c1, c1f, c2, c3) ; HEAD d262591.
- Référence v2 tel quel : `results/M2` (gelé) ; table de correspondance : `article/R10_plateau/c/table_v2_plateau.md`.
- `ks_reconstruction.npz` copié de `results/M2` (indépendant de M ; md5 identique), `ved_analysis.npz` refait (P-c2).
- md5 : `MD5SUMS_2026-09-30.txt` (33 fichiers).
