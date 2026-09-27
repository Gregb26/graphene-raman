#!/bin/bash
# R6 étape 3 — ordre de soumission après le test d'or (3.0b PASS). À exécuter depuis la racine du dépôt (cd $GRAPHENE_RAMAN).
# A. niveau 1 (3.3) : six tailles, grilles 60/120/240, eta 0.05/0.02/0.01, R_cut 0..4 ; tag prod
for S in 5x5 6x6 7x7 8x8 9x9 12x12; do GRIDS=60,120,240 ETAS=0.05,0.02,0.01 RCUTS=0,1,2,3,4 sbatch --job-name=specwd_$S scripts/submit_spectral_wannier_dense.sh $S prod; done
# B. N_k^int 150..600 et Re Sigma par R_cut (9x9)
sbatch scripts/submit_nkint_check.sh 9x9
sbatch scripts/submit_rcut_resigma.sh 9x9 0,1,2,3
sbatch scripts/submit_rcut_resigma.sh 9x9 4
# C. 3.1 / 3.2 / 3.6 (indépendants du niveau 1)
sbatch --job-name=post_loc scripts/submit_post.sh locality
sbatch --job-name=post_ana scripts/submit_post.sh analyze
sbatch --job-name=post_ks  scripts/submit_post.sh ksrec
# --- STOP 3.3 : porte de niveau 1 (r6_level1_gate.py) sur specwd_9x9_prod + nkint + resigma ; si échec, ne pas lancer D ---
# D. 3.4 / 3.5 (après la porte)
sbatch --job-name=post_res9  scripts/submit_post.sh resonance 9x9
sbatch --job-name=post_res6  scripts/submit_post.sh resonance 6x6
sbatch --job-name=post_res12 scripts/submit_post.sh resonance 12x12
sbatch --job-name=post_c14   scripts/submit_post.sh c14
# E. figures et tables (après tout)
sbatch --job-name=post_fig scripts/submit_post.sh figures
