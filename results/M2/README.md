# results/M2 — M corrigés (v2, campagne R6)

M2 = N_cells·M^L + M^NL, en Hartree, M^L et M^NL tous deux en norme `unit_cell` (sidecar `M_normalization = v2`). Remplace
`results/M/` (v1, gelé : M^L sans le facteur N_cells). Rapport : `graphene/qe/defects/R6_production_corrigee/R6_rapport.md`
(copie `article/R6_production_corrigee/`).

**Pas de miroir** (décision de Greg, 2026-09-27) : les 43 fichiers M (`*.npy`) se reconstruisent en environ 11 min à partir de
fichiers déjà sur `/project`.

## Dépendances à conserver

- `results/M/` : parties L et NL v1, `M_ed_*` de juin, `M_L_dense_9x9_coarsecheck.npy`. Les `M_NL_dense_*` d'ici sont des **liens
  symboliques** vers `results/M/` (20,8 Go) : ne pas déplacer, archiver ni supprimer `results/M/` sans remplacer d'abord ces liens
  par des copies.
- `graphene/qe/defects/R6_production_corrigee/ml/M_L_{5x5,6x6,8x8}_v2.npy` : M^L grossiers recalculés avec le noyau corrigé (24 Mo).
- `graphene/qe/defects/R5_base_vs_M/b/M_{L,NL}_9x9_nb128.npy` : 128 bandes (3,4 Go).

## Reconstruction

```
cd graphene/qe/defects/R6_production_corrigee && sbatch submit_r6.sh assemble    # scripts/assemble_M2.py ; 16 cœurs, ≈ 11 min (job 21820490)
cd $GRAPHENE_RAMAN/results/M2 && md5sum -c MD5SUMS_2026-09-25.txt                # attendu : 43 OK
```

`assemble_M2.py` refuse d'écraser un fichier existant (ou son sidecar) sans `--force` : ne reconstruire que ce qui manque.
43 OK = fichiers identiques au bit à ceux qui ont passé la porte A.2 (job 21820491) et servi à toute la production R6. Sinon,
repasser la porte : `sbatch submit_r6.sh gate` (demande les `.save` de maille du scratch).

Produits dérivés (`specwd_*`, `resonance_*`, `*.csv`…) : suivis par git (exceptions `.gitignore`) ; ils se refont avec
`scripts/submit_post.sh` et les autres `submit_*` (enchaînement : `R6_production_corrigee/etape3/runbook_3.sh`). Non suivis :
`resigma_9x9_*.npz` (`submit_rcut_resigma.sh`, `submit_nkint_check.sh`) et `logs/`.
