# Données de la chaîne 5×5 (export du 2026-10-08)

Copie des données de rorqual (rien n'a été modifié à la source). Le code et `config/production.json` viennent du dépôt `graphene-raman` (git).
`MD5SUMS.txt` : md5 de chaque fichier, calculés sur les sources ; vérifier avec `md5sum -c MD5SUMS.txt` depuis ce dossier.

## Contenu
| Chemin | Source sur rorqual | Rôle |
|---|---|---|
| `data/graphene/supercell/qe/defect_5x5_{p,d}.save/` | `qe_tmp/defect_5x5_{p,d}/…save` + `super_cell/5x5/{pristine,defective}/Vks_5x5_{p,d}` | paire de super-cellules (Γ) et leurs potentiels pp.x |
| `qe_tmp/defect_uc_dense_25/defect_uc_dense_25.save/` | `qe_tmp/defect_uc_dense_25` (nbnd 20, 625 k) | maille unitaire dense de la chaîne de production pour la 5×5 (`config.dense["5x5"]` : p = 5, D = 25) |
| `wannier/25x25/` | `graphene-raman/wannier/25x25` (= `unit_cell/25x25`, construit sur `defect_uc_dense_25`, 20 bandes) | Wannier de la chaîne de production (manifeste sha256) |
| `data/graphene/unit_cell/qe/defect_unit_cell_5x5_wann.save/` | `qe_tmp/defect_unit_cell_5x5_wann` (nbnd 16, 25 k) | maille unitaire 5×5 sur laquelle `wannier/5x5` a été construit |
| `wannier/5x5/` | `graphene-raman/wannier/5x5` | Wannier grossier 5×5 |
| `data/graphene/unit_cell/qe/defect_5x5.save/` (+ lien `defect_unit_cell_5x5.save`) | `qe_tmp/defect_unit_cell_5x5` + `unit_cell/5x5/Vks_uc_5x5` | maille unitaire vers laquelle pointe `data/` sur rorqual (M grossier de production `M_ed_5x5`, `test_ks_reconstruction`) |

## Deux ensembles cohérents (même jauge de Wannier)
1. **Production** (niveau 1, matrice T) : super-cellules + `defect_uc_dense_25` + `wannier/25x25`. M dense par `scripts/compute_M_dense_stages.py`
   (étapes ml, nl, combine), puis les pilotes qui lisent `config.dense_paths` (par défaut `scratch=/home/gregb26/links/scratch/qe_tmp` :
   passer le chemin du `qe_tmp/` de cet export).
2. **Grossier** : super-cellules + `defect_unit_cell_5x5_wann.save` + `wannier/5x5`.

**Attention** : `defect_5x5.save` (= `defect_unit_cell_5x5`) et `defect_unit_cell_5x5_wann` sont deux runs distincts (fonctions d'onde
différentes). `wannier/5x5` n'est pas dans la jauge de `defect_5x5.save` : ne pas les combiner pour une interpolation de Wannier.
`wannier/5x5` et `wannier/25x25` sont aussi suivis par git (fichiers identiques).
