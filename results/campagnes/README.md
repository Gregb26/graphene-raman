# results/campagnes — données de campagne lues par les scripts de figures

Copies (`cp -p`, 2026-09-30) des json/npz/csv de `campagnes/` que lisent `scripts/fig/make_figures_{controles,em,kaasbjerg}.py`,
au même chemin relatif que dans la campagne, pour que tout ce qui produit une figure du mémoire se lise dans `results/`
(`config.campaigns_dir()`, clé `campaigns_results_dir` de `config/production.json`). Rien n'est calculé ici ; la campagne
(répertoire, README, pilote, rapport) reste la source, et on la recopie si elle change. Sommes : `MD5SUMS_2026-09-30.txt`
(`md5 -r`, vérifiées égales aux originaux à la copie).

| Sous-dossier | Source (`campagnes/`) | Fichiers | Lu par | Figures |
|---|---|---|---|---|
| `R7_tailles_3m/d1_8pts/` | `R/R7_tailles_3m/d1_8pts/` | `d1_results.json` | `make_figures_controles.py` | `fig_size_3m`, `fig_localized_3m` |
| `R9_controles/{a,b,c}/` | `R/R9_controles/{a,b,c}/` | `a3_results{,_plateau}.json`, `b_results{,_plateau}.json`, `c_results.json` | `make_figures_controles.py` | `fig_resonance_vs_nkint{,_plateau}`, `fig_rcut_aligned`, `fig_folded_vs_R7` |
| `R10_plateau/{a,c}/` | `R/R10_plateau/{a,c}/` | `a1_results.json`, `profiles_<S>.npz` (13 tailles), `a2_results.json`, `c2_kaasbjerg_maps.npz` | `make_figures_controles.py` | `fig_offset_profiles_13`, `fig_levels_vs_invN`, `fig_kaasbjerg_plateau_ws` |
| `R8_kaasbjerg/{dos,spec,sens,7a,sigeff}/` | `R/R8_kaasbjerg/out/{dos,spec,sens,7a,sigeff}/` | `dos_9x9.npz`, `spectral_GKM_9x9.npz` (14 Mo, rapatrié de rorqual le 2026-09-30 ; absent de la copie de campagne), `sens_results.json`, `fig13_VA_interp_1meV.csv`, `7a_results.json`, `sigma_K_9x9.npz`, `sigeff_results.json` | `make_figures_kaasbjerg.py` | `fig_kaasbjerg_{dos_c,spectral_GKM,sensibilites,superposition,sigma_K}` |
| `EM/{M4_sigma,EM3,EM2}/` | `EM/{M4_sigma,EM3,EM2}/` | `em_scalars.json`, `em_sigma_{full,centres_only,no_berry}_N1200_eta0.04.npz`, `em3_ring_2p33.npz`, `em2_postw90_sigma_ti.npz`, `em2_A.npz` | `make_figures_em.py` | `fig_em_coupling` |

Les autres produits des campagnes (rapports, tables, gates, sorties intermédiaires, inputs QE) ne sont pas copiés : voir
`campagnes/README.md` et le README de chaque campagne.
