# Réécriture de l'historique git — 2026-09-29

L'historique de `main` a été réécrit le 2026-09-29 (décision de Greg) :

- lignes `Co-Authored-By: … Claude …` retirées des messages de commit (91 commits) ;
- courriels `gregb26@rorqualN.rorqual.calcul.quebec` remplacés par `gregoire.barrette@umontreal.ca` (auteur et committeur, 99 commits) ;
- contenu des commits (arbres), noms, dates, parents et reste des messages inchangés, vérifiés commit par commit ;
- les 5 commits faits par l'interface web de GitHub perdent leur signature (entête `gpgsig`) ;
- branches supprimées : `cleanup-M-driver` (locale, jamais poussée, 5 commits non fusionnés) et `scripts-improvements` (GitHub, entièrement contenue dans `main`).

238 commits, dont 203 changent de numéro. Un ancien numéro désigne exactement le même contenu que le nouveau numéro en regard.

Ancien historique complet (toutes les branches, GitHub et Rorqual) :
`/nearline/rrg-cotemich-ac/gregb26/graphene-raman_avant_reecriture_2026-09-29.bundle` (md5 `2f6e0a1664511daa642b537c75f0e155`).
Pour le consulter : `git clone <bundle> ancien` ; les branches de Rorqual y sont sous `refs/rorqual/heads/*`.

Citations mises à jour dans le dépôt (et dans les originaux identiques hors dépôt) : documents (`*.md`) et code vivant
(`scripts/`, `requirements.txt`). Citations laissées telles quelles, à traduire avec cette table : sorties de programme
(champs `head` des json, journaux) et code des campagnes déjà exécutées (`article/R6_production_corrigee/etape*`,
pilotes R7 et R9), dont des md5 sont enregistrés ailleurs ; journaux SLURM et `submitted/` hors dépôt.

## Table ancien → nouveau (ordre chronologique)

| ancien | nouveau | date | sujet |
|---|---|---|---|
| `6e371685af068beb4677e5d4721f30630238cd2d` | (inchangé) | 2025-09-19 | fresh start |
| `3bd9a3d2604e3a81c77707d9b8c084f5ad352158` | (inchangé) | 2025-09-19 | added utilities functions |
| `3f13e0bd81680fd4fd3dba5694c55faee9d35e2f` | (inchangé) | 2025-09-19 | added defect matrix elements |
| `1ab25fab5e655e2764082a13304f5d4d72d9f301` | (inchangé) | 2025-09-19 | added wfk.py |
| `23ab9d7e3173f9e76116987c206e9bbbc9ef24e1` | (inchangé) | 2025-09-19 | added wavefunctions |
| `dbb7ba0c1b3ec5014647ee01c2c0db9f2885cec1` | (inchangé) | 2025-09-19 | added wavefunctions |
| `2d3d03b24df2fe6e68e99b8ad8cd355106cc2104` | (inchangé) | 2025-09-19 | small corrections |
| `57885950e05de569b32625425d3a65ffcabc3a49` | (inchangé) | 2025-09-19 | small corrections |
| `844d41ee514ddbe8a6bda6109627b9f0745f4385` | (inchangé) | 2025-09-19 | small corrections |
| `534d6a34825c94af85e073a01834ab48ebc22236` | (inchangé) | 2025-09-19 | small corrections |
| `e4241d5aff54b22b3900911d9cba1e409228bd78` | (inchangé) | 2025-09-19 | fresh commit |
| `6b7c9e0b0cb4ebdce0b38a8b9e9ad0d9d2881147` | (inchangé) | 2025-10-29 | added fft utils functions |
| `2b770951b62aa322dc05d1e3ac49d4d803b57e3b` | (inchangé) | 2025-10-29 | added utils functions |
| `6c9a3794345f6014fec920c82967f760522a6584` | (inchangé) | 2025-10-29 | added map from G vectors to FFT indices |
| `b3c99e727e70eb001e4e389f90cb6f2404f94228` | (inchangé) | 2025-10-29 | added map from G vectors to FFT indices |
| `03fb111ecebe6fd3a590b2dafaeec55af7769998` | (inchangé) | 2025-10-29 | added a function to compute the wavefunctions from pw coeffs |
| `8c5bd8b56e96b487c8e22eda68eebf94fa7a686e` | (inchangé) | 2025-10-31 | removed deleted files |
| `75a77bf1a4dcc9605ad033e3191e213489ee1010` | (inchangé) | 2025-10-31 | added utils functions |
| `1a9bf0f636fc3af858c118b450707c2005aa4a79` | (inchangé) | 2025-10-31 | added analysis |
| `3468ba93ff61ea1d8bde078ca210208e1be7d567` | (inchangé) | 2025-10-31 | removed deleted files |
| `51054740c44a66f3101c2f54d20ecf3bfcbe41f2` | (inchangé) | 2025-10-31 | fresh commit |
| `6341fafa35d179bdfd369323d387e87e12d63b61` | (inchangé) | 2025-10-31 | added kinetic energy function |
| `473d0a4bbf344c6d8f1910790331503bb754a879` | (inchangé) | 2025-11-02 | adapted read_psps to any kind of .psp8 file |
| `8d7895e0223fe623f772d3f8edaf3bced061e9b7` | (inchangé) | 2025-11-02 | added a function to plot radial projectors |
| `e20599616db219ae378f483ee41216b436aa38e7` | (inchangé) | 2025-11-13 | optimized local matrix elements calculation |
| `a434d4503b8e69c91eb7f394e3b44a2c288c6540` | (inchangé) | 2025-11-13 | streamlined pseudo reader and added Wannier reader |
| `7cda2f5174161e8cc526bcabe0a808f21431bbe9` | (inchangé) | 2025-11-13 | added plotting functions |
| `67a558d97d66d4977a93651c1fef04d4c4ad5584` | (inchangé) | 2025-11-13 | fixed bugs |
| `395061ebac52f3fed986cb93dabfe516df5b7590` | (inchangé) | 2025-11-13 | added functions to do Wannier interpolation |
| `7056f7f13a31331bab11d85ca2d2ccc1a79581bf` | (inchangé) | 2025-11-13 | renamed non local file |
| `ae9d82aa9ad9dc5d7a9b4ffb968508d72918d664` | (inchangé) | 2025-11-13 | fixed a bug |
| `dc6cbf9e7fd56a03aed51698651fa489c54c6ac7` | (inchangé) | 2025-11-13 | refactored |
| `b1a461264859d761b2c32ba3bee14fe32cde9cc2` | (inchangé) | 2025-11-13 | added many body quantities for single defect case |
| `1d63c794da22a96124abce248a6aad0307f6f0e1` | (inchangé) | 2025-11-26 | modified source folder |
| `be45124a808668656d62ee2bb8b3d1d61beaf214` | (inchangé) | 2025-11-26 | added Wannier-interpolated band structure |
| `1b167e8f36c4c01b4e08adc26937aacdd7ca1050` | `4f65abb8a6073be4b10cd232dc51507ef2c43735` | 2025-11-26 | Add project description to README |
| `99d1adab0758e266b29c35dc04477b3ce7e2d777` | `106a26ea850d6b998321c1a8560de0978bd5ee2f` | 2025-11-26 | Merge branch 'main' of https://github.com/Gregb26/ab-initio-electron-defect-interaction i… |
| `7709199bf066208c3a3d7988be331aadc2aa5f79` | `4339ad01f0af107d85467f3a8a0f454922164d24` | 2025-11-26 | added comments |
| `9ce0fda3ee3a1c513e22636e7053b17dc6e484c6` | `d3aead025a9e227787b0e09109b173cb4bd00396` | 2025-11-26 | fixed a bug in get_Cnk |
| `87461a455994f50026b4e51af045b602ebb1b6f8` | `2c969956e8aa593172b5409265ae5aa1d8f1a790` | 2025-11-26 | cleaned up fft_utils,py |
| `d4d794e8fccb488526fff6dff5b28f4e957b05f1` | `4dfe83174fdf3e0dfe380c205f9e1f5416f73f2d` | 2025-11-26 | added progress bars to computation of wavefunctions |
| `9f4248944fa3075662e1935a3fc46e06e0049562` | `dace96892b5e75a59dd7c753bcd1135ca1da8422` | 2025-11-26 | fixed a bug in local_R.py |
| `919a4797d13bcef5cc7e017d018405ac562992ba` | `f296432369075a14de73baa833de9585c00407f0` | 2025-11-26 | added requirements.txt for venv |
| `6f8d9cb840c1a0a0ebbfe33f41120c3439c9369c` | `5c1bb5d2f999ce51ba7b3c0d3d33c63478b51a22` | 2025-11-26 | added a small script to compute local matrix |
| `6b7ffcac033364041c671d77dce790c866348254` | `9b7ea8fdbd3adaf58185f65b6071c503d61ffdd2` | 2025-11-27 | modified requirements.txt and .gitignore |
| `457b93d578304dbac28df8c6caa8f2ebf542dd57` | `ba25c96629f1c7e1614ceeb6b270be2ac4392dc9` | 2025-11-27 | modified requirements.txt and .gitignore |
| `7182c63b64174a5b1b2ef8c16fe7a2aa12e371b2` | `439314eb8ae27325c1a4c9ae7b0ba22e6eb38bb3` | 2025-11-27 | modified requirements.txt and .gitignore |
| `1c188a252a3820fd821aa0f3e3a925ed5163a742` | `d6f8a3f4c6857fdafaebdca271e50c6dbd428f52` | 2025-11-27 | cleaning |
| `75afeac27775be48b55e4482b626eb59467fb7d0` | `3fd144d6befadbf845b237dfa7c1d598aa87c663` | 2025-11-27 | added first python script |
| `680ee46ec6c0660c13651169786d0c9600717843` | `cea32a23b0a2628cc039b29a1dc0a3deae1fa404` | 2025-11-27 | modified .gitignore |
| `1d81d04d87d4de3f94cae13af5eeb7153e83c4b9` | `d3d1186cb791f63f4508633b58247657f9d91383` | 2025-11-28 | added a script to read the hamiltonian from wannier90 |
| `f9866f22cc511ff4cffc2158edbc84b2eac866b1` | `9debfe8fd097a9ea49a9c9ca00a195594eca9676` | 2025-11-28 | added TODO |
| `82fd16781abbe42600817b38d99c0ae705e68a94` | `2fb523a39d0fc586687b3eb4bbe4628f919d41e0` | 2025-11-28 | added a double FT function |
| `d7b390d1327e69471c8bc9700f202efd7c9247c1` | `6d4d0d0e4da09fd5807e61602eecd2fa9c28405f` | 2025-11-28 | added a function to generate a MP grid |
| `9bf3fa50a4ab00cef54bc66160cf3b054d76bdeb` | `9d880a82369dbc708d7e54c84f8488d869989853` | 2025-12-09 | added a function to computet the potential in the unit cell |
| `08ad080d86c7c2853d3883fe76d2add24e13cba2` | `256866428f5bf4813482fe9e49a5218f3a555c15` | 2025-12-09 | added Wannier Hamiltonian interpolation |
| `342d6d24f07dfc351c6c317f054ca76089223011` | `c1e6571ddfd2065a3534e6f2a908615ad5867905` | 2025-12-09 | added functions to interpolate matrix elements using Wannier functions |
| `4cbdcbee0f81da9d9e4eb3db50fee5c3c66ebddd` | `f86dbf9e3d18e5b663a668fc9e9e690bab713906` | 2025-12-09 | added normalization check |
| `09eebd91bba2b8229ba6ef97e4287aa9577a1454` | `b38bba205853e5cb647448604c052c28d923149e` | 2025-12-15 | added a script to compare DFT and Wannier90 bands |
| `34c4a0a535e057a34ad0b2cb5e0b70501eb9e666` | `3fe7ef41c4ab1960e47c95747182a66fe14aa5a9` | 2025-12-15 | fixed docstring |
| `a6a5cbfb2504ba7ea6e93793a2697dd3ec7b478d` | `0d4b1d60f74fe54c0ff0987dcdbe54cd8d09867c` | 2026-01-23 | efficient computation in real space with MPI |
| `7e7712a4cfd9f143312019093c27313a68a681cd` | `9dc141ce4cb50677e981f72ce709916a64354b2c` | 2026-01-23 | fixed get_eig() function |
| `32dbc4bbc97e70d5c4f115ab4086d559737bee34` | `ea45be7ec8fbbd7165a4f0addb1bf4979092562b` | 2026-01-23 | changed to double complex precision |
| `a5e55599b0b707ba345476a45598fee553fc2e48` | `abda79f0c3b3130a4956532285cf3d0e90f8e2cb` | 2026-01-23 | added script to use mpi |
| `8906fcac5fd97106cdc2cc5148e55cdadf10111c` | `40bf163bdb628f6960e9c7e300767c6145535f74` | 2026-01-23 | fixed mpi script |
| `88e447b9d2369b105dad27b4cd8f5bdf7dae7733` | `ce659923ad8514125df7690bfd9fd5d061876a08` | 2026-01-23 | refactored non mpi script |
| `b56f699f5a5487e891ee002ccae66a255d3a32bd` | `9f0a8344b791e65eb402db26c327cacaae34e05b` | 2026-01-23 | Remove obsolete script |
| `7fd11a23a27f59623775037ea2a0bac8a8d880e9` | `2a90814bdf47747d45ad5559be1215d8edea2e75` | 2026-01-23 | fixed bugs |
| `073a0ec82f4b2eb0c0838ada89d449b9d5511986` | `fb5d9bbacc5efaee0df73c839f08b1957c704378` | 2026-01-23 | fixed double complex MPI attribute |
| `966129761b4d39c2c06ab99c6a5032083fb63f00` | `894ed7c790aabd623d806b18e2ddf818d1a334c6` | 2026-01-23 | added prep on rank 0 |
| `41356499ae38f6b00fcb72ff5889566d54a166c7` | `fc7d970a291dad9e66d943c201c0c4ae4232968a` | 2026-01-23 | added timing |
| `94ac4dcea1e35abd6454c9adcb0411e1040e5364` | `48d25985c275c99296cba20cf47652087ab53059` | 2026-01-23 | added timing |
| `14425ac8c28835fb995d0342c996c337a974f323` | `0f7952bcca4eb6043b22ca9d13584cbf8e84250f` | 2026-01-23 | added timing |
| `f2fa50763e3be3dc76af791b7ca1ec5474bf7c3c` | `3188bb793bfe7ccaf24ebf1258ca91c7074e09bc` | 2026-01-23 | removed timing |
| `d51e0c4eaec8196d2a8b44e25582faaf349383d5` | `a3d747c88c7b9f2787a273071fc539d94f16ab61` | 2026-01-23 | test |
| `2323d840479b735f98bc02aa9b02792c2dfb3934` | `d6313f6b163f793304dcea0a312f7e5051a5470c` | 2026-01-23 | optimized non-local part |
| `0ba863d42ba57f29a025426c33e54438354457a2` | `0365298d9301036deaaff2b27b5109b792687f40` | 2026-01-23 | optimized non-local part |
| `294776208512db028d66685e6f1a7e5d3cdbff24` | `1bdd8f0d3942fe7459bbc8b03251cbbe76f3b8c3` | 2026-01-23 | memory friendly |
| `570f606afdc0a31e004c2d00bfd435bc227b1965` | `766328421b5bb76bbd2ab1467b64f3106543c21b` | 2026-01-26 | added a function to build a path in kspace |
| `5c5bc56eb9d2de5efdfcf4754afbada5a34e9078` | `119bf320d6da336ebb7638bdfbd14f3ffed42da5` | 2026-02-02 | added trilinear interpolation for fourier transformed potential |
| `1f34397a9ff0cb529c242ee07250ec4cccfa5a97` | `e4e36e9c0db21d848a307dca7f146d095a5ea2bb` | 2026-02-02 | added trilinear interpolation |
| `487880ce017993126a352f7b8aba8d4871f1d9a3` | `5aed0040867bf2c1a111bc386216f1adeadab821` | 2026-02-02 | refactored |
| `8ebcf5ff93b2afb1845388bd526bbdc350cf380c` | `fd5a323b3578164e7ce9305d2a161eb3b314f3c3` | 2026-02-02 | added MPI for interpolation case |
| `b9b703e48859d3a05b5318ab85bf542ebf59ddfa` | `24a679708b028d69b17fa0b05ec07932588956c9` | 2026-02-03 | started to add MPI for non-local part |
| `21dfa6955b90062c2d8729243953fbc97cde47d6` | `85980f9d75e6c4c2dee0d73ed33ebd3131716d8a` | 2026-02-03 | fixed a bug |
| `b8ddce8784e3f76be8a97d0be90cc2ccf0be39d9` | `4a4915c39fae5d4a508c77d436e5e7c986e42c4d` | 2026-02-03 | added periodic cubic spline interpolation of 3d array |
| `1d43ae616b68ffb54b0cac6fe2b519f738b497c1` | `c61aa822ec259f827faf841e0c994d9ce73d5e14` | 2026-02-03 | replaced trilinear for more accurate cubic spline |
| `cf77b9bbe0ce1714cec2b64ae78fc26b95b6c638` | `d8a051befa159c1eb6db7a1fa9724c0563bdb09d` | 2026-02-03 | fixed a bug with cubic spline interpolation |
| `59d374ac3408b2de68936cfbafd89848f70b04c5` | `a84fa35fb35a7a21234aa30f5090305a45f54d79` | 2026-02-03 | fixed a bug with cubic spline interpolation |
| `f1470ee807d677cf65d397973451eac73204e031` | `ef658205edc391c0c31b5302f43f4ea45c5f0fae` | 2026-02-04 | fixed a bug |
| `69d2773e50f876419ff17762e921a966999f6c84` | `85f167e8ed73d99621b38c5ff231c5dfbb1d27f6` | 2026-02-04 | fixed a bug |
| `bc1e598c58c91c2d46ddbfb5a9b91b97aa4183a5` | `fd1c6f618b95372da172d638b44d88490a0a6373` | 2026-06-17 | Port electron-defect matrix-element pipeline to Quantum ESPRESSO |
| `90500a07feb8e853942709b9ff35e25abf768235` | `a167c6082667634b565f7654268d00ec2c6f182a` | 2026-06-17 | Add standalone QE pipeline validation script |
| `eeed43656d2b2b3e8859adda1fb3a5a7a614af39` | `173fa36a51c6c27bbdc9b477dcd736c50f3d7f0e` | 2026-06-17 | Add cluster driver to compute full M = M^L + M^NL (MPI) |
| `266794b4d4f8cd45fe0aa66326f92c08d38eb5d0` | `ecb11b97bdd5144eb7a287c4bf8f1d01782e2cd9` | 2026-06-17 | Add grid-distributed MPI real-space M^L and --method switch |
| `e73a8a8096b73d566e0a514930244eef2930c668` | `3cd1fd80c21c50696f4605d307607e821ed4d1a6` | 2026-06-17 | Parallelize M^NL over k' and add band restriction |
| `025a3287e0e80e521e3464c43a49277dcac6584b` | `dae9849eb826cf1aaa4d6efc70747e6e94707bd6` | 2026-06-18 | up to date |
| `cab1e3ce8603a57560266924aa5aab1432d55a09` | `f51a514f30492916ed80d1395cd4c8bb2f7d4459` | 2026-06-18 | refactored local_R to use QE inputs |
| `56c5eb852c91727f319804bbd9a3f9c7a8ef1d00` | `2d04b95e7095faf1f95876dadcaf9d87fc1b6b08` | 2026-06-18 | refactored local_G to use QE inputs |
| `a1e53e5e0852508c4c27389dbb91f94da1e5b6b6` | `7a6266eb82d6099a3d929c55a9722b9587087116` | 2026-06-18 | refactored non_local to use QE inputs |
| `8e167398acfdd13a38f8b69ed70a1387dc62bc3c` | `4720211bae0be6ef7d07403f46066093db3a7459` | 2026-06-18 | removed reference to abinit |
| `615283396974d752e0a6f5cff72e8a0013e9d09a` | `b8dedc1259d297f0c30c4e0c1d278a5e576089d2` | 2026-06-18 | refactored to remove reference to ABINIT |
| `ccc5a9ace65237f8cbee08fcd36fbaca205369af` | `5e7e41eba342cc2069c024a61db9c46b9da3f3ae` | 2026-06-18 | finished refactoring |
| `3ee4efae6cdde0000d6f75cc9075d882791f596c` | `59a783b0f41003b978b4a20c86f7a5167ac2e296` | 2026-06-18 | remove deleted files |
| `5c7f507d8bd4c1b511e933bb4b225d30f8cab46d` | `ec68532eb3981957bdb46f85d3a89e2b14929ac5` | 2026-06-18 | cleaned docstrings |
| `2a3060037a102361bdf1afeee954568d8b331399` | `51a569401dcac7e02f9e891d8c4caddff59aa600` | 2026-06-18 | added readme |
| `c30f2ebc66b8029afa1926690ac67c60d9ea7cb6` | `7831cdbae6a4e1bbb6d6709860e0beab3a0fa1cb` | 2026-06-18 | added readme |
| `d895ddda8bd4f87c9d3af48f6d0c8f825360fa1f` | `17c76c8931d5c9461b37882e3fe745abc81ef4f2` | 2026-06-18 | Stop tracking M_ed.npy and ignore computed outputs |
| `948a865bd88657133401dc7f2d10bfed64de351c` | `eff2c95e49afb4ba961e476705167e542aaf0d6b` | 2026-06-19 | added test for pristine supercell |
| `9fbf3691d97b245240613ad04133c083ecb3c9d1` | `7d75ede1c3a7aa6643aac47c88791cfb26ff134a` | 2026-06-19 | Add exact zero-padding densification of M^L |
| `d7657f84fe9a72779d486046a5682f0cfdc41f04` | `429fc6d017107b96084fe9997fb5ce80bd6d95ad` | 2026-06-19 | Clean up scripts: centralize paths, dedup band helpers, drop run_mpi |
| `a706af57e1c144d84100fad2beb82818d1816eb1` | `f8a5f1a70adeec3b58ebd24d3b5871dfbc455ce0` | 2026-06-19 | signed convention for lattice vectors |
| `f7c465716898a26bba9911f77c4b35614f81473e` | `cf05e92fae7c16ece908e6e138ad294c329f81d6` | 2026-06-19 | Consolidate M = M^L + M^NL into a single two-stage driver |
| `beafada4ea331dcabf17c2da7c4f5840e7f14199` | `fe49fd602c2d4026e5e713514fab39df2cafab90` | 2026-06-19 | Split M driver into ml/nl/combine stages |
| `ed2690fdb5cdf7091e0ce6cb23ba6bd45438049f` | `f317ce86186c6ad20bb01c79d7a3fdc2db547a39` | 2026-06-19 | Add per-stage timing to compute_M.py |
| `a101ac4cb6668519a163e1e0b12e678f564685fe` | `7211862d545435143cf22ac94bc42a71ad49ebe7` | 2026-08-29 | Fix heap corruption in compute_M_NL_mpi via full-signature einsum |
| `e61d18b7e021a629e009e27d8e82abeeba3b30ec` | `078cd1670f7cfeeea6ba4dedf6a993f766436a56` | 2026-08-29 | Add single-defect T-matrix/DOS post-processing; drop unused netCDF4 import |
| `c38e47a9160fdf05f476dcf9c05836d0d6de8180` | `7ed44916f786fc5e83191405dd172c84af4bad4b` | 2026-08-31 | Add single-defect scattering rate/lifetime; make compute_M.py QE-only |
| `05fa1e83109fb0269e98a194f3d20d818306cb12` | `024713db7c8617e896df1f263845c63c4b2cb403` | 2026-08-31 | Add supercell-size convergence synthesis of single-defect observables |
| `c8e2613c114e6503c97ad4108a33df93a08a1476` | `143cbca5cd6fcd645dafd01bbc8410628fdde7e9` | 2026-09-02 | Single-defect T-matrix: supercell Bloch norm, structural k-pairing, local scheme + Wannie… |
| `c682468cd027a56efd9d3d6892e3821b77ff26e2` | `071dffa887d3e63dc149f2986e89b51dc019b312` | 2026-09-03 | Wannierize 5x5/7x7/8x8 unit cells; fast grid-cached local rate; Level-1 driver on dense g… |
| `2f26d96ebfd3854e3a1a84159e08e9fce514fd83` | `424cfff3a74f5d4698c2d807e3a68f3bb4b0069d` | 2026-09-03 | Local t-matrix: recenter defect to R0=0, feed intensive M_raw; REAL golden test PASSES (1… |
| `f3d81339368fedb3344f049ecba5bbeaa7b076aa` | `bd1636d82439a5cf20240d8709b350d39dcc8a97` | 2026-09-04 | Dense M by exact zero-padding (option ii-b): 3-stage driver, --dense mode, boundary guard… |
| `5289bc034890232cde74a325e31114891f140620` | `f6cf5956734ea3530a408ca93d60a042b925dbec` | 2026-09-04 | Add dense Level-1 submit (--dense, per-size dense wannier manifest, R_cut 0-3) |
| `c991ccf86890d592ec9254c28f8f40d63204f2e3` | `496792d4f63b51e4940f882b3d9d2abf99dce7e0` | 2026-09-04 | Dense unit-cell wannierizations 25/27/28/32 (froz = E_D+2.5 eV): pi and pi* frozen, manif… |
| `edb7c24ef4dd3ed4333cf9ea9c18df2df07c4238` | `9b58a024b8099228ee9addc1038fda98ea886e14` | 2026-09-04 | Dense M^L via real-space kernel with node-shared u_nk (exact zero-padding, BLAS): replace… |
| `e2c0380faee6207b84a89ba6c70dbb21613be462` | `975313a0acec7b300962c44976a86406efbe7876` | 2026-09-04 | check dense M vs coarse: gauge-closed band subsets + optional dense file arg |
| `20142ecc29a97098a537e5ab2588f0185fe88525` | `13f33bc51297b3da8d30b681b580cf6275c65b76` | 2026-09-04 | dense chain nbnd=20: on-site pz-pz / M^L vs M^NL check, nb20-vs-nb16 SV check, output tag… |
| `4e1e6b4540d463d018f0a8522ad7362b7f5e2890` | `30c30f66bc3a6a81ceb47a9b7047a79f32fea731` | 2026-09-04 | local_green_batch: exact batched g0 (translation invariance + eigenbasis zgemm) for scatt… |
| `bdd876c5cd235ec9cec4b9e0f072c6a22f876bf3` | `ed0b08b8ce677905cd30a92c9afbef85aff4a994` | 2026-09-04 | summarize_level1_maps: R_cut / grid / eta tables from specwd logs |
| `29a6ba238754ad00036bff7ed426e3565fa5a639` | `690234cd22b03fe6088b77fe920ce1066b08ce02` | 2026-09-05 | UNITS FIX + frozen production config |
| `cdfb053505b26eb6989b4e095cdd13058edf16e0` | `efb5556f4ee359fc023b02cb00fadf7504a5c0cc` | 2026-09-05 | Production data (units-fixed, frozen config): Level-1 maps x4, resonance metrics 9x9, Mwr… |
| `d8b26fee43cfa51c6619e3d3edf53a7886cac289` | `f8766059eaaf85b95c15097ca8d020bd0c8542dd` | 2026-09-07 | matrix_io: mandatory units tag ('hartree') on M sidecars, explicit units at load, single … |
| `2c9f0b7d51c80ae5df32079ce34ea4d13b23e1ae` | `9be768f983906e55080fdc43899198b48aa53123` | 2026-09-07 | resonance criteria 9x9 (det/eigenvalue, full-band sum rule, Gamma at c=0.1%) + fig_spectr… |
| `ea9b169d7c18818cc74e1ff34ec23a874f40f859` | `32b180a951234afb59092b04d80e094781859d15` | 2026-09-07 | §4.1.5 : figures de M (carte \|M\| à K, V_ed^L site→frontière, scaling intensif), tableau… |
| `5771f574cf2055f6f30cd1ca94eb7a85388f41cc` | `05bc0732a1103698bc145a6b31e190b978c8fef0` | 2026-09-07 | Style memoire.mplstyle (version de l'auteur), 8 figures régénérées ; fig_M_map : panneaux… |
| `1ce2295ab643e519b1afdfeff84e7c684fea8e4c` | `d757421ca38387b471407e6d09eff91079533863` | 2026-09-07 | Reconstruction KS avec V_p (coarse + dense, 4 tailles, restriction de Fourier si grilles … |
| `268c198bf34f52ba2b0982c8162598bd8c2f7f8e` | `a087a556adc0722a2a89ea11f02bdaad98460521` | 2026-09-07 | Noyaux M^L réels : rééchantillonnage de Fourier de V_ed sur une grille commensurable (7x7… |
| `b1fb837acc73c181f539081a6e8e2ba1d45cea9f` | `7250f24a6e4645d413dee03406d204e03707711e` | 2026-09-07 | Tailles 6x6 et 12x12 (famille N = 3m) : config (p=4/2, D=24), cartes de tailles dans les … |
| `fe7abae59fb343db7cf2da8d74c9cd90c5ee5a74` | `a30a637fcd0d2c606c4f2c8f1624c3d779aacbf0` | 2026-09-07 | Six tailles dans les scripts d'analyse/figures (familles N=3m marquées, on-site pz du sou… |
| `b18fbbc62c99497959d78e09095d8837c2a60020` | `a934f731df2ba80e849bd11561ce4ce0fee187c8` | 2026-09-07 | Équivalence A/B documentée (CLAUDE.md, tag vacancy_sublattice dans les manifests de M) ; … |
| `cde8aefbb408730f5f938bbf6f5623f598770179` | `d9716ee075a00cd70a9f075da137baf8c1fe9391` | 2026-09-07 | Table d'échantillonnage : 10x10 et 11x11 (V_ed frontière et radial) |
| `7071541d424a4d74db2271bae7e97a1764849299` | `6e46ffb8e195f2ef5ec32bd70c3d023c073380e9` | 2026-09-07 | Reconstruction KS 6x6/12x12 (coarse + dense 24x24) ; table d'échantillonnage complète |
| `12be9e454ade7109abbd15ceafe569e7c575246d` | `b08b861802ab4b6fba2f18d771d8c2618dba0028` | 2026-09-07 | check k-coïncidents : troncature au nombre de bandes commun ; localité 6x6 dense |
| `f4c4cb32c09a851cf1c364f1e9f2bb9ea14d38b8` | `032b06246c13f2b709ddc960a328691d5353de0c` | 2026-09-07 | 7x7 dense recalculé (noyau corrigé) : quantification de l'ancienne erreur, check k-coïnci… |
| `f4676ede825cdb5ae86fe0250f65a5dadaee4fef` | `2c87c461d4778b728f46ae834c90a7406315ba5a` | 2026-09-07 | 12x12 dense : post-traitement (tag, k-coïncidents, localité) ; correctif N pour les taill… |
| `3bd8a316a5de71d7dcec80c3ad184d1368c3763c` | `9a1ae2c7d58b06db86af982bfbdec492edc0d104` | 2026-09-07 | Cartes niveau 1 : 6x6, 12x12 (120/240 × 0.02/0.01) et 7x7 recalculé ; niveau 2 par famill… |
| `9886f80cdf23f925a35946489758da8a21bd32f8` | `69f22ff37a94728b016721aa2eecadab649190df` | 2026-09-07 | check_M_dense_vs_coarse : troncature au nombre de bandes commun |
| `1c464a478f0d0ea9a1b359067f2995ea047f5ae7` | `1f4a4a0dcaaf6471fd459fff5c919a8b5465ab6f` | 2026-09-08 | Résonance 6x6 et 12x12 (métriques), critères 6x6 ; tests d'or 6x6 (6.1e-9) et 12x12 (3.7e… |
| `dd0c626d09d118d9d638b845678ea4f60635c8e4` | `32b18686f5994bef26babadda80ecc6f886be3df` | 2026-09-08 | rapport <\|M^NL\|>/<\|M^L\|> sous-espace pi complet pour 6x6/9x9/12x12 |
| `add458689c492bc14e983ed37375e85942419cf5` | `bdfb116130d5af9f6f55c048282d38c2d8eee9c9` | 2026-09-08 | Figures finales du mémoire (make_figures_memoire.py) : fig_convergence, fig_locality_fina… |
| `f11ff2c74fd93fad3bc6200157f4943eec7a16f3` | `46e5e90905cdaa8b77c2163d8b0e2f488b23c04f` | 2026-09-09 | Convergence de la troncature R_cut au niveau de M (9x9, 12x12, grille fine 60x60) -> m_rc… |
| `14e07db1bd7e232c36f89409b5d9c5a2345d3ac3` | `d617f4bef43e5962182b059f53c6d73a637f22bc` | 2026-09-09 | Ratio Frobenius L/NL sur le sous-espace pi/pi* (toutes paires k',k) pour six tailles (lnl… |
| `a9855af8d119aea897249092441e3d7aa3ec39af` | `ff39c7f3bd48ae0d486274d2aeddec9e21e1b57d` | 2026-09-09 | Complément R_cut 9x9 : Re/Im Sigma_nk sur couche (t-matrice locale, grille 240, eta 0.02)… |
| `cb94d4b96651bcf33c70c8fe15ab57e8d0b94a99` | `e8fe798586b19182c10f2afaf43ce182cc8eba16` | 2026-09-09 | m_rcut_convergence.csv : R_cut 4, 5, 6 pour 9x9 et 12x12 (49, 81, 113 mailles) |
| `d2d0bacf3fa085e88b269dc4362c87157efb09b5` | `9e0d1a96d603f669ae4107d356ee9eadc9d668e0` | 2026-09-09 | rcut_resigma.py : option --npe (test de convergence du pas d'énergie eta/n_e) |
| `d49f5bb200767c2cc55112eef3f38fe565aa4712` | `2be8d79e6aa20ac5b0f7a3a6a8f9af6b1237efcc` | 2026-09-10 | EPW P1 : extracteur prtgkk (reproduit les fichiers 16k-8q octet pour octet), lecture dire… |
| `305f4ef1f4ae0995573d731000f839450cced55b` | `3458d806b8980a04663aff70fa0fd967f88b023e` | 2026-09-10 | Chapitre 5 : pointeur CLAUDE.md vers graphene/qe/epw/NOTES_EPW.md ; validation P1 24k-24q… |
| `ba151b375c265d69486e6f53bea4c41475b7d4f0` | `95ee83c6046846e8e92fca0d24e8de9ecc647766` | 2026-09-10 | Chapitre 5 P2 : elecselfen EPW (module electron_phonon/selfen.py, post-traitement epw_sel… |
| `bc83f91e3ec92ab784f36d4ef68f5186cee61d80` | `26a9956a95b76751735fd5e75482619ebf091b81` | 2026-09-10 | Chapitre 5 P2–P5 : production elecselfen 240², degaussw 0.02, T = 300 et 10 K (npz tagués… |
| `133dbcc828483ee18a731319dd858999ec65bd30` | `97c7b9862e7cb6587f1155baf228d650f8b76e79` | 2026-09-10 | fig_epw_gamma : titre du panneau (b) raccourci (coupé à droite) |
| `0afbd840be7ebb747858de07a279c35f472895db` | `07ffb57d1c1f5033cd56aee8823ef411a86d8169` | 2026-09-11 | Chapitre 5 P6 : largeurs de phonons EPW phonselfen (module electron_phonon/phself.py, epw… |
| `b0b91629990f4cc6e7d5b0c19030e022d1e64906` | `10af03a7337b84d42ab70d18b9b44a4924c1556a` | 2026-09-11 | Chapitre 5 P7 : extraction de <D²_Γ>, <D²_K> (epw_d2_extract.py, deux routes) — vertex à … |
| `38ba454c25ab3a521e15c9126af333d25eef27d6` | `0e06daf63d8256f12c0a11499869f3a0455fdc64` | 2026-09-11 | Chapitre 5 P9 : fig_epw_decay (décroissance de H, D, g en représentation de Wannier, chaî… |
| `e28166e831656310efb5edac22a0fb0279785749` | `d3732a649a0ce9e88031ce4757d1e515c44f670c` | 2026-09-14 | Chapitre 5 P10 : fig_epw_decay — étiquettes en convention du mémoire (C au lieu de D en e… |
| `7ba3bf4944e0445c6ac91662b6c65c69730da78a` | `04c19f46201b8d5f003980371ffbe222a5cb2364` | 2026-09-14 | fig_epw_decay : exposant (αμ) entre parenthèses sur g, comme g^{(ν)} dans le mémoire |
| `e3e2f6f963854bbf4dae6f2e5f909f8ae7bf034d` | `b8f6390320629f0a47f068b360dd5089897e3470` | 2026-09-14 | Chapitre 5 P11 : test prtgkk sur l'anneau résonant (epw_ring_check.py) — gamma___ d'EPW =… |
| `a55d7c3c971e64aa015759ca62f3a4911d4bcb9d` | `c3c00aa20c6207cacdfe570b40c9ed4910b26dcd` | 2026-09-14 | fig_epw_decay : suppression de la ligne et de l'étiquette « apothème WS » |
| `c44045cb01f80b15e81ad8b43106476276645f7f` | `44fba52b6478ced0e33dedf63835930743b176d9` | 2026-09-15 | Chapitre 4 : suppression de toute mention de la famille N = 3m dans les figures (titres, … |
| `48cff0cc6a4a03c67ca97e94a784d23b6de173db` | `2987b0de0c452a7820335b6cde921799c45944db` | 2026-09-15 | Chapitre 4 : fig_Ved sans les repères verticaux (1er voisin, R_cut = 3 mailles) ; fig_loc… |
| `65b2f40455270796e32bd767b204ab351fc30bb9` | `034a4879f08edb84d0826731371b15fdcf803aca` | 2026-09-15 | fig_Ved : titre du panneau (c) « Moyenne du potentiel local dans le plan du graphène », o… |
| `4ebec9f995928d283278cfd74de783f169f1e2ca` | `54bcbf461f7ec4bc78db1f612154fea5d02f3fb1` | 2026-09-15 | fig_locality_final : légende réduite aux types de grilles (valeurs pz–pz retirées) |
| `a839fca638c3a57a5c1f16167274e0c67b4b6b97` | `2d6aec1ac99673abba4d5e26306a58e4fd1845d0` | 2026-09-15 | fig_M_scaling_final : légende « Grille élargie » / « Grille grossière », ordonnée limitée… |
| `ab98f6b85b3c3a1a000e102af92ef67b7fbcbcda` | `3f2ac7a63bb418b4203154aab91666a31bc9a2ef` | 2026-09-15 | fig_M_scaling_final : titre « Vérification de la convention intensive » |
| `2a0f084ef20d8cb1fa95a675bc5b86b46f26d6ac` | `e2374ed75f69a9258039faebd0ff28b54ba4978a` | 2026-09-16 | fig_locality_final : R_0 → 0 (\|R\|, M_ij(R,0)) ; fig_M_map_final : titres (a) π, (b) π*,… |
| `ec0fc43880bb99ca331a08640728fb5f8fed5272` | `b1b91c030cc99b03d60d01a98bf04300199eaa0b` | 2026-09-18 | P14 : message du test d'or — seuil réel (rel < 1e-8) et écart atteint affichés (test_loca… |
| `26b20a98185bddc4d874f628d2cd459a33406db0` | `d1b78776afd5dee44c38617b31b2c5147598c22b` | 2026-09-18 | P13 : balayage de la grille interne N_k^int (150/300/450/600, 9x9, R_cut 3) — option --nk… |
| `327355845632e1044a04bf07c7625e1897de51c1` | `20c42db09677e827af78012c0b36aaa001cf246f` | 2026-09-18 | Created README.md |
| `cd1b697c7591e4428f31c0325825ca8022df4884` | `464da6d3b5268c65dd153bebe4498642023bb98f` | 2026-09-18 | Updated README.md |
| `02764365b60c5678c68364aa68c54449cef63fe9` | `cf6acb8fc33abe0c7af0f578d32822af57b1d988` | 2026-09-18 | Updated README.md |
| `b77386df85d5b5780bb3c19b706584934d504f13` | `5ea739110992b520197148c8ec258129ab92dda9` | 2026-09-18 | Created LICENSE |
| `a7363da2f5cbdf87d00ccd6dee11d235e7dfde4d` | `84f66cdb1f28504b1c65869ebe5b644f42880a96` | 2026-09-21 | P18 : chaîne EPW complète refaite à degauss 0.02 Ry (24k-24q_mv0.02) — phonons/bands/nscf… |
| `722125786b4890cfca6ba61ccfe7f02fb765e8f5` | `c95035cb8b6bcf3920f15704e9af5b0a856926e6` | 2026-09-21 | P18 compléments : selfen 120² dg 0.01/0.05 sur 24k-24q_mv0.02 (fig_epw_gamma (a) à quatre… |
| `b5855d731abcd109acee49a346afd4a0d810d8e0` | `9a6ca0f946672f62b3f54424d750576e53799fbd` | 2026-09-21 | P18 complément : contrôle \|g\| sans filtre fsthick — toutes les sommes G > 20 meV, q = M… |
| `73917d151dda47a12e706a9ab620994b86244993` | `de95d96fb0c739b4eaf08e4368322bfe23e0297a` | 2026-09-21 | Rapatriement des notes EPW dans le dépôt : graphene/qe/epw/NOTES_EPW.md → NOTES_EPW.md, g… |
| `c867791b28c9114c7d743b2960b31fd045fec034` | `e2c1301359681a0751587922d9e92f333b2462a0` | 2026-09-21 | Figures EPW : étiquettes mises en cohérence avec le texte — validation (légendes DFT/EPW,… |
| `9c54d994c6fe09879a9036be307eb55a51a41102` | `2cbb6918f06c3f7091f19d1220fac7875eadba61` | 2026-09-21 | Palette harmonisée des figures : bleu marine #000080 (darkblue des hyperliens) comme coul… |
| `800150ee5ecd2024e269837cf8ee0c32ceb7e70c` | `fd58021155939185a2c9cae254bda87629e51988` | 2026-09-21 | fig_epw_phonons (préparation) : bloc 7 de make_figures_epw.py — (a) dispersion matdyn Γ–K… |
| `229d029d2932794332c988fae555684541eaf226` | `003734efaa1493dfa02eec1586b5d59d28d86b31` | 2026-09-21 | fig_epw_phonons : dispersion matdyn Γ–K–M–Γ + DOS matdyn 500² (ΔE 1 cm⁻¹, σ 3 cm⁻¹) sur l… |
| `78ed94af3a898eee5d094c637bae3b3f7b905cd2` | `24a44de1d38c099a5344200fe7029b240eb89468` | 2026-09-22 | R1/R1b : relaxation 9x9 lacune (nspin1/2, Γ) et contrôles 3×3×1, inputs, scripts et rappo… |
| `44ceaa45c7f387227dbb74ef6e84dae3ea90e132` | `7a225ca78d9803dda3e05a03fdb2813e6494a5e1` | 2026-09-23 | cluster: état des lieux et plan de ménage |
| `11d169558ae12fa52d47145f7a57fd2b102ca9e4` | `856d13ea6c51a2ee69668687a6315baf7295725b` | 2026-09-23 | R1: relaxation de la lacune 9×9, Γ, nspin 1/2, contrôles 3×3 |
| `9fbccd6abec164d228bd15fb48ce281484fd1e4e` | `0a248ca595681ee4b5f87dd06d43b91db662f42d` | 2026-09-23 | CLAUDE.md : campagnes — répertoire de travail (graphene/qe/…) vs copie versionnée (articl… |
| `e19538c6005eb4eea942ba4b718b19df17244c64` | `91afafa1be6dcb01cbbcc342110779364e501480` | 2026-09-23 | cluster: étages 1, 3, 2, 6, 4 et 7 exécutés (1 927 Gio), journal dans CLEANUP.md |
| `6eb0a4e4a2a57a5355903c5256b12cc47db4f099` | `bb8b0e39222d5df6cfe5cedec0f0bb54425b9ede` | 2026-09-25 | gitignore |
| `830d5544a303b2a1211ddf579294f3ad4072d3b8` | `7a21450d88f82ccbee3ce5295263c3ec88f9f21a` | 2026-09-25 | nouveau module optics |
| `ba43b5ec774d7368dae86a9247805569a37d4848` | `669f99504caee5630d7d1bc8c086364916072073` | 2026-09-24 | R2 : série en taille 5x5–12x12 de la lacune relaxée (inputs, scripts, rapport) |
| `0beb3b854c2af7a1ebb9de266cfd580c3aad89f6` | `d21965f75674ce3f00b02615e5357c5c326a5e61` | 2026-09-24 | cleanup: journal du 2026-09-24 (miroir R2, piège setgid rsync) |
| `eaba955a641f261accb17a13e00fcf901177aa19` | `c5f6bc78a85be065c19a9f0ac3419516eff0f537` | 2026-09-24 | cleanup: étage 5 exécuté (scratch qe_conv + wfc en vrac, 60,2 Gio) — plan P13 complet |
| `7e863586e3ea0abe2d82429a0c69eb7592d559f8` | `cb2b845602088502946e0897ce32c4b22273f91f` | 2026-09-24 | Renommage du dossier local en graphene-raman : plus de chemin de dépôt codé en dur |
| `d92014003bb319aacdf3c452d70e1886532d0b2c` | `2c82101fb85b3032b112862fde50f5325171e628` | 2026-09-24 | cleanup |
| `75c8656e0dfbcf92adcb05ce44daed06afed6356` | `112617d4552a0294f4e5cc53734730187d18b1d0` | 2026-09-25 | R4 terminé |
| `f4b7ec3a4825627dc25e170af71c98064b33aecb` | `a5afdd42294d7826b023b3df89dc64a351d0e148` | 2026-09-25 | EM1 : éléments de position r(R) de la wannierisation de référence 27×27 (restart = plot) |
| `9fa76f057f03f6289e59f60da76840284e8ebd4f` | `fbf6bc5bd3898c0dac3525f4a35cf924616d7138` | 2026-09-25 | R5 terminé R6 commencé |
| `94f90b2155182de2788be148d794777c9f3c666d` | `18430a7d0549fd94f245f0e3dc2227ebe2c7be9a` | 2026-09-25 | Merge branch 'main' of https://github.com/Gregb26/ab-initio-electron-defect-interaction |
| `ee051ec52e59a6056e1c467163e0fbae473d6c80` | `99d64da6061862b32a0aa50fa753c18901be773a` | 2026-09-25 | correction de la normalisation du noyau local |
| `b9f36b0ac37eea083486aefed3663b0e9c571980` | `7187266195028fa6f42a034c01a339be0a8ebdb5` | 2026-09-25 | R6-R7 |
| `b3ad54b6824ef475985585a1df01f3987e47d392` | `54ac2fa377bc24d060a5cea5999739c0cde6b5c6` | 2026-09-25 | fourier transform tests |
| `09280bbbbcb0f9b3927bde39e4135aa9c8616780` | `a223687194f39251166037a7ad0fe97bb35517a5` | 2026-09-25 | fix normalisation v2 |
| `b500a346a9746a7f459a889b4200ce719431c93b` | `7faffbe62c7268658f3d1ac891968f4605aaf5db` | 2026-09-25 | f4 terminé |
| `1fd2cca24d872a5ad278c91e03658a3dd53f402c` | `841043fedbc4fabe0f654ba8e437a5677427a822` | 2026-09-26 | R6-R7 checkpoint |
| `ba35525e076108197d9eda5cf34350cb227863da` | `a8c16cf59573a5f5334474506f8909a50e265f3a` | 2026-09-26 | gitignore and diff |
| `2e49ffaefa9f75cd910a4cbf758f45b588077e71` | `edd6ee68d1489b99a9aaaa06316611820754cc2b` | 2026-09-26 | r6 stage 3 |
| `fc12a29fe7ea0b73a2aa1179439d2e23f81538c9` | `887b3059437175fe801e86bbc01c4dd1dd074d3b` | 2026-09-27 | R6 terminé |
| `3f364d7fd76b3ae75c7c07855e56cf8679794b8b` | `4fd6575f559a0d9338cf8b3c2f6cd819ced8955e` | 2026-09-27 | post r6 |
| `e4122abb57cbbd42e99c719bdb6422c0507a9cbb` | `d90c929f77f20659671213e9d19394e0ccac002a` | 2026-09-27 | R6 clos |
| `e176c48a4e36a4ab227af590708f7e8b68f468ac` | `aae4f6de529cabba3bfac786f7669e5af5605bc6` | 2026-09-27 | tests velocity operator |
| `00fecbb65dfe8591cd1cc5913763f7320e545c1e` | `4c37abf7609294c19267a05db8fd34d9770b17b0` | 2026-09-27 | velocity operator tests |
| `ac6d787492b636aff184d45fc9a0c6ac74898990` | `3026de157789ff584cb700fc58f97e8414c77e70` | 2026-09-28 | velocity operator f1-f4 and testing completed |
| `5589c99c2c43a4750821cf6fbf4b8d9b038520ff` | `780ad41a096da258e4899e54700fc9790737ce82` | 2026-09-28 | electron-photon series EM moved under mémoire/EM |
| `67e3147bf81da1122ed281dda22b2c973d1a4d65` | `b17ca6b0c1f572e6654d53df78e284db79b0dc0a` | 2026-09-28 | added plan for EM series |
| `7b72280bc6c4306e30cb1cf55d46cbf774dad091` | `9ecd9bbfde6d6ec72ef4b3f74845c14763c36603` | 2026-09-28 | f6 and tests done |
| `b257e7b0418eb8f531af7aed08cbe41fa10ddad4` | `f2518d7015365d5e846843272f013c942196ea46` | 2026-09-28 | R9 checkpoint |
| `5a4bc9478f6d01832b8a595bd1abc260cae801bd` | `cb7241d40ee893ba0743c332b481953aed681fec` | 2026-09-28 | R9 checkpoint |
| `3f32ee3951ea6cec19a61f25501a2b717c403efb` | `e4da0d842fc2fa7892ef0198ae9f0acb9f26ffb6` | 2026-09-28 | completed and tested f6 ring |
| `4eee646a61cae8f4c0178c6579e8dbe5c1ad2d08` | `a25a9199ae6799b041c224e163092ec6fe5d8bde` | 2026-09-28 | finished and tested Kubo, refactored into electron_photon module |
| `fa143e4ac61606846cd799199b3c14e3afcc9177` | `4da063670d279b38da42455dd43e04560c4cf8e0` | 2026-09-28 | finished M0 |
| `603afb9af529911acd0b99cbce9c42fcb45b6994` | `fe7b20b0d9709cf5913ac5594061650be86970aa` | 2026-09-28 | added function to create WannierTB data class from Wannier90 outputs |
| `7dd87319f38c9c9c27c52720368555352a3bd564` | `c8eea59ebecf0605614d92b2140d702deba3aaf0` | 2026-09-28 | added function to keep only Wannier centers in position operator |
| `e09846c8e97fc9ba94dfe19a4eeae37230c89d58` | `81dabc4fa580197290dec9e93a3e3855603fa4f6` | 2026-09-28 | added and tested a function to test hermicity by block |
| `ab2d88427fa68b25783c8b86ff39d4754b6212a9` | `18838614c626a27a0f6cf0dc09e10c05e20313a7` | 2026-09-28 | R9 finished |
| `e4a6c7a3833578600fd86cd7a1fcebc155881e08` | `d1d1ab71b0f8d10693cfffab2507fd3c037c0014` | 2026-09-28 | added and tested a function that checks the symmetries of Wannier functions in graphene. … |
| `d1161f40c1cd4f3f157b87178c1947c9ad92b303` | `cff65cdd446328d0b530bca14dc1f64e0517feb9` | 2026-09-28 | R9 done |
| `3156ad1e757014470cb5ac127757853453a91c56` | `de5eaf597f7a8503f8442d18800d9af81aab1d8c` | 2026-09-28 | M2 done |
| `c178ddf5bcd070356d8b07e74e3996e8e1a3eced` | `27d067df00bfa76a96c9782e3b34d86c5b779479` | 2026-09-29 | finished coding M3: symmetries, nodes and resonant rings |
| `82b46fd1b80a63d897ff1457133b317d640e98c1` | `52bbd57381070500baca2e95e4baf6d4f6678c1d` | 2026-09-29 | DOS and spectral functions from R8 (coarse) |
| `b1d2da3c752902e42bd820db43df64a278aa03ea` | `ea8c039904d2385255f8b35967188929e1d9a1e2` | 2026-09-29 | smooth DOS and spectral functions, yay! |
| `49ae2e316bacb55bb1c68ace9bf72d9f39e1f961` | `89491ac3c2856207edc6f54261e617b2a90e9b83` | 2026-09-29 | added and tested a function to find frozen window limit |

## Commits de la branche supprimée `cleanup-M-driver` (dans le bundle seulement)

| numéro | date | sujet |
|---|---|---|
| `53148f6172386a6db4b4858063694fb1bbcdda3c` | 2026-08-29 | Add single-defect T-matrix/DOS post-processing; drop unused netCDF4 import |
| `5c3d130680c4acbc2d06d728b4209a11d79d8d65` | 2026-08-29 | Fix heap corruption in compute_M_NL_mpi via full-signature einsum |
| `a17109876812aae48e1a1a6f950d81edb7f8a005` | 2026-06-19 | Add per-stage timing to compute_M.py |
| `b6cfcd782e8dec25f2754668738946c8edef0fa3` | 2026-06-19 | Split M driver into ml/nl/combine stages |
| `070b235d4525970c31533d98cea96ca3e8bb66ac` | 2026-06-19 | Consolidate M = M^L + M^NL into a single two-stage driver |
