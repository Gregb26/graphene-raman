## R — Préparation des géométries relaxées (J3, `r`, 6 s ; rien codé)

Lecture du code (fichiers et lignes de la chaîne actuelle) :

- **Partie non locale** (`defects/non_local.py`, `compute_M_NL`, l. 109–201) : les positions de **tous** les atomes sont lues
  séparément dans `sc_p` (`tau_s_p`) et `sc_d` (`tau_s_d`), les phases e^{∓iK·τ} sont calculées atome par atome et M^NL = M_d − M_p :
  une géométrie relaxée est traitée par construction (il suffit de passer le `.save` relaxé comme `sc_d`). Aucun indice de spin ;
  les projecteurs KB (PP_DIJ, β_ℓ) ne dépendent pas du spin : **ΔV^NL est identique pour ↑ et ↓** à géométrie donnée.
- **Partie locale** (`defects/local_R.py` l. 15–97, `local_G.py`) : ΔV = V_d − V_p lus par `qe_io.get_pot` (un fichier filplot,
  aucune notion de spin, l. 393–453) ; positions non utilisées. nspin = 2 exige deux potentiels pp.x (spin_component 1 et 2, fournis
  par J0) et **deux passages de la partie locale seulement** ; la somme M^L_σ + M^NL par spin reste à assembler par un pilote.
- **Lecteur de wfc** (`qe_io._read_all_wfc`, l. 342–391) : `npol = 1` exigé ; les fichiers `wfcup1`/`wfcdw1` d'un run lsda (même
  attribut `ik`) ne sont pas distingués ; gamma_only non géré. Sans effet sur la chaîne M (seules les wfc de maille sont lues).
- Manquant : pilote à deux passages locaux + un passage non local avec `sc_d` relaxé, somme par spin, sidecar `matrix_io` portant le
  spin, et deux potentiels pp.x de la parfaite si l'on voulait aussi une référence spin-polarisée.

Décroissance de ΔV au bord pour R1 (potentiel de la géométrie relaxée moins V_p de la 9×9 parfaite non relaxée, même grille) :

%%TABLES_R%%

Sur les plans frontière, max |ΔV| passe de 34 meV (non relaxée) à 2,59 (nspin1) et 3,67 eV (nspin2) : les plans coupent des cœurs
atomiques déplacés de ~1e-3 Å entre les deux géométries. Le bord de la ligne a₁ vaut +89 / +92 / +89 meV (nspin1, ↑, ↓) contre
−24 meV (non relaxée) ; la ligne passe de −4,3 eV (t = 0,125) à +0,22 eV (t = 0,375) pour nspin1 contre −0,21 / −0,04 eV pour la non
relaxée. V_↑ − V_↓ : max 10,13 eV (près des cœurs), moyenne 3D −2,64 meV, au site de la lacune +0,704 eV, plans frontière 57 meV.
