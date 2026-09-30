# Spin-polarized electron scattering by the magnetic monovacancy in graphene: DFT magnetism, spin-dependent scattering and transport, and whether a spin-resolved ab initio T-matrix exists

Scope note: searches done 2026-09-29 (web search + arXiv/APS/publisher pages). "Verified" means I saw authors, journal, volume, article number and year on an arXiv, publisher, or institutional repository page during this session. Items marked UNVERIFIED were not confirmed from a source. I could not open Google Scholar "cited by" lists, so the citing-list sweep of Yazyev & Helm 2007, Kochan et al. 2014 and Kaasbjerg 2020 was done by keyword search, not by walking the full lists.

## Q1. DFT magnetism of the graphene monovacancy: moments, σ vs π, Jahn–Teller, supercell/functional dependence, the π-moment controversy, Kondo physics, and experiments

### Takeaway
All spin-polarized DFT studies agree that the bare monovacancy is magnetic. The σ dangling bond on the one unreconstructed carbon atom (the other two neighbours bond after a Jahn–Teller distortion, giving a 5-9 geometry) carries about 1 μB. The quasi-localized π (Vπ) state adds a partial, strongly method-dependent contribution, so reported totals run from about 1.04 to 2.0 μB. That π contribution is the controversial part: it depends on vacancy concentration (supercell size), smearing and k-sampling, the functional (hybrids give an integer 2 μB), and out-of-plane distortion. At experimentally relevant low concentrations it may vanish. Experiments see spin-½ paramagnetism with roughly equal σ and itinerant (π) contributions, an STM spin-split Vπ peak, and gate-tunable Kondo screening. The Kondo interpretation of the transport data is disputed.

### Cited Findings
**DFT: magnitude of the moment and its σ/π origin**
- Yazyev & Helm, "Defect-induced magnetism in graphene", Phys. Rev. B 75, 125408 (2007), arXiv:cond-mat/0610638 (verified). First-principles study of H chemisorption and vacancy defects, finding itinerant magnetism from defect-induced extended states. The moment is 1 μB per H defect and **1.12–1.53 μB per vacancy depending on defect concentration**. The coupling between moments is FM or AFM depending on whether the defects sit on the same sublattice or on different ones. — [arXiv](https://arxiv.org/pdf/cond-mat/0610638); [EPFL Infoscience](https://infoscience.epfl.ch/record/101414?ln=en)
- Lehtinen, Foster, Ma, Krasheninnikov, Nieminen, "Irradiation-induced magnetism in graphite: a density functional study", Phys. Rev. Lett. 93, 187202 (2004) (verified; authors and year from the APS abstract page). Spin-polarized DFT finds both the vacancy and the vacancy–H complex magnetic. **H adsorption on one vacancy dangling bond doubles the moment of the naked vacancy.** The abstract gives no numerical moment. — [APS abstract](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.187202); [Aaltodoc record](https://aaltodoc.aalto.fi/items/4f1d09dd-edea-43cd-baef-91281fc9a1ce)
- Nanda, Sherafati, Popović, Satpathy, "Electronic structure of the substitutional vacancy in graphene: density-functional and Green's function studies", New J. Phys. 14, 083004 (2012), arXiv:1105.1129 (verified). All-electron spin-polarized LAPW combined with tight-binding and impurity Green's functions. The three sp²σ dangling bonds give localized Vσ mid-gap states, split by the crystal field and a **Jahn–Teller distortion**, and pzπ gives a sharp Vπ resonance. **Hund's coupling aligns the localized spins into S = 1 (2 μB). AFM polarization of the itinerant π band reduces this by about 0.3 μB, giving a net 1.7 μB.** Lippmann–Schwinger analysis gives a ~1/r Vπ tail with two-valley (K, K') interference. — [arXiv PDF](https://arxiv.org/pdf/1105.1129); [search record incl. abstract text](https://api.emergentmind.com/papers/1105.1129)
- Sci. Rep. 2017, "Controlling magnetic transition of monovacancy graphene by shear distortion" (authors UNVERIFIED). States that an isolated vacancy has about 1.5 μB, **with ~1 μB from unsaturated σ states and ~0.5 μB from π electrons**, and that the ground state is FM with 1–1.7 μB depending on vacancy concentration. — [Scientific Reports](https://www.nature.com/articles/s41598-017-01881-3); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5431955/)

**DFT: sensitivity and controversy about the π moment**
- Palacios & Ynduráin, "Critical analysis of vacancy-induced magnetism in monolayer and bilayer graphene", Phys. Rev. B 85, 245443 (2012), DOI 10.1103/PhysRevB.85.245443, arXiv:1203.6485 (verified). Spin-resolved DFT study that finds the **vacancy-induced extended π moments, which interact at long range and could order, vanish at any experimentally relevant vacancy concentration**. — [arXiv](https://arxiv.org/pdf/1203.6485); [UA record](https://observatorio-cientifico.ua.es/documentos/63d5b49ff851ee1ba3ea039d?lang=en)
- Bin Wang & Sokrates T. Pantelides, "Magnetic moment of a single vacancy in graphene and semiconducting nanoribbons", Phys. Rev. B 86, 165438 (2012) (verified). The moment comes from spin-up/down asymmetry within a ~1 eV window that contains a sharp resonance near E_F. **The computed moment stays sensitive to the smearing width even with thousands of k-points per unit cell.** The paper explains the scatter in the literature as a result of these numerical choices. — [APS abstract](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.165438)
- Ana M. Valencia & Marilia J. Caldas, "Vacancy in graphene: insight on magnetic properties from theoretical modeling", Phys. Rev. B 96, 125431 (2017), DOI 10.1103/PhysRevB.96.125431, arXiv:1611.08246 (verified). **DFT results in the literature span μ = 1.04–2.0 μB.** Hybrid functionals (with Fock exchange) are essential for the π states at E_F and give **an integer 2 μB for the isolated vacancy**. Periodic vacancy arrays show "diffuse spin–spin interactions". — [arXiv](https://arxiv.org/abs/1611.08246)
- Padmanabhan & Nanda, "Intertwined lattice deformation and magnetism in monovacancy graphene", Phys. Rev. B 93, 165403 (2016), arXiv:1605.03921 (verified). DFT finds two competing structures. The **planar ground state has a saturated 1.5 μB**. A **metastable non-planar structure has a vanishing moment and costs only ~50 meV**. — [arXiv](https://arxiv.org/pdf/1605.03921)
- Paz, Scopel, Freitas, "On the connection between structural distortion and magnetism in graphene with a single vacancy", Solid State Commun. (2013), DOI 10.1016/j.ssc.2013.05.004, arXiv:1211.2695 (volume and pages UNVERIFIED). The planar, locally distorted structure is the magnetic ground state: two bonds reconstruct and one dangling bond remains. It has σ and π spin contributions. Displacing the adjacent atom out of plane gives **metastable non-magnetic states**. — [arXiv](https://arxiv.org/abs/1211.2695)
- Casartelli, Casolo, Tantardini, Martinazzo, "Spin coupling around a carbon atom vacancy in graphene", Phys. Rev. B 88, 195424 (2013), arXiv:1303.1924 (verified). Uses magnetization-constrained periodic DFT plus multireference PT2 on a cluster. **Two local moments, one π-like and one σ-like, are decoupled from the π band and coupled to each other. The ground state is a triplet** with a planar geometry, in which the apical C atom sits opposite a pentagon. A **bare vacancy is predicted to be spin-1**. Spin-½ can occur with foreign-species binding, ripples, a substrate or doping. — [arXiv](https://arxiv.org/pdf/1303.1924); [UNIMI record](https://air.unimi.it/handle/2434/235190)
- Goswami & Gahlot, "Local moment formation by vacancies in mono-layer graphene", arXiv:1205.1230 (2012; journal UNVERIFIED). Green's function model in which next-nearest-neighbour hopping suppresses the zero mode and Hund's coupling with the σ states polarizes the π state. — [arXiv](https://arxiv.org/abs/1205.1230)

**Kondo physics (theory)**
- Cazalilla, Iucci, Guinea, Castro Neto, "Local moment formation and Kondo effect in defective graphene", arXiv:1207.3135 (2012; journal version UNVERIFIED). A model that includes vacancy reconstruction and non-planarity (strain or temperature). It finds several impurity phases that control the moment magnitude and whether a Kondo effect occurs, depending on U, Hund's J, doping and particle-hole symmetry breaking. — [arXiv](https://arxiv.org/pdf/1207.3135)
- Mitchell & Fritz, "Kondo effect with diverging hybridization: possible realization in graphene with vacancies", Phys. Rev. B 88, 075104 (2013), arXiv:1212.2631 (verified from the arXiv listing and the UCD-hosted PDF). NRG for a host with a diverging DOS, as for the vacancy. It finds a spin-½ (doublet) Kondo phase. — [arXiv](https://arxiv.org/pdf/1212.2631); [UCD PDF](https://www.ucd.ie/nanoelectronics/t4media/Phys_Rev_B_88_075104_(2013).pdf)
- May, Lo, Deltenre, Henke, Mao, Jiang, Li, Andrei, Guo, Anders, "Modeling of the gate-controlled Kondo effect at carbon point defects in graphene", Phys. Rev. B 97, 155419 (2018), DOI 10.1103/PhysRevB.97.155419, arXiv:1803.03196 (verified). An effective **two-orbital (unbound σ + vacancy-induced bound π) single-impurity model**. — [arXiv](https://arxiv.org/pdf/1803.03196)

**Experiments**
- Nair, Sepioni, Tsai, et al., "Spin-half paramagnetism in graphene induced by point defects", Nat. Phys. 8, 199–202 (2012), arXiv:1111.3775 (verified). Magnetization shows that F adatoms and irradiation vacancies **carry spin-½ moments. There is paramagnetism but no ordering down to liquid-He temperatures**, and at most ~1 moment per ~1000 C atoms. — [arXiv](https://www.arxiv.org/pdf/1111.3775); [Loughborough record](https://vufind.lboro.ac.uk/PrimoRecord/cdi_proquest_miscellaneous_1031326963)
- Nair, Tsai, Sepioni, Lehtinen, Keinonen, Krasheninnikov, Castro Neto, Geim, Grigorieva, "Dual origin of defect magnetism in graphene and its reversible switching by molecular doping", Nat. Commun. 4, 2010 (2013), DOI 10.1038/ncomms3010, arXiv:1301.7611 (verified). **Vacancy magnetism has a dual origin with two roughly equal contributions, one itinerant (π) and one from dangling bonds (σ).** Doping switches the itinerant part on and off. — [arXiv](https://arxiv.org/pdf/1301.7611); [Manchester record](https://research.manchester.ac.uk/en/publications/dual-origin-of-defect-magnetism-in-graphene-and-its-reversible-sw/)
- Yu Zhang, Si-Yu Li, Huaqing Huang, Wen-Tian Li, Jia-Bin Qiao, Wen-Xiao Wang, Long-Jing Yin, Ke-Ke Bai, Wenhui Duan, Lin He, "Scanning tunneling microscopy of the π magnetism of a single carbon vacancy in graphene", Phys. Rev. Lett. 117, 166801 (2016), arXiv:1604.06542 (verified). The **Vπ state splits into two spin-polarized DOS peaks separated by several tens of meV**. A magnetic field increases the splitting with an **effective g ≈ 40**. — [arXiv](https://arxiv.org/pdf/1604.06542); [APS DOI link](https://link.aps.org/doi/10.1103/PhysRevLett.117.166801)
- Ugeda, Brihuega, Guinea, Gómez-Rodríguez, "Missing atom as a source of carbon magnetism", Phys. Rev. Lett. 104, 096804 (2010), arXiv:1001.3081 (verified). STM on graphite vacancies shows a sharp resonance at E_F, associated with local moment formation and with reduced mobility. — [arXiv](https://arxiv.org/pdf/1001.3081)
- Jian-Hao Chen, Liang Li, William G. Cullen, Ellen D. Williams, Michael S. Fuhrer, "Tunable Kondo effect in graphene with defects", Nat. Phys. 7, 535–538 (2011), arXiv:1004.3373 (verified). Irradiated graphene shows a Kondo-like resistivity with **T_K = 30–90 K tunable by gate**. — [arXiv](https://arxiv.org/pdf/1004.3373). **Contested** by Jobst & Weber, "Origin of logarithmic resistance correction in graphene", Nat. Phys. 8, 352 (2012), DOI 10.1038/nphys2297 (verified), who attribute the log correction to electron–electron interaction in the diffusive regime ([FAU record](https://cris.fau.de/publications/114266944)). Chen et al. published a reply ([Nature Physics](https://www.nature.com/articles/nphys2306)). See also Jobst, Kisslinger, Weber, Phys. Rev. B 88, 155412 (2013), arXiv:1308.6671 (verified) ([arXiv](https://arxiv.org/pdf/1308.6671)).
- Yuhang Jiang, Po-Wei Lo, Daniel May, Guohong Li, Guang-Yu Guo, Frithjof B. Anders, Takashi Taniguchi, Kenji Watanabe, Jinhai Mao, Eva Y. Andrei, "Inducing Kondo screening of vacancy magnetic moments in graphene with gating and local curvature", Nat. Commun. 9, 2349 (2018), DOI 10.1038/s41467-018-04812-6 (verified). STM signature of Kondo screening, with a quantum phase transition between screened and unscreened vacancy moments that is controlled by gating and curvature. — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6002358)
- Jing Jing Chen, Han Chun Wu, Da Peng Yu, Zhi Min Liao, "Magnetic moments in graphene with vacancies", Nanoscale 6, 8814–8821 (2014), DOI 10.1039/c3nr06892g (verified). Magnetotransport plus DFT (1.2–1.8 μB, sharp Vπ resonance). Transitions between hopping transport and Kondo behaviour, and between giant negative MR and positive MR, as carrier density, T, B and vacancy density change. — [BIT record](https://pure.bit.edu.cn/en/publications/magnetic-moments-in-graphene-with-vacancies/)
- Review: Yu Zhang, Liangguang Jia, Yaoyao Chen, Lin He, Yeliang Wang, "Recent advances of defect-induced spin and valley polarized states in graphene", DOI 10.1088/1674-1056/ac70c4 (Chin. Phys. B, 2022, arXiv:2206.05936; volume UNVERIFIED). An STM review of vacancy, N dopant and H chemisorption states. — [arXiv](https://arxiv.org/abs/2206.05936)

### Inferences
- The robust part of the moment is the σ dangling bond (~1 μB). The π part (0 to ~1 μB) is what matters most for transport scattering, because the Vπ resonance sits at E_F and is what electrons scatter off. It is also the least robust part. A spin-polarized ab initio T-matrix will therefore mainly probe the exchange splitting of the Vπ resonance, and that splitting depends on supercell, smearing and functional.
- Casartelli's spin-1 prediction and Nair's spin-½ measurement conflict. Casartelli reconcile them through environment effects (doping, ripples, substrate). A collinear DFT calculation cannot settle this. It gives a mean-field S_z, not the total spin.
- The STM splitting of "several tens of meV" (Zhang 2016) sets a scale against which to compare the student's ↑/↓ Vπ resonance splitting in the nspin=2 supercell.
- Kondo screening (T_K up to tens of K claimed, and disputed) means a static spin-polarized picture is appropriate only above T_K, or where the moment is unscreened.

### Gaps
- Numerical moment in Lehtinen et al. 2004 (the abstract gives only the doubling ratio). Full text not accessed.
- Authors of the 2017 Sci. Rep. shear-distortion paper, the Paz et al. SSC volume and pages, and the Chin. Phys. B volume were not confirmed.
- No systematic benchmark was found comparing DFT+U, hybrid and GGA for the vacancy at fixed large supercell. Valencia & Caldas is the closest.
- Haase et al. (Kondo near vacancies, PRB 2011) and El-Barbary et al. (vacancy in graphite, PRB 2003) are commonly cited here but were not verified in this session.

## Q2. Spin-dependent scattering and transport with magnetic vacancies or magnetic adatoms: fully ab initio vs DFT-fitted models vs NEGF devices

### Takeaway
Nearly all spin-dependent scattering work in graphene treats **hydrogen (or other) adatoms**, not the vacancy, using **model Hamiltonians fitted to spin-polarized DFT**. Examples are the Kochan–Gmitra–Fabian spin-flip T-matrix with a DFT-fitted exchange J, and mean-field-Hubbard Kubo transport (Soriano, Roche). Vacancies appear only as a "could also be a source" remark. Fully ab initio spin-dependent treatments either are finite-device NEGF calculations of nanoribbons or concern non-magnetic adatoms with spin–orbit coupling (Fedorov et al., relativistic). No work found builds a spin-resolved ab initio T-matrix for the vacancy in bulk graphene.

### Cited Findings
**A. DFT-fitted model Hamiltonian + T-matrix (magnetic resonant scatterers)**
- Kochan, Gmitra, Fabian, "Spin relaxation mechanism in graphene: resonant scattering by magnetic impurities", Phys. Rev. Lett. 112, 116602 (2014), DOI 10.1103/PhysRevLett.112.116602, arXiv:1306.0230 (verified). At resonance, local moments act as **spin hot spots: spin-flip rates are as large as spin-conserving ones when the exchange exceeds the resonance width**. About 1 ppm of moments, smeared by e–h puddles, reproduces the ~100 ps spin lifetimes. — [Regensburg record](https://pred.uni-regensburg.de/id/eprint/10461/); [ar5iv](https://ar5iv.labs.arxiv.org/html/1306.0230)
  - Method details (from ar5iv full text): the paper says "if the local moments sit at resonant scatterers, **such as vacancies** and adatoms, they can act as spin hot spots" and "while the local moments can come from a variety of sources, **we specifically focus on hydrogen adatoms**." The orbital parameters (ε_h = 0.16 eV, T = 7.5 eV) come from an earlier DFT fit by Gmitra et al. The exchange comes from a least-squares fit to supercell spin-polarized DFT: **J_h = −0.82 eV, J1 = 0.69 eV, J2 = −0.18 eV**. The rates are then computed with a **generic J = −0.4 eV**, and the authors say the rates are "hardly influenced by the precise value and the sign of J". It is a T-matrix calculation. — [ar5iv](https://ar5iv.labs.arxiv.org/html/1306.0230)
- Kochan, Irmer, Gmitra, Fabian, "Resonant scattering by magnetic impurities as a model for spin relaxation in bilayer graphene", Phys. Rev. Lett. 115, 196601 (2015), arXiv:1504.03898 (verified). Adatoms on dimer vs non-dimer sites; only dimer sites give narrow resonances at the charge neutrality point. Explains the opposite density dependence of the rate in bilayer compared with monolayer. — [Regensburg record](https://pred.uni-regensburg.de/id/eprint/4474/); [arXiv](https://arxiv.org/pdf/1504.03898)

**B. DFT-derived mean-field Hubbard + Kubo / Landauer (large disordered samples)**
- Soriano, Leconte, Ordejón, Charlier, Palacios, Roche, "Magnetoresistance and magnetic ordering fingerprints in hydrogenated graphene", Phys. Rev. Lett. 107, 016602 (2011), DOI 10.1103/PhysRevLett.107.016602, arXiv:1105.1005 (verified). **Hydrogen, not vacancies.** A self-consistent mean-field Hubbard model derived from first principles gives spin-dependent on-site energies, which feed Kubo–Greenwood transport for paramagnetic, AFM and FM configurations. **MR up to ~7% at ~0.25% H.** — [arXiv](https://arxiv.org/abs/1105.1005)
- Soriano, Van Tuan, Dubois, Gmitra, Cummings, Kochan, Ortmann, Charlier, Fabian, Roche, "Spin transport in hydrogenated graphene", 2D Mater. 2, 022002 (2015), arXiv:1504.01591 (verified). A theory review of spin transport in hydrogenated graphene. — [arXiv](https://arxiv.org/pdf/1504.01591); [Regensburg PDF](https://epub.uni-regensburg.de/41389/1/Soriano_2015_2D_Mater._2_022002.pdf)
- Thomsen, Ervasti, Harju, Pedersen, "Spin relaxation in hydrogenated graphene", Phys. Rev. B 92, 195408 (2015) (verified). Landauer–Büttiker spin transport with magnetic H adatoms. Spin relaxation is not always Markovian, and the inverse spin relaxation length and the sheet resistance scale roughly linearly with concentration. — [DTU CNG PDF](https://cng.dtu.dk/-/media/centre/cng_nanostructured_graphene/research/publications/2015/physrevb-92-195408.pdf)

**C. Fully first-principles spin-dependent scattering, but not for the magnetic vacancy**
- Fedorov, Gradhand, Ostanin, Maznichenko, Ernst, Fabian, Mertig, "Impact of electron-impurity scattering on the spin relaxation time in graphene: a first-principles study", Phys. Rev. Lett. 110, 156602 (2013), DOI 10.1103/PhysRevLett.110.156602 (verified). **Relativistic ab initio** momentum and spin relaxation times for **non-magnetic C and Si adatoms**, where spin relaxation comes from SOC. The specific method (a KKR Green function / T-matrix, as is usual for the Mertig group) is UNVERIFIED from the abstract. — [Regensburg record](https://pred.uni-regensburg.de/id/eprint/16825/)
- Wehling, Yuan, Lichtenstein, Geim, Katsnelson, "Resonant scattering by realistic impurities in graphene", Phys. Rev. Lett. 105, 056802 (2010), arXiv:1003.0609 (verified). A first-principles-based theory of resonant (midgap, within ±0.03 eV) impurities such as organic groups and H, explaining sublinear σ(n). Not spin-resolved. — [arXiv](https://arxiv.org/pdf/1003.0609)

**D. NEGF / finite devices (nanoribbons)**
- Topsakal, Aktürk, Sevinçli, Ciraci, "First-principles approach to monitoring the band gap and magnetic state of a graphene nanoribbon via its vacancies", Phys. Rev. B 78, 235435 (2008), arXiv:0808.1468 (verified). Periodic vacancies metallize and magnetize semiconducting ribbons through spin-polarized defect states. Whether spin-resolved transport is computed there is UNVERIFIED. — [arXiv](https://arxiv.org/pdf/0808.1468)
- A SIESTA/TranSIESTA study (arXiv:1703.06504, porphyrin detection on GNRs) reports that, for the **clean single-vacancy ribbon, the current is strongly spin-polarized and dominated by the majority-spin channel**. This comes from a search snippet. Authors and journal were not confirmed (UNVERIFIED details). — [arXiv](https://arxiv.org/pdf/1703.06504)
- Wang & Pantelides 2012 (above) also treat vacancy moments in semiconducting nanoribbons and find a 2 μB ground state. — [APS abstract](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.165438)

### Inferences
- Classification:
  1. **Fully ab initio bulk T-matrix**: Kaasbjerg 2020 (spin-unpolarized; see Q3) and Fedorov 2013 (non-magnetic adatoms, with SOC).
  2. **DFT-fitted TB + T-matrix or Kubo**: Kochan/Fabian and Soriano/Roche (hydrogen). These are spin-resolved, and they capture spin-flip only through an explicit exchange J s·S operator.
  3. **NEGF-DFT of finite ribbons**: spin-polarized, but these give device conductance, not bulk disorder-averaged rates.
- The vacancy is harder to fit with a single-orbital TB model than hydrogen. It has two coupled local moments (σ + π, with Hund's coupling) and a Jahn–Teller distortion. This may be part of why the DFT-fitted community stayed with hydrogen. This is my inference; I found no source that says so.
- A collinear nspin=2 ab initio T-matrix gives **spin-conserving ↑ and ↓ channels only**. The Kochan-type spin-flip "hot spot" physics needs a quantum spin operator or a non-collinear treatment and is outside a collinear pipeline.

### Gaps
- I did not find any spin-polarized Kubo or Landauer transport study for random **vacancies** in bulk graphene based on DFT-derived spin-dependent parameters (the vacancy analogue of Soriano 2011). The keyword searches found only hydrogen versions. This is absence of evidence, not proof.
- The method of Fedorov et al. 2013 and whether the Mertig group later treated vacancies with spin-polarized KKR were not verified.
- The follow-up by Irmer, Kochan, Fabian on adatom positions (PRB 2018) was not verified.

## Q3. Does a spin-resolved ab initio T-matrix (or beyond-Born rate, spin-dependent resistivity, optical conductivity) exist for the graphene vacancy in bulk, disorder-averaged graphene? Closest work, and what would be new

### Takeaway
I found **none**. The closest work is Kaasbjerg's atomistic DFT + T-matrix theory, Phys. Rev. B 101, 045433 (2020). It treats exactly the graphene vacancy in the bulk disorder-averaged setting, going beyond Born and computing DOS, linewidth and conductivity. However, it uses a **spin-independent (spin-diagonal ŝ0) defect potential with no treatment of the vacancy moment**. Other ab initio electron-defect pipelines (Bernardi's group, Wannier-based) are Born-level and non-magnetic. Ab initio T-matrix work on optics (exciton–defect, MoS₂, 2025) is also non-magnetic. A spin-resolved (↑/↓) ab initio T-matrix for the magnetic monovacancy, with spin-resolved linewidths, a two-current resistivity and spin-resolved optical conductivity, appears to be genuinely new, provided it is framed as "no such study found" and includes the sensitivity caveats of Q1.

### Cited Findings
- Kristen Kaasbjerg, "Atomistic T-matrix theory of disordered two-dimensional materials: bound states, spectral properties, quasiparticle scattering, and transport", Phys. Rev. B 101, 045433 (2020), DOI 10.1103/PhysRevB.101.045433, arXiv:1911.00530 (verified). A DFT + T-matrix framework giving the DOS (including bound states), quasiparticle spectrum, spectral linewidth (scattering rate) and conductivity/mobility, with a **strong breakdown of the Born approximation in 2D**. It is demonstrated on **graphene and TMDs with vacancies and substitutional impurities**. — [arXiv abs](https://arxiv.org/abs/1911.00530)
  - Method details (ar5iv full text read through a summarizing fetch; check exact wording against the PDF): GPAW with a DZP LCAO basis in 11×11 supercells with 10 Å vacuum. The **defect potential is the difference between the crystal potentials of the defective and pristine lattices**, which is the same construction as the student's Ṽ_ed. The defect potential is written as **"V^i = V_i(r̂) ⊗ ŝ0, where V_i(r) is the scalar spin-independent defect potential"**, and the paper says "we thus neglect defect-induced changes in the spin-orbit interaction". SOC enters only through the host band structure. **The fetch found no spin-polarized DFT and no discussion of the vacancy magnetic moment.** — [ar5iv](https://ar5iv.labs.arxiv.org/html/1911.00530)
- I-Te Lu, Jin-Jian Zhou, Marco Bernardi, "Efficient ab initio calculations of electron-defect scattering and defect-limited carrier mobility", Phys. Rev. Materials 3, 033804 (2019), arXiv:1901.03449 (verified). Elastic electron–defect matrix elements and relaxation times for **neutral vacancies and interstitials in Si**. — [arXiv](https://arxiv.org/abs/1901.03449); [Caltech record](https://resolver.caltech.edu/CaltechAUTHORS:20190328-093304967)
- "Ab initio electron-defect interactions using Wannier functions", npj Comput. Mater. (2020) (Bernardi group; exact author list and volume UNVERIFIED in this session). Wannier interpolation of electron–defect matrix elements for **neutral vacancies in Si and Cu**. This is the closest methodological analogue of the student's Wannier-interpolated M, and it is non-magnetic. — [npj Comput. Mater.](https://www.nature.com/articles/s41524-020-0284-y)
- Yang-hao Chan, Jonah B. Haber, Mit H. Naik, Diana Y. Qiu, Felipe H. da Jornada, "Exciton-defect interaction and optical properties from a first-principles T-matrix approach", arXiv:2505.15523 (2025; journal UNVERIFIED). A disorder-averaged Green's function in the T-matrix approximation for exciton–defect bound states, giving absorption and PL spectra for **MoS₂**. No spin or magnetism appears in the abstract. — [arXiv](https://arxiv.org/abs/2505.15523)
- S. Yuan, H. De Raedt, M. I. Katsnelson, "Modeling electronic structure and transport properties of graphene with resonant scattering centers", Phys. Rev. B 82, 115448 (2010), arXiv:1007.3930 (verified bibliographically). A large-scale TB model of resonant scatterers, vacancies included. From memory it computes the optical conductivity, but the fetched snippet did not confirm this (UNVERIFIED). It is non-magnetic TB, not ab initio. — [arXiv](https://arxiv.org/pdf/1007.3930)
- Kochan et al. 2014 name vacancies as possible magnetic resonant scatterers but compute only hydrogen, with DFT-fitted TB. — [ar5iv](https://ar5iv.labs.arxiv.org/html/1306.0230)
- Soriano et al. 2011 do spin-dependent Kubo transport with DFT-derived mean-field Hubbard parameters, for hydrogen only. — [arXiv](https://arxiv.org/abs/1105.1005)
- Practical enabler (QE docs): in `pp.x`, `plot_num=1` is "total potential V_bare + V_H + V_xc", and for this plot_num **`spin_component` = 0 gives the spin-averaged potential (default), 1 the spin-up potential, 2 the spin-down potential**. Spin-resolved defect potentials Ṽ^σ = V^σ_defective − V_pristine can therefore be extracted directly from nspin=2 runs. — [QE INPUT_PP](https://www.quantum-espresso.org/Doc/INPUT_PP.html)

### Inferences
- **What would be new**, as a combination:
  1. Spin-resolved defect matrix elements M^σ (local plus KB non-local) built from a relaxed, spin-polarized vacancy supercell in plane-wave DFT.
  2. A per-spin T-matrix beyond Born, T^σ = V^σ + V^σ G0 T^σ, in the bulk disorder-averaged setting, giving spin-resolved linewidths Γ^σ(ε) and the exchange splitting of the Vπ resonance.
  3. Spin-dependent resistivity in the two-current picture and spin-resolved optical conductivity σ↑ + σ↓.
  4. Wannier interpolation of M^σ.
  Each piece has a precedent: Kaasbjerg for (2) without spin, Kochan/Soriano for spin with fitted models of hydrogen, Bernardi for (4) without spin or T-matrix. **I found no work that combines them for the vacancy.**
- **A useful direct comparison**: rerun the unpolarized (Kaasbjerg-like) vacancy T-matrix and the ↑/↓ versions on the same supercell. The unpolarized vacancy is a near-E_F resonant scatterer, and exchange splits the resonance by tens of meV, as STM finds (Zhang 2016). This should change the energy dependence of Γ near the Dirac point strongly, while at energies far from the resonance (optical energies of ~1–2 eV) the ↑/↓ difference is likely small. This is a hypothesis to test, not a sourced claim.
- **Caveats the student should state**:
  1. The π moment and its splitting depend on supercell, smearing and k-points (Wang & Pantelides) and on the functional (Valencia & Caldas: GGA < 2 μB, hybrid = 2 μB). They may vanish at low concentration (Palacios & Ynduráin) or with out-of-plane distortion (Padmanabhan & Nanda; Paz et al.).
  2. A collinear treatment is static mean-field. It misses Kondo screening, whose T_K is contested (Chen 2011 vs Jobst & Weber 2012; Jiang 2018), and spin-flip scattering (Kochan 2014).
  3. The spin-1 vs spin-½ question (Casartelli vs Nair) is outside collinear DFT.
  4. The ↑ and ↓ defect potentials both use the same (non-magnetic) pristine reference, so Ṽ^↑ − Ṽ^↓ = V^↑_d − V^↓_d is purely the exchange-correlation magnetization potential of the defective cell.
- Project-specific note (inference): the repo's CLAUDE.md records that Γ-only N×N supercells fold K onto Γ only for N = 3m. For the vacancy moment, whose π part is sensitive to the states at E_F, the honest comparison family is again N = 6, 9, 12. The existing relaxed nspin=1/2 runs are 9×9 (per the project's campaign R1).

### Reference verification ledger
| Reference | Status |
|---|---|
| Yazyev & Helm, PRB 75, 125408 (2007) | Verified |
| Lehtinen et al., PRL 93, 187202 (2004) | Verified (no numeric moment in abstract) |
| Nanda et al., NJP 14, 083004 (2012), arXiv:1105.1129 | Verified |
| Palacios & Ynduráin, PRB 85, 245443 (2012) | Verified |
| Wang & Pantelides, PRB 86, 165438 (2012) | Verified |
| Valencia & Caldas, PRB 96, 125431 (2017) | Verified |
| Padmanabhan & Nanda, PRB 93, 165403 (2016) | Verified |
| Paz, Scopel, Freitas, SSC 2013, DOI 10.1016/j.ssc.2013.05.004 | Partially verified (volume/pages UNVERIFIED) |
| Casartelli et al., PRB 88, 195424 (2013) | Verified |
| Cazalilla et al., arXiv:1207.3135 | Verified as arXiv; journal UNVERIFIED |
| Mitchell & Fritz, PRB 88, 075104 (2013) | Verified |
| May et al., PRB 97, 155419 (2018) | Verified |
| Nair et al., Nat. Phys. 8, 199 (2012) | Verified |
| Nair et al., Nat. Commun. 4, 2010 (2013) | Verified |
| Zhang et al., PRL 117, 166801 (2016) | Verified |
| Ugeda et al., PRL 104, 096804 (2010) | Verified |
| Chen et al., Nat. Phys. 7, 535 (2011) | Verified |
| Jobst & Weber, Nat. Phys. 8, 352 (2012) | Verified |
| Jobst, Kisslinger, Weber, PRB 88, 155412 (2013) | Verified |
| Jiang et al., Nat. Commun. 9, 2349 (2018) | Verified |
| Chen, Wu, Yu, Liao, Nanoscale 6, 8814 (2014) | Verified |
| Kochan, Gmitra, Fabian, PRL 112, 116602 (2014) | Verified |
| Kochan, Irmer, Gmitra, Fabian, PRL 115, 196601 (2015) | Verified |
| Soriano et al., PRL 107, 016602 (2011) | Verified |
| Soriano et al., 2D Mater. 2, 022002 (2015) | Verified |
| Thomsen et al., PRB 92, 195408 (2015) | Verified |
| Fedorov et al., PRL 110, 156602 (2013) | Verified (method UNVERIFIED) |
| Wehling et al., PRL 105, 056802 (2010) | Verified |
| Topsakal et al., PRB 78, 235435 (2008) | Verified |
| Kaasbjerg, PRB 101, 045433 (2020) | Verified (spin-unpolarized, via ar5iv) |
| Lu, Zhou, Bernardi, PRMaterials 3, 033804 (2019) | Verified |
| npj Comput. Mater. 2020 Wannier electron-defect (Bernardi group) | Title/venue seen; authors/volume UNVERIFIED |
| Chan et al., arXiv:2505.15523 (2025) | Verified as arXiv; journal UNVERIFIED |
| Yuan, De Raedt, Katsnelson, PRB 82, 115448 (2010) | Verified; optical-conductivity content UNVERIFIED |
| Zhang, Jia, Chen, He, Wang, review DOI 10.1088/1674-1056/ac70c4 | Verified DOI; volume UNVERIFIED |
| Sci. Rep. 2017 shear-distortion vacancy paper | Venue verified; authors UNVERIFIED |
| Goswami & Gahlot, arXiv:1205.1230 | Verified as arXiv; journal UNVERIFIED |
| arXiv:1703.06504 (TranSIESTA GNR vacancy spin current) | Snippet only; details UNVERIFIED |

### Gaps
- It is not proven that nobody has done this. Searches (standard and extended, 2020–2026 keywords) found no spin-resolved ab initio T-matrix, beyond-Born spin rate, spin-dependent resistivity, or spin-resolved optical conductivity for the graphene monovacancy. Full citing-list sweeps of Kaasbjerg 2020, Yazyev & Helm 2007 and Kochan 2014 were not possible, so a check in Google Scholar or Web of Science "cited by" lists before claiming priority is recommended.
- Places that might hide a precedent:
  - KKR Green-function impurity studies (Mertig/Zeller-type codes naturally give spin-polarized single-impurity T-matrices).
  - DFT + DMFT or NRG studies of vacancy Kondo resistivity.
  - Theses.
  - Non-English journals.
- No ab initio optical-conductivity study of graphene with magnetic vacancies, spin-resolved or otherwise, was found. The nearest items are the non-magnetic TB work of Yuan et al. 2010 (content not confirmed) and the MoS₂ exciton–defect T-matrix of Chan et al. 2025.
