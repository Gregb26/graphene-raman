# Model-level theory of σ(ω) in graphene with vacancies / resonant scatterers

Scope: tight-binding (TB) and continuum (Dirac) theory of the IR-to-visible optical conductivity of graphene with
vacancies or resonant scatterers (H adatoms, organic groups), the approximations used, the predicted features, vertex
corrections and T-matrix validity, and whether anyone has fed a DFT-derived defect T-matrix into a Kubo σ(ω).

Verification method (2026-09-29): every bibliographic entry below was checked against the Crossref API (authors,
title, journal, volume, article number, year, DOI). Physics content was read from the arXiv preprint full text
(PDF text extraction) unless marked "title only" (content not checked in this session). arXiv versions can differ
slightly from the published versions.

## Q1. Which papers computed σ(ω) of graphene with vacancies / resonant scatterers?

### Takeaway
σ(ω) of graphene with vacancies or resonant scatterers has been studied since 2006. Two families: (i) analytic or
semi-analytic Dirac/TB models with a self-consistent, k-independent self-energy (full self-consistent Born/CPA)
inserted in a Kubo bubble: Peres–Guinea–Castro Neto 2006, Peres–Stauber–Castro Neto 2008, Stauber–Peres–Castro Neto
2008. (ii) Numerically exact real-space Kubo on π-band TB lattices of 10^7–10^8 sites (time evolution or Chebyshev/KPM):
Yuan–De Raedt–Katsnelson 2010, Yuan–Roldán–De Raedt–Katsnelson 2011, Cysne et al. 2016. Most of the other papers the
user named (Ferreira 2011, Ostrovsky 2006/2010, Wehling 2010, Skrypnyk–Loktev 2010, Pellegrino 2009, Kaasbjerg 2020)
treat only DOS or dc transport, not σ(ω).

### Cited Findings

**Papers that compute σ(ω) with vacancies / resonant scatterers**

- **Peres, Guinea, Castro Neto, "Electronic properties of disordered two-dimensional carbon", PRB 73, 125411 (2006)**, doi:10.1103/PhysRevB.73.125411, arXiv:cond-mat/0512091. Honeycomb lattice with vacancies. The self-energy is computed in the "full Born approximation" (FBA, single-vacancy T-matrix, non-self-consistent) and the "full self-consistent Born approximation" (FSBA: bare propagator in the T-matrix replaced by the full one, Σ(ω) = −n_i Ḡ⁰_AA(ω − Σ(ω))). The FSBA self-energy is momentum independent. σ(ω,T) comes from a Kubo bubble with the current expanded around K, built only from Im G_AA · Im G_AA (the paper states that products G_AB G_BA do not contribute). No vertex correction is discussed. Results: a universal T=0 dc value σ₀ = (2/π)(e²/h)·(1 − [ImΣ(0)]²/(D² + [ImΣ(0)]²)) ≈ (2/π)e²/h, independent of dilution when ImΣ(0) ≪ D (spin/valley and h vs ħ conventions are ambiguous in the preprint: the Fig. 8 caption says units e²/(πħ), the axis says πσh/e²). At low T, σ(ω) develops a maximum at a frequency that depends on n_i and "almost shows scaling" in ω/√n_i. At higher T it becomes Drude-like. Parameters: D = 1–2 eV cutoff, n_i = 10⁻⁴–10⁻¹, T = 10–300 K. The paper also gives the magneto-optical σ_xx(ω) in the FSBA. — [arXiv:cond-mat/0512091](https://arxiv.org/abs/cond-mat/0512091); [DOI](https://doi.org/10.1103/PhysRevB.73.125411)
- **Peres, Stauber, Castro Neto, "The infrared conductivity of graphene on top of silicon oxide", EPL 84, 38002 (2008)** (arXiv title "The infrared conductivity of graphene"), doi:10.1209/0295-5075/84/38002, arXiv:0803.2816. IR σ(ω) at finite μ and T with charged impurities (RPA screening combined with CPA), unitary scatterers ("edge defects, cracks, vacancies") and phonons (2nd-order perturbation theory). Essentially one free parameter, the charged-impurity density n_i^C, with 1 ≫ n_i^C ≫ n_i → 0. In their model the unitary scatterers only produce a finite DOS at the Dirac point. They find "an anomalous enhancement of the conductivity in a frequency region that is blocked by Pauli exclusion and an impurity broadening of the conductivity threshold", plus phonon Stokes/anti-Stokes lines above σ₀ = (π/2)e²/h. — [arXiv:0803.2816](https://arxiv.org/abs/0803.2816); [DOI](https://doi.org/10.1209/0295-5075/84/38002)
- **Stauber, Peres, Castro Neto, "Conductivity of suspended and non-suspended graphene at finite gate voltage", PRB 78, 085418 (2008)**, arXiv:0809.2578. dc and optical σ with midgap states (unitary scatterers) in the CPA, charged impurities (2nd-order self-energy) and phonons. Parameters: t = 3 eV, cutoff D = 7 eV. The CPA DOS is benchmarked against a numerically exact method for n_i = 0.005 and 0.01. Their fit to experiment uses T = 45 K, n_i(unitary) = 4.0×10⁻⁵, n_i^C = 1.3×10⁻⁴. With midgap states ten times fewer than Coulomb scatterers, the sub-2μ conductivity is "mainly controlled by phonons and charged impurities". They note that experiments show a plateau σ ≈ σ₀/3 for k_BT ≪ ħω ≪ 2μ and an unusually broad 2μ edge. On dc: σ(μ=0) → 4e²/(πh) when the broadening limit Γ→0 is taken last, and the same value is "also obtained without the limit Γ→0 within the self-consistent CPA"; π e²/(2h) results in the other order of limits. — [arXiv:0809.2578](https://arxiv.org/abs/0809.2578); [DOI](https://doi.org/10.1103/PhysRevB.78.085418)
- **Yuan, De Raedt, Katsnelson, "Modeling electronic structure and transport properties of graphene with resonant scattering centers", PRB 82, 115448 (2010)**, arXiv:1007.3930. Noninteracting π-band TB with vacancies (site removal) and "hydrogen" resonant impurities (Anderson impurity with ε_d = −t/16, V = 2t). Method: time-dependent Schrödinger equation with Chebyshev propagation and random superposition states, on samples of 4096² to 8192² sites. Computes DOS, quasieigenstates, and ac and dc conductivity (Kubo, regular part only; the Drude weight is omitted). Current operator J = −(ie/ħ) Σ t_ij (r_j − r_i) c†_i c_j, i.e. diagonal position operator, no intra-atomic dipole. — [arXiv:1007.3930](https://arxiv.org/abs/1007.3930); [DOI](https://doi.org/10.1103/PhysRevB.82.115448)
- **Yuan, Roldán, De Raedt, Katsnelson, "Optical conductivity of disordered graphene beyond the Dirac cone approximation", PRB 84, 195418 (2011)**, arXiv:1109.3485. This is the key model paper. Full π-band TB with Kubo σ(ω) by the same time-evolution method: 8192×8192 sites for σ, 4096×4096 for the DOS, T = 300 K, periodic boundary conditions, regular part only. Disorder types: random or Gaussian on-site potentials, random or Gaussian hoppings, random vacancies, random H adatoms, and vacancy/H clusters. Features are listed in Q2. The H parameters V ≈ 2t, ε_d ≈ −t/16 are "obtained from ab initio DFT calculations" (i.e. Wehling et al. 2010). — [arXiv:1109.3485](https://arxiv.org/abs/1109.3485); [DOI](https://doi.org/10.1103/PhysRevB.84.195418)
- **Cysne, Rappoport, Ferreira, Lopes, Peres, "Numerical calculation of the Casimir–Polder interaction between a graphene sheet with vacancies and an atom", PRB 94, 235405 (2016)**, arXiv:1608.04368. Complex σ(ω) (real and imaginary parts) of NN-TB graphene (t = 2.7 eV) with compensated random vacancies (both sublattices), n_i = 0.4%, on a 3200×3200 lattice. Method: exact Chebyshev expansion of the Kubo–Bastin-type formula at T = 0 (8000² Chebyshev moments, 5000 random vectors, broadening η ≈ 8 meV, a single disorder realization). Results in Q2. — [arXiv:1608.04368](https://arxiv.org/abs/1608.04368); [DOI](https://doi.org/10.1103/PhysRevB.94.235405)
- **Yuan, Wehling, Lichtenstein, Katsnelson, "Enhanced screening in chemically functionalized graphene", PRL 109, 156601 (2012)**, arXiv:1205.2782. Kubo-formula dielectric and polarization function (not σ(ω) per se) with resonant H adatoms. Midgap states increase ε, but the static ε does not diverge metallically as q→0 ("bad metal"). A new length l_c beyond which screening is suppressed is identified with the Anderson localization length. — [arXiv:1205.2782](https://arxiv.org/abs/1205.2782); [DOI](https://doi.org/10.1103/PhysRevLett.109.156601)

**Related papers checked: dc/DOS only, or only a σ(ω) formula**

- **Ferreira, Viana-Gomes, Nilsson, Mucciolo, Peres, Castro Neto, PRB 83, 165402 (2011)**, arXiv:1010.4026. **dc only.** Shows that partial waves, Lippmann–Schwinger and T-matrix give equivalent results for resonant scatterers. KPM Kubo near neutrality gives a minimum-conductivity plateau of "about e²/h (per layer)". — [arXiv:1010.4026](https://arxiv.org/abs/1010.4026); [DOI](https://doi.org/10.1103/PhysRevB.83.165402)
- **Pereira, Lopes dos Santos, Castro Neto, "Modeling disorder in graphene", PRB 77, 115109 (2008)**, arXiv:0712.0806. DOS only. Random, compensated and uncompensated vacancies, local and substitutional impurities give "localized zero modes, strong resonances, gap and pseudogap behavior, and non-dispersive midgap zero modes". — [arXiv:0712.0806](https://arxiv.org/abs/0712.0806); [DOI](https://doi.org/10.1103/PhysRevB.77.115109)
- **Pereira, Guinea, Lopes dos Santos, Peres, Castro Neto, "Disorder induced localized states in graphene", PRL 96, 036801 (2006)**. Title only. — [DOI](https://doi.org/10.1103/PhysRevLett.96.036801)
- **Ostrovsky, Gornyi, Mirlin, "Electron transport in disordered graphene", PRB 74, 235443 (2006)**, arXiv:cond-mat/0609617. dc only, self-consistent T-matrix approximation (SCTMA) in the Dirac model. Strong scatterers give conductivity linear in carrier density; weak ones give a logarithmic dependence. At half filling σ ~ e²/h if the disorder preserves a chiral symmetry. The vertex-correction results are in Q3. — [arXiv:cond-mat/0609617](https://arxiv.org/abs/cond-mat/0609617); [DOI](https://doi.org/10.1103/PhysRevB.74.235443)
- **Wehling, Yuan, Lichtenstein, Geim, Katsnelson, "Resonant scattering by realistic impurities in graphene", PRL 105, 056802 (2010)**, arXiv:1003.0609. dc only (Boltzmann plus numerically exact Kubo). DFT-derived impurity parameters, see Q4. — [arXiv:1003.0609](https://arxiv.org/abs/1003.0609); [DOI](https://doi.org/10.1103/PhysRevLett.105.056802)
- **Skrypnyk & Loktev, "Electrical conductivity in graphene with point defects", PRB 82, 085436 (2010)**, arXiv:1004.4606. dc at T = 0. Binary-alloy (Lifshitz) model. The "modified propagator" self-energy is "practically indistinguishable" from the CPA at c ≪ 1 and is inserted into a bubble Kubo formula σ = (4e²/πh){1 + […] arctan[…]}. The minimal conductivity "does not have a universal character and corresponds to the impurity resonance energy rather than to the Dirac point", with a gate-voltage asymmetry. They state that the Kubo formula and the propagator method are not reliable in Ioffe–Regel localized ranges. — [arXiv:1004.4606](https://arxiv.org/abs/1004.4606); [DOI](https://doi.org/10.1103/PhysRevB.82.085436)
- **Skrypnyk & Loktev, "Impurity effects in a two-dimensional system with the Dirac spectrum", PRB 73, 241402 (2006)**; **"Local spectrum rearrangement in impure graphene", PRB 75, 245401 (2007)**; **Pershoguba, Skrypnyk, Loktev, "Numerical evidence of spectrum rearrangement in impure graphene", PRB 80, 214201 (2009)**. Title only. — [DOI 73.241402](https://doi.org/10.1103/PhysRevB.73.241402); [DOI 75.245401](https://doi.org/10.1103/PhysRevB.75.245401); [DOI 80.214201](https://doi.org/10.1103/PhysRevB.80.214201)
- **Pogorelov, Loktev, Kochan, "Impurity resonance effects in graphene vs impurity location, concentration and sublattice occupation"**, arXiv:2005.12643 (2020); journal version UNVERIFIED. Lifshitz and Anderson models, top/bridge/hollow adatoms: resonances, quasi-gaps, impurity sub-bands, Mott-like transitions (Ioffe–Regel–Mott criterion). Optics is only mentioned as a possible implication. — [arXiv:2005.12643](https://arxiv.org/abs/2005.12643)
- **Pellegrino, Angilella, Pucci, PRB 80, 094203 (2009)**, arXiv:0909.1903. Site, bond or hollow impurity, T-matrix/full Born approximation (dilute, n_imp² ≪ n_imp). Writes the Kubo σ(μ,T;ω) but reports only the lifetime and the **static** conductivity (sublinear in carrier density). — [arXiv:0909.1903](https://arxiv.org/abs/0909.1903); [DOI](https://doi.org/10.1103/PhysRevB.80.094203)
- **Stauber, Peres, Guinea, PRB 76, 205423 (2007)**, "Electronic transport in graphene: A semiclassical approach including midgap states". Title only. — [DOI](https://doi.org/10.1103/PhysRevB.76.205423)
- **Peres, "Colloquium: The transport properties of graphene: An introduction", RMP 82, 2673–2700 (2010)**. Title only; whether it treats σ(ω) with vacancies was not checked. — [DOI](https://doi.org/10.1103/RevModPhys.82.2673)
- **Stauber, Peres, Geim, PRB 78, 085432 (2008)**, arXiv:0803.1802. **Clean** graphene beyond the Dirac cone (NN and NNN hopping), visible range. Corrections to σ₀ are "surprisingly small (a few percent)". No disorder, so it is a clean-limit reference only. — [arXiv:0803.1802](https://arxiv.org/abs/0803.1802); [DOI](https://doi.org/10.1103/PhysRevB.78.085432)
- **Gusynin & Sharapov, PRB 73, 245411 (2006)**, arXiv:cond-mat/0512157. Analytic ac and Hall conductivities of Dirac quasiparticles in a B field. Not specific to resonant scatterers. — [arXiv](https://arxiv.org/abs/cond-mat/0512157); [DOI](https://doi.org/10.1103/PhysRevB.73.245411)
- **Gusynin, Sharapov, Carbotte, "Sum rules for the optical and Hall conductivity in graphene", PRB 75, 165407 (2007)**. Title only. — [DOI](https://doi.org/10.1103/PhysRevB.75.165407)

**Other explanations of the sub-2μ IR absorption (not resonant scatterers)**

- Grushin, Valenzuela, Vozmediano, "Effect of Coulomb interactions on the optical properties of doped graphene", PRB 80, 155417 (2009). Title only. Yuan 2011 summarizes the candidate mechanisms as "disorder, electron-electron interaction or excitonic effects". — [DOI](https://doi.org/10.1103/PhysRevB.80.155417); [Yuan 2011](https://arxiv.org/abs/1109.3485)
- Kechedzhi & Das Sarma, "Plasmon anomaly in the dynamical optical conductivity of graphene", PRB 88, 085403 (2013). Charged impurities with plasmon-dressed screening give a frequency-dependent scattering rate and "a broad peak in the frequency-dependent graphene optical conductivity". — [arXiv:1305.4940](https://arxiv.org/abs/1305.4940); [DOI](https://doi.org/10.1103/PhysRevB.88.085403)
- Principi, Vignale, Carrega, Polini, "Impact of disorder on Dirac plasmon losses", PRB 88, 121405 (2013). Plasmon lifetime ≠ τ_tr, charged impurities. — [arXiv:1307.7371](https://arxiv.org/abs/1307.7371); [DOI](https://doi.org/10.1103/PhysRevB.88.121405)
- Experimental benchmark: Li, Henriksen, Jiang, Hao, Martin, Kim, Stormer, Basov, "Dirac charge dynamics in graphene by infrared spectroscopy", Nat. Phys. 4, 532–535 (2008). — [DOI](https://doi.org/10.1038/nphys989)

**Minimal / dc conductivity with resonant scatterers (for the σ(ω→0) end)**

- Titov, Ostrovsky, Gornyi, Schuessler, Mirlin, PRL 104, 076802 (2010): in undoped ballistic samples, each resonant impurity adds 16e²/(π²h) to the conductance. — [arXiv:0908.3793](https://arxiv.org/abs/0908.3793); [DOI](https://doi.org/10.1103/PhysRevLett.104.076802)
- Ostrovsky, Titov, Bera, Gornyi, Mirlin, PRL 105, 266803 (2010): for smooth resonant impurities σ grows logarithmically with concentration (class DIII). For vacancies or strong on-site impurities it saturates at a constant that depends on the sublattice distribution (class BDI). — [arXiv:1006.3299](https://arxiv.org/abs/1006.3299); [DOI](https://doi.org/10.1103/PhysRevLett.105.266803)
- Ferreira & Mucciolo, PRL 115, 106601 (2015): Kubo dc of vacancy-induced zero modes on lattices of more than 10⁹ sites equals 4e²/(πh) within 1%, over a wide window of broadenings and concentrations. — [arXiv:1507.00488](https://arxiv.org/abs/1507.00488); [DOI](https://doi.org/10.1103/PhysRevLett.115.106601)
- Cresti, Ortmann, Louvet, Van Tuan, Roche, PRL 110, 196601 (2013): the Dirac-point conductivity is "supermetallic" (robust) or suppressed by a gap or algebraic localization, depending on the sublattice balance of the vacancies and the geometry. Weak and Anderson localization at higher energies. — [arXiv:1304.8061](https://arxiv.org/abs/1304.8061); [DOI](https://doi.org/10.1103/PhysRevLett.110.196601)
- Yuan 2010: a plateau of order 4e²/(πh) near neutrality. At a few % impurities the impurity band conducts and can give a maximum of σ in the midgap region. Away from neutrality, a Boltzmann resonant form σ ∝ (e²/h)(n_e/n_i) ln²|E_F/D| is consistent with the numerics. — [arXiv:1007.3930](https://arxiv.org/abs/1007.3930)

### Inferences
- The main model benchmarks for σ(ω) with vacancies or H adatoms are Yuan et al. 2011 (full spectrum, doped and undoped, exact real-space Kubo) and Cysne et al. 2016 (complex σ, dilute vacancies, KPM). The CPA/FSBA papers (2006–2008) give the analytic Dirac-cone picture with a bubble and a k-independent self-energy, which is the closest formal analogue of the user's "bubble + Σ = cT̄" scheme.
- Real-space exact Kubo methods (time evolution, KPM) contain all disorder diagrams (vertex corrections, crossed diagrams, localization) up to finite size and broadening. Comparing to them therefore tests the no-vertex, non-self-consistent approximation directly.

### Gaps
- Content of Skrypnyk–Loktev 2006/2007, Pershoguba 2009, Pereira 2006, Stauber–Peres–Guinea 2007 and the Peres RMP 2010 was not checked (bibliographic data only). The RMP may contain a σ(ω) discussion with midgap states.
- I did not find a Low Temp. Phys. or other paper by the Loktev/Pogorelov group computing σ(ω) explicitly. The journal version of arXiv:2005.12643 is UNVERIFIED.
- A systematic citation-forward search of Yuan 2011 and Peres 2006 (e.g. post-2016 KPM/KITE optical studies with vacancies) was not completed within the budget. KITE (João et al., R. Soc. Open Sci. 7, 191809 (2020), [DOI](https://doi.org/10.1098/rsos.191809)) is a verified tool paper for such response functions, but its demo cases were not checked.

## Q2. What features are predicted?

### Takeaway
Resonant scatterers (vacancies, H/CH₃ adatoms) are predicted to produce:
- a finite IR background inside the Pauli-blocked window 0 < ħω < 2μ (up to ≈0.4σ₀ at 0.5% impurities, T = 300 K);
- an extra peak at ħω ≈ μ (impurity band → Fermi level);
- a new peak at ħω ≈ t ≈ 2.7 eV (midgap states ↔ van Hove singularities), growing with concentration;
- smearing and damping of the ħω ≈ 2t van Hove peak (fully smeared by ~1% H);
- broadening of the 2μ edge, and wash-out of the Drude peak below a critical μ_c (≈0.15 eV at 0.4% vacancies).

At μ = 0, the FSBA predicts a low-frequency maximum scaling as ω/√n_i. Non-resonant (Anderson-type or hopping) disorder barely changes the IR part.

### Cited Findings
- **Sub-2μ background (resonant scatterers).** Yuan 2011: for μ = 0.1t and 0.2t, all disorder types give a peak near ω = 0 (intraband broadening). At slightly higher ω, σ "drops to almost zero for the case of non-resonant scatterers", while vacancies and H keep "a non-zero background". This is attributed to "transitions between the newly formed impurity band and the conduction band". — [Yuan et al. 2011](https://arxiv.org/abs/1109.3485)
- **Concentration and μ scaling of the background.** Yuan 2011: σ in 0 < ω < 2μ grows with n_x(i) "from σ(ω)=0 for a clean sample to σ(ω) ≈ 0.4σ₀ for the larger concentration of impurities considered (0.5%)". σ(ω < 2μ) increases as doping decreases. About 0.25% resonant impurities reproduces the background of Li et al. (graphene on SiO₂), and about 0.1% reproduces Chen et al. (ion-gel gating). Correlated Gaussian potential clusters (puddles) give similar IR backgrounds, so the IR background alone does not identify the defect type. — [Yuan et al. 2011](https://arxiv.org/abs/1109.3485)
- **Peak at ħω ≈ μ.** Yuan 2011: "the peak observed in σ(ω) for the case of resonant impurities at the energy ω ≈ μ is associated to transitions between the above discussed impurity band and states at the Fermi level." — [Yuan et al. 2011](https://arxiv.org/abs/1109.3485)
- **New peak at ħω ≈ t (undoped).** Yuan 2011: vacancies or H at n_x = 1–10% produce "a new peak at an energy ω ≈ t … associated to optical transitions between the newly formed midgap states (E≈0) and the states of the Van Hove singularities (E≈t)". Unlike the 2t peak, "the height of this ω≈t peak grows with the strength of disorder". — [Yuan et al. 2011](https://arxiv.org/abs/1109.3485)
- **Van Hove peak at ħω ≈ 2t.** Smeared by every disorder type, proportionally to disorder strength (Yuan 2011). With H adatoms, "the Van Hove singularity of clean graphene is smeared out completely for concentrations as small as 1%" (Yuan 2010). The graphane limit (100% H) opens a 2t gap with σ = 0 for |ω| < 2t. — [Yuan 2011](https://arxiv.org/abs/1109.3485); [Yuan 2010](https://arxiv.org/abs/1007.3930)
- **Undoped IR modulation.** Yuan 2011: at μ = 0 the midgap DOS peak modulates the IR σ, with "lower peaks … σ ≈ 0.9σ₀" attributed to excitations involving states around the zero modes. Random on-site or hopping disorder without midgap states leaves the low-energy σ ≈ σ₀ ("almost no effect … unless the disorder is too large"). — [Yuan et al. 2011](https://arxiv.org/abs/1109.3485)
- **Vacancy vs H equivalence.** H adatoms (V ≈ 2t, ε_d ≈ −t/16) give DOS and σ(ω) similar to vacancies. The central DOS peak is shifted by ε_d < 0. Clusters of vacancies or H give fewer midgap states than isolated defects at equal total concentration. — [Yuan 2011](https://arxiv.org/abs/1109.3485); [Yuan 2010](https://arxiv.org/abs/1007.3930)
- **Drude wash-out and sign of Im σ.** Cysne 2016 (n_i = 0.4% compensated vacancies, T = 0, η ≈ 8 meV):
  - At μ = 0.5 eV there is a well-defined step at 2μ, and the Drude part is fit with a single parameter ħ/τ ≈ 0.07 eV.
  - As μ decreases the Fermi step blurs; at μ = 0.2 eV "there is no trace of the Fermi step".
  - Below a critical μ_c ≈ 0.15 eV the Drude peak is "completely washed out", Re σ at low ω is strongly suppressed, and Im σ < 0 over the whole frequency range (no TM plasmons).
  — [Cysne et al. 2016](https://arxiv.org/abs/1608.04368)
- **μ = 0 low-frequency behaviour (FSBA, vacancies).** Peres 2006: at low T, σ(ω) has a maximum at a frequency set by n_i, with approximate scaling in ω/√n_i. At higher T it becomes Drude-like. The dc limit is universal, ≈ (2/π)e²/h in the paper's convention. — [Peres et al. 2006](https://arxiv.org/abs/cond-mat/0512091)
- **Edge smearing.** CPA+RPA gives "impurity broadening of the conductivity threshold" at 2μ and an "anomalous enhancement" below it (Peres EPL 2008). Experiments show a broad 2μ edge and a plateau of ≈σ₀/3 below it (quoted in Stauber 2008). — [Peres EPL 2008](https://arxiv.org/abs/0803.2816); [Stauber 2008](https://arxiv.org/abs/0809.2578)
- **Minimal dc values.** 4e²/(πh):
  - numerics for vacancies within 1% (Ferreira–Mucciolo 2015);
  - plateau (Yuan 2010);
  - CPA (Stauber 2008);
  - about e²/h per layer (Ferreira 2011 KPM);
  - non-universal and pinned at the impurity-resonance energy when the resonance is shifted from E_D (Skrypnyk–Loktev 2010).
  — [Ferreira–Mucciolo](https://arxiv.org/abs/1507.00488); [Yuan 2010](https://arxiv.org/abs/1007.3930); [Stauber 2008](https://arxiv.org/abs/0809.2578); [Ferreira 2011](https://arxiv.org/abs/1010.4026); [Skrypnyk–Loktev 2010](https://arxiv.org/abs/1004.4606)
- **Position of the resonance.** DFT for CH₃, C₂H₅, CH₂OH, H, OH puts the midgap states within ±0.03 eV of neutrality (Wehling 2010). Resonant scattering is generic for short-range impurities whose potential is of order the bandwidth or larger (Basko 2008). — [Wehling 2010](https://arxiv.org/abs/1003.0609); [Basko 2008](https://arxiv.org/abs/0806.2785)

### Inferences
- A calculation that is correct at the model level should reproduce, at n ≈ 0.1–1% and μ ≈ 0.1–0.5 eV:
  - a sub-2μ background of order 0.1–0.4σ₀;
  - a feature at ħω ≈ μ;
  - a Drude width of order ħ/τ ~ 0.1 eV at 0.4%;
  - (undoped, n ≳ 1%) a peak near ħω ≈ t and a damped 2t peak.

  These are natural sanity checks for the user's ab initio σ(ω).
- The ħω ≈ μ peak and the IR background come from transitions into or out of the quasi-localized impurity band near E_D. They are therefore controlled by exactly the energy region where a non-self-consistent T-matrix self-energy is least reliable (see Q3).
- Yuan's numbers are at T = 300 K with random vacancies on both sublattices. Cysne's are at T = 0 with compensated vacancies. The user should match T, sublattice distribution and broadening before comparing magnitudes.

### Gaps
- Yuan 2011 gives no closed-form scaling of the ħω ≈ μ peak height or of the background with n. Only the discrete concentrations 0.05–0.5% (doped) and 1–10% (undoped) were read.
- No source found giving an analytic (T-matrix) prediction of the ħω ≈ μ peak; it appears only in the exact numerics.

## Q3. Vertex corrections; validity of non-self-consistent T-matrix vs CPA/SCTMA; pitfalls

### Takeaway
For atomically sharp (intervalley) short-range scatterers in graphene, the ladder vertex correction to the current is absent within SCTMA (Ostrovsky et al. 2006). Long-range (intravalley) disorder carries a factor-2 vertex correction from chirality. All model σ(ω) papers with a self-energy (Peres 2006/2008, Stauber 2008, Skrypnyk–Loktev 2010) use a bare bubble. I found no paper that quantifies the vertex correction to σ(ω) for resonant scatterers; exact real-space Kubo includes it implicitly. The non-self-consistent T-matrix (Σ = c⟨T⟩) is exact to first order in c (Kaasbjerg 2020) but breaks down near the Dirac point / resonance. There the vacancy T-matrix gives Im Σ ∝ c/(|E| ln²(D/|E|)); self-consistency (FSBA/SCTMA/CPA) regularizes this below a scale Γ_c ∝ √c (up to logs), and at still lower energies even SCTMA fails (zero modes, Gade-type singularities).

### Cited Findings
- **Vertex correction, short- vs long-range (Dirac/SCTMA).** Ostrovsky 2006: "The sum of the ladder diagrams … give the correction to the current vertex … In the limit of the short-range potential disorder, we have V = 1. In the opposite long-range case, the summation of ladder diagrams yields V = 1/(1 − n_imp U²|1−Ug|⁻² Π_RA)". For unitary scatterers: σ_SCUA(ε) = 4e²ε²/(π²ηΔ²) log²(Δ/|ε|) for long-range disorder; "if the disorder is short-range, the vertex correction is absent and the resulting conductivity is twice smaller." — [Ostrovsky et al. 2006](https://arxiv.org/abs/cond-mat/0609617)
- **Bubble-only treatments in σ(ω) papers.**
  - Peres 2006 uses a bubble built only from G_AA products with FSBA Σ, with no vertex discussion.
  - Skrypnyk–Loktev 2010 insert a CPA-like Σ into a closed-form bubble.
  - Stauber 2008 and Peres EPL 2008 use Dyson-equation self-energies (CPA for midgap states) inside the Kubo formula.
  — [Peres 2006](https://arxiv.org/abs/cond-mat/0512091); [Skrypnyk–Loktev 2010](https://arxiv.org/abs/1004.4606); [Stauber 2008](https://arxiv.org/abs/0809.2578)
- **Vertex corrections with chirality (electron–phonon, dc).** Cappelluti & Benfatto, "Vertex renormalization in dc conductivity of doped chiral graphene", PRB 79, 035419 (2009). Explicit Kubo vertex corrections keeping the chiral matrix structure: "at least in the regime of large chemical potential the Boltzmann picture is justified", and this is robust to small sublattice inequivalence. This is for phonons, not defects. — [arXiv:0809.4215](https://arxiv.org/abs/0809.4215); [DOI](https://doi.org/10.1103/PhysRevB.79.035419)
- **Validity of the non-self-consistent T-matrix.** Kaasbjerg 2020: the T-matrix ("full Born") self-energy "takes into account multiple scattering off defects to all orders in the defect potential, and is therefore exact to lowest order in the disorder concentration", and "essentially exact for dilute disorder (c_dis ≪ 1)". Born approximation breaks down for graphene and TMD defects: it overestimates TMD scattering rates "by up to several orders of magnitude". — [arXiv:1911.00530](https://arxiv.org/abs/1911.00530); [DOI](https://doi.org/10.1103/PhysRevB.101.045433)
- **Singularity at E_D (single vacancy).** Peres 2006: the vacancy correction to the DOS is ≈ −(2/N)/(|ω| log²(D/|ω|)) as ω→0, "singular in the low frequency regime". The FBA (non-self-consistent) DOS peak appears at an energy ~n_i D/4. The FSBA includes some multiple-site scattering and yields Σ(ω) = −n_i Ḡ⁰_AA(ω−Σ). — [Peres et al. 2006](https://arxiv.org/abs/cond-mat/0512091)
- **Crossover scale for self-consistency (unitary limit).** Ostrovsky 2006: for ε ≫ Γ_η, Σ(ε) ≈ (ηΔ²/2ε)[1/log(Δ/|ε|) ∓ iπ sgnε/(2 log²(Δ/|ε|))] (this is the non-self-consistent regime). For ε ≪ Γ_η, Σ → ∓iΓ_η and the DOS saturates to a constant. Γ_η depends on η as a power law (∝ √η up to logarithms), unlike the exponentially small Born scale. Beyond first order in η, crossed diagrams matter. SCTMA is "not quantitatively justified in the Born regime" but qualitatively correct far from the degeneracy point. — [Ostrovsky et al. 2006](https://arxiv.org/abs/cond-mat/0609617)
- **SCTMA vs exact numerics at very low energy.** Häfner et al. 2014, compensated vacancies (class BDI): numerics are compatible with a Gade-type singularity ρ(E) ~ |E|⁻¹ exp(−|log E|^{−1/x}) in a pre-asymptotic regime. At even lower energies ρ ~ E⁻¹|log E|^{−𝔵} with 1 ≤ 𝔵 < 2 appears, beyond SCTMA. — [arXiv:1404.6138](https://arxiv.org/abs/1404.6138); [DOI](https://doi.org/10.1103/PhysRevLett.113.186802)
- **CPA benchmark.** Stauber 2008 compares the CPA DOS for midgap states with a numerically exact method at n_i = 0.005 and 0.01; Im Σ increases near E = 0. Skrypnyk–Loktev: at c ≪ 1 the modified-propagator method ≈ CPA, and the Kubo bubble fails in Ioffe–Regel localized windows (around renormalized phase ≈ π/6, 5π/6). — [Stauber 2008](https://arxiv.org/abs/0809.2578); [Skrypnyk–Loktev 2010](https://arxiv.org/abs/1004.4606)
- **Localization and sublattice pitfalls.** At the Dirac point the conductivity depends on the chiral symmetry class and sublattice imbalance (BDI saturation vs DIII log growth; supermetallic vs insulating). This physics is not in any single-impurity T-matrix. — [Ostrovsky 2010](https://arxiv.org/abs/1006.3299); [Cresti 2013](https://arxiv.org/abs/1304.8061); [Ferreira–Mucciolo 2015](https://arxiv.org/abs/1507.00488)

### Inferences
These are my own reasoning, not taken from a source.
- **Symmetry argument for the ladder vertex.** For a scatterer whose T-matrix acts on a single site (s-wave, on-site), the ladder vertex correction at any (ω, z₁, z₂) involves Λ_α = Σ_k ⟨site|G(k,z₁) v_α(k) G(k,z₂)|site⟩. This is a single number that must transform as a vector under the C₃ site symmetry, so it vanishes. The bubble with a local, k-independent Σ is then conserving at the ladder level (the same mechanism as in single-site CPA). This is consistent with Ostrovsky's V = 1 for short-range disorder.
- **The user's extended T-matrix.** A 29-cell cluster, multi-orbital T-matrix from DFT includes p-like (anisotropic) channels, so the ladder vertex need not vanish. It is expected to be small at low energy, where s-wave scattering dominates (k_F R ≪ 1), and to grow at optical energies ~1–3 eV, where the electron wavelength approaches the cluster size. Worth estimating once, e.g. by a first ladder term, or by comparing to Yuan/Cysne exact numerics on an NN-TB vacancy test case.
- **Where non-self-consistent Σ = cT̄ fails.** Near E_D for vacancies (a resonance exactly at or near E_D), the non-self-consistent Σ = cT̄_k diverges as ~c/(|E| ln²). My own estimate of the self-consistent crossover scale, from the SCTMA equation at ε = 0, is Γ_c ~ Δ·√(c/ln(1/c)): about 0.1 eV for c = 0.1% and Δ ~ 7 eV. This agrees in order of magnitude with Cysne's μ_c ≈ 0.15 eV at 0.4% and with Peres's √n_i scaling. For |μ|, ħω ≲ Γ_c, the user's non-self-consistent bubble is expected to be quantitatively unreliable: the IR background, the ħω ≈ μ peak and the Drude width at low doping.
- **Where the scheme should be fine.** For |μ| and ħω well above Γ_c (e.g. μ ≳ 0.2–0.3 eV at c ≲ 0.1%, and the visible range), Σ = cT̄ is controlled at first order in c.
- **Pitfall: magnetism.** A DFT vacancy is magnetic in spin-polarized DFT. The resonance position, and hence where the ħω ≈ μ feature and the minimal conductivity sit (Skrypnyk–Loktev: pinned to the resonance energy), depends on nspin 1 vs 2. This follows from Skrypnyk–Loktev's resonance-energy result, not from a source on DFT vacancies.

### Gaps
- No paper found that computes the ladder vertex correction to the **optical** (interband, finite-ω) conductivity for resonant scatterers in graphene, nor one that quantifies the bubble-vs-exact discrepancy for σ(ω) with vacancies.
- The exact prefactor and log power of Ostrovsky's Γ_η could not be read unambiguously from the PDF text extraction; the √(c/ln(1/c)) form above is my own derivation and should be checked against Eq. (56) of the paper.
- The classic single-site CPA result that vertex corrections vanish (usually attributed to Velický, Phys. Rev. 184, 614 (1969)) was not verified in this session: UNVERIFIED.

## Q4. Has anyone combined a DFT-derived (not fitted) vacancy T-matrix with a Kubo σ(ω)?

### Takeaway
I found no such paper for graphene. The closest pieces are:
- DFT-**fitted** TB adatom parameters (Wehling 2010: V ≈ 2t, ε_d ≈ −t/16) used in exact real-space Kubo σ(ω) (Yuan 2011) and ε(q,ω) (Yuan 2012);
- a first-principles DFT **T-matrix** for vacancies in graphene and TMDs, giving DOS, linewidths and Boltzmann dc (Kaasbjerg 2020), with no σ(ω);
- ordered-supercell DFT σ(ω) for hydrogenated graphene (Putz, Gmitra, Fabian 2014), with no disorder averaging;
- a first-principles T-matrix for **exciton**–defect optics in MoS₂ (Chan et al., Nano Lett. 2026).

The user's combination (DFT vacancy T-matrix → Σ = cT̄ → Kubo bubble σ(ω) with the full-Berry-connection Wannier velocity) appears new in the literature I could find.

### Cited Findings
- **Wehling et al. 2010 (DFT → fitted TB, dc only).** DFT supercells for CH₃, C₂H₅, CH₂OH, H, OH, mapped onto a non-interacting Anderson impurity on NN-TB graphene (t ≈ 2.6 eV). "The band structure of graphene with a methyl group is well fitted with V ≈ 2t = 5.2 eV and ε_d ≈ −t/16 = −0.16 eV". Other groups have |V| ≳ 2t, |ε_d| ≲ 0.1t. Conductivity by Boltzmann theory and numerically exact Kubo (dc). — [arXiv:1003.0609](https://arxiv.org/abs/1003.0609); [DOI](https://doi.org/10.1103/PhysRevLett.105.056802)
- **Yuan et al. 2011 (DFT-fitted parameters → exact Kubo σ(ω)).** Uses "V ≈ 2t and ε_d ≈ −t/16 … obtained from the ab initio density functional theory (DFT) calculations" for H/organic adatoms. Vacancies are modeled by simple site removal (no DFT). — [arXiv:1109.3485](https://arxiv.org/abs/1109.3485)
- **Kaasbjerg 2020, "Atomistic T-matrix theory of disordered two-dimensional materials", PRB 101, 045433.**
  - Built on "realistic density-functional theory (DFT) descriptions of the defects and their scattering matrix elements".
  - Computes the DOS with bound states, the quasiparticle spectrum, the linewidth/scattering rate, and conductivity/mobility via the Boltzmann equation with a T-matrix collision integral, for graphene and TMDs with vacancies and substitutional atoms.
  - Treats the self-energy at the Born and non-self-consistent T-matrix ("full Born") level.
  - Argues the T-matrix route beats Kubo approaches for dc at meV resolution in dilute systems.
  - Lists extensions (spin-dependent potentials, charged defects needing self-consistency) but **no optical conductivity**. A search of the full text for "optical" returns only introductory mentions.
  — [arXiv:1911.00530](https://arxiv.org/abs/1911.00530); [DOI](https://doi.org/10.1103/PhysRevB.101.045433)
- **Putz, Gmitra, Fabian, "Optical conductivity of hydrogenated graphene from first principles", PRB 89, 035437 (2014).** WIEN2k LAPW DFT on graphene supercells of various sizes, each with one H atom (uniform, ordered coverage). The complex σ(ω) in the IR, visible and UV "has different characteristic features depending on the degree of hydrogen coverage". This is ordered, periodic-defect σ(ω), not a disorder-averaged self-energy. — [arXiv:1309.1016](https://arxiv.org/abs/1309.1016); [DOI](https://doi.org/10.1103/PhysRevB.89.035437)
- **Chan, Haber, Naik, Qiu, da Jornada, "Exciton-Defect Interaction and Optical Properties from a First-Principles T-Matrix Approach", Nano Lett. 26, 961–966 (2026).** First-principles T-matrix for exciton–defect scattering in monolayer MoS₂. Disorder-averaged exciton Green's function; absorption and photoluminescence; bound exciton–defect states captured. This is the methodological analogue for optics, but excitonic (BSE) and for a TMD, not a graphene Kubo bubble. — [arXiv:2505.15523](https://arxiv.org/abs/2505.15523); [DOI](https://doi.org/10.1021/acs.nanolett.5c04479)
- **Lherbier, Dubois, Declerck, Niquet, Roche, Charlier, "Transport properties of graphene containing structural defects", PRB 86, 075402 (2012).** Title only; to my knowledge this is ab initio-parametrized TB for structural defects with real-space Kubo dc, but the content was not verified here. — [DOI](https://doi.org/10.1103/PhysRevB.86.075402)
- **Current operator in model σ(ω) work.** Model works use the TB current J = −(ie/ħ) Σ t_ij (r_j − r_i) c†_i c_j (Peierls, diagonal position). Intra-cell / Berry-connection position matrix elements are therefore absent. — [Yuan 2010](https://arxiv.org/abs/1007.3930)

### Inferences
- Novelty, as far as this search goes:
  - The individual ingredients exist: DFT T-matrix (Kaasbjerg), DFT-fitted resonant-adatom TB plus exact Kubo σ(ω) (Yuan 2011), and ordered-supercell DFT σ(ω) (Putz 2014).
  - I found no paper combining a first-principles, unfitted vacancy T-matrix self-energy with a Kubo σ(ω) for graphene.
  - I also found none using a Wannier velocity with the full Berry connection in a disordered-graphene σ(ω).
  - The defensible claim is "first-principles extension of the Peres/Stauber CPA-bubble and Kaasbjerg T-matrix approaches to the optical conductivity", with Yuan 2011 and Cysne 2016 as model benchmarks.
- Limitations a referee would raise, from Q3:
  - non-self-consistency near E_D;
  - no vertex corrections, which cannot be neglected a priori for an extended multi-orbital cluster T-matrix;
  - no localization physics at the Dirac point.

  An NN-TB vacancy cross-check (T-matrix → bubble vs Yuan/Cysne exact numerics) would address the first two directly.

### Gaps
- No systematic forward-citation search of Kaasbjerg 2020 was done. A later paper applying the atomistic T-matrix to optics (graphene or TMD, electron picture) could exist and was not found.
- Not checked: DFT-parametrized (e.g. Gmitra–Kochan–Fabian H-adatom, PRL 110, 246602 (2013), [DOI](https://doi.org/10.1103/PhysRevLett.110.246602), title verified only) TB models later fed into KPM optical conductivity.
