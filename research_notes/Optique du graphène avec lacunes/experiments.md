# Experiments on the optical conductivity of disordered or defective graphene (residual sub-2E_F absorption, irradiated or functionalised graphene, the van Hove/excitonic peak, links to Raman defect metrics)

Scope: which experiments an ab initio Kubo–Greenwood + T-matrix calculation of σ(ω) for graphene with dilute
vacancies could be compared against. The calculation would predict (i) absorption inside the Pauli-blocked
window 0 < ħω < 2|μ|, (ii) broadening of the 2|μ| edge, and (iii) changes near the M-point (van Hove) peak.

Verification legend, used on every bibliographic entry:
- **[V-full]**: I read the arXiv full text and extracted numbers from it. Bibliographic data were checked against Crossref, OpenAlex or the arXiv metadata.
- **[V-abs]**: I checked the bibliographic data (Crossref, OpenAlex or arXiv API) and read the abstract only.
- **[V-meta]**: I checked the bibliographic data only.
- **UNVERIFIED**: I could not check this against a source.

Units: σ0 = πe²/2h = e²/4ħ ≈ 6.08×10⁻⁵ S is the universal sheet conductivity, and the absorbance of freestanding
graphene is πα ≈ 2.3 %. 1 cm⁻¹ = 0.12398 meV. All conversions from cm⁻¹ to meV below are mine.

---

## Q1. Residual (sub-2E_F) infrared absorption in gated graphene: magnitudes, proposed explanations, and whether vacancies were considered

### Takeaway
Every gated IR experiment on single-layer graphene (exfoliated or CVD) finds that Pauli blocking below 2E_F is
incomplete. The residual σ1 is of order 0.1–0.3 σ0 and depends only weakly on gate voltage. The 2E_F edge is
far broader than thermal smearing (~1400 cm⁻¹ ≈ 170 meV in Li et al. 2008). By the f-sum rule, the Drude weight
is reduced by ~20–45 %. The authors propose charged impurities, "unitary scatterers (edge defects, cracks,
vacancies)", electron–phonon coupling, electron–electron interactions, localisation and doping
inhomogeneity. Li et al. explicitly tested the one disorder theory that included vacancy-like unitary scatterers
and found it too weak and with the wrong gate dependence. None of these experiments varied the vacancy density.

### Cited Findings

**Li, Henriksen, Jiang, Hao, Martin, Kim, Stormer, Basov, "Dirac charge dynamics in graphene by infrared spectroscopy", Nature Physics 4, 532–535 (2008), DOI 10.1038/nphys989, arXiv:0807.3780 [V-full; arXiv journal-ref "Nature Physics 4, 532 (2008)"; DOI, pages 532–535 and author list confirmed via OpenAlex]**
- Samples: exfoliated monolayer (Kish graphite) on 300 nm SiO2/Si back gate. Mobility up to 8700 cm² V⁻¹ s⁻¹ at ~2×10¹² cm⁻². Synchrotron IR micro-spectroscopy (ALS) at 700–8000 cm⁻¹ (≈ 87–990 meV), with R(ω) and T(ω) measured versus V_g at 45 K — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- At the charge-neutrality voltage, σ1 matches the universal value πe²/2h to within ±15 % over 4000–6500 cm⁻¹ — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- The 2E_F threshold in σ1(ω,V) has a width of **~1400 cm⁻¹ (~174 meV)**. The width does not depend on gate voltage or carrier density, even though N changes seven-fold between 10 V and 71 V. The thermal estimate at 45 K is only ~500 cm⁻¹ (~62 meV). The authors cite theory saying that "disorder effects and electron-phonon coupling are needed" to explain the width. They also invoke spatial inhomogeneity of E_F (electron–hole puddles), averaged over the few-µm IR spot — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- The authors "registered significant conductivity below 2E_F", which "has not been anticipated by theories developed for Dirac Fermions". The Drude scattering rate is 1/τ = 30 cm⁻¹ (~3.7 meV) at 71 V (from transport), far below 2E_F, so the Drude tail cannot explain it — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- **Magnitude (supplement):** "An assumption of the universal value for σ1(ω,V_CN) implies that the residual conductivity is as strong as **0.3·πe²/2h**." The value depends on the σ1(V_CN) reference below 4000 cm⁻¹. Other experiments support σ1(V_CN) = σ0 across the mid-IR, which "implies strong residual absorption below the 2E_F cut-off that is **nearly independent of the applied voltage**" — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- **Vacancies were explicitly considered and judged insufficient.** Quote: "charged impurities and unitary scatterers (edge defects, cracks, vacancies, etc) were shown to induce considerable residual conductivity below 2E_F [Ref. 11]. However, the theoretical residual absorption in Ref. [11] is systematically suppressed with voltage, whereas this suppression was not observed in our data. In addition, the magnitude of the theoretical residual absorption is smaller compared to experimental values… it is likely that other mechanisms are also responsible." The authors favour many-body effects: a frequency-dependent 1/τ(ω) from electron–electron and electron–phonon interactions — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- The Fermi velocity rises from v_F ≈ 1.10×10⁶ m/s at high doping to ~1.25×10⁶ m/s at low bias. The authors interpret this as a many-body renormalisation — [arXiv:0807.3780](https://arxiv.org/abs/0807.3780)
- Citation mismatch: Li's Ref. 11 is printed as "Peres, Stauber & Castro Neto, PRB 78, 085418 (2008), preprint arXiv:0803.2816". These are two different papers:
  - arXiv:0803.2816 is Peres, Stauber, Castro Neto, "The infrared conductivity of graphene", EPL 84, 38002 (2008), DOI 10.1209/0295-5075/84/38002 — [arXiv abs page](https://arxiv.org/abs/0803.2816)
  - PRB 78, 085418 is Stauber, Peres, Castro Neto, "Conductivity of suspended and non-suspended graphene at finite gate voltage" — [Crossref DOI 10.1103/PhysRevB.78.085418](https://doi.org/10.1103/PhysRevB.78.085418)
  - The Li text (vacancies, the gate dependence of the residual) matches the EPL preprint.
- The EPL/arXiv:0803.2816 theory that Li compared against has the following ingredients. Unitary scatterers "due to structural disorder (edge defects, cracks, vacancies)" at **n_i = 4.0×10⁻⁵ per carbon**. Charged impurities at n_i^C = 1.3×10⁻⁴ per carbon. Phonons, treated with RPA + CPA. It gives "an anomalous enhancement of the conductivity in a frequency region that is blocked by Pauli exclusion and an impurity broadening of the conductivity threshold". In their fits the edge broadening is controlled by the charged impurities, not by the unitary scatterers — [arXiv:0803.2816](https://arxiv.org/abs/0803.2816) [V-full]

**Wang, Zhang, Tian, Girit, Zettl, Crommie, Shen, "Gate-variable optical transitions in graphene", Science 320, 206–209 (2008), DOI 10.1126/science.1152793 [V-abs]**
- First demonstration of gate-tunable interband absorption (Pauli blocking) in monolayer and bilayer graphene by IR spectroscopy. For the monolayer it gives the linear Dirac dispersion directly — [Crossref/OpenAlex DOI](https://doi.org/10.1126/science.1152793)
- I did not access the full text, so I did not extract a residual-absorption magnitude from this paper (see Gaps).

**Mak, Sfeir, Wu, Lui, Misewich, Heinz, "Measurement of the optical conductivity of graphene", Phys. Rev. Lett. 101, 196405 (2008), DOI 10.1103/PhysRevLett.101.196405, arXiv:0810.1269 [V-full]**
- Exfoliated monolayers on SiO2, reflectance and transmission over 0.2–1.2 eV at room temperature. Over **0.5–1.2 eV the absorbance is flat at A = (2.3 ± 0.2) % = (1.0 ± 0.1) πα**, i.e. σ = σ0 — [arXiv:0810.1269](https://arxiv.org/abs/0810.1269)
- Below ~0.5 eV there is "a significant deviation from the value of the universal absorbance" and the results "differ from sample to sample". The absorbance of sample 1 drops strongly below 0.4 eV; sample 2 shows only a modest decrease. The authors attribute this to different spontaneous doping, finite temperature and intraband (Drude) transitions — [arXiv:0810.1269](https://arxiv.org/abs/0810.1269)
- The authors list many-body effects and non-linear band dispersion as limits of the free-Dirac model. Because excited-state lifetimes are very short, "significant spectral broadening may be present, leading to enhanced absorption even at lower photon energies" — [arXiv:0810.1269](https://arxiv.org/abs/0810.1269)

**Horng, Chen, Geng, Girit, Zhang, Hao, Bechtel, Martin, Zettl, Crommie, Shen, Wang, "Drude conductivity of Dirac fermions in graphene", Phys. Rev. B 83, 165113 (2011), DOI 10.1103/PhysRevB.83.165113, arXiv:1007.4623 [V-full]**
- Samples: large-area CVD graphene on SiO2/Si, back-gated. Hole mobility ~2700 cm²/V·s; initial hole doping 1.05×10¹² cm⁻² at V_g = 0; CNP at 14 V. Transmission measured over **30–6000 cm⁻¹ (≈ 3.7–744 meV)**, i.e. THz to mid-IR — [arXiv:1007.4623](https://arxiv.org/abs/1007.4623)
- The gate-induced Drude absorption reaches >15 % at THz frequencies. The Drude fits reproduce the DC conductivity measured directly. The scattering rate γ is roughly constant for holes but larger and density-dependent for electrons. This "cannot be described by pure unitary scattering or charge impurity scattering of Dirac fermions. A combination of different scattering mechanisms seems to be necessary" — [arXiv:1007.4623](https://arxiv.org/abs/1007.4623)
- **Drude weight:** in the strongly doped region (|V_g − V_cnp| > 20 V) the measured D is **"lower than the theoretical value [D = v_F e² √(πn)/ħ-type Boltzmann value] by 20–45 %"**, with experimental uncertainty <10 %. D is also electron–hole asymmetric — [arXiv:1007.4623](https://arxiv.org/abs/1007.4623)
- **Pauli window:** "the observed interband optical conductivity decrease in doped graphene is appreciably less than σ0". Equivalently, residual absorption persists below 2E_F, and the deviation is larger for electron doping. The integrated interband loss equals the integrated Drude gain "almost perfectly, just as the sum rule requires". So the missing Drude weight is the same thing as the residual sub-2E_F absorption. The authors note the same reduction in exfoliated graphene (Li et al.) and conclude that the "anomalous Drude weight reduction is general for both exfoliated and CVD samples" — [arXiv:1007.4623](https://arxiv.org/abs/1007.4623)
- Proposed explanations:
  - Electron–electron interactions: a pseudospin-driven Drude weight reduction (Polini et al.), but that theory predicts ">80 %", far more than observed.
  - Electron–impurity interactions, which could cause the e–h asymmetry.
  - Carrier localisation, which "clearly goes beyond impurity scattering of Dirac fermions considered in current quantum transport theory". Vacancies are not named specifically.
  - [arXiv:1007.4623](https://arxiv.org/abs/1007.4623)

**Yan, Xia, Zhu, Freitag, Dimitrakopoulos, Bol, Tulevski, Avouris, "Infrared spectroscopy of wafer-scale graphene", ACS Nano 5, 9854–9860 (2011), DOI 10.1021/nn203506n, arXiv:1111.3714 [V-full]**
- Samples: epitaxial graphene on SiC, and CVD graphene transferred onto quartz. For CVD on quartz, "with decreasing photon frequency in the mid-IR, the absorption decreases **from 1.9% to 0.4%** due to Pauli blocking". In the far-IR the absorption rises to 25 % at 40 cm⁻¹ (Drude) — [arXiv:1111.3714](https://arxiv.org/abs/1111.3714)
- The Pauli-blocking step is fit by a Gaussian-broadened step, centred near 5500 cm⁻¹, giving |E_F| = 341 meV and n = 7.1×10¹² cm⁻². The FWHM appears as "2001 cm⁻¹" in the extracted arXiv text. This may be a garbled "200 ± 1"; the value to check in the published version is ~2000 cm⁻¹ ≈ 250 meV. Raman G-peak maps give a mean E_F = −321 meV and show **doping inhomogeneity that "has considerable contribution to the broadening of the absorption step at 2|E_F|"** — [arXiv:1111.3714](https://arxiv.org/abs/1111.3714)
- Drude fit (CVD): D = 1.24×10⁴ e²/h cm⁻¹, Γ = 102 cm⁻¹ (~12.6 meV, τ ~50 fs), R_s = 670 Ω/□, hole mobility 1300 cm²/V·s. The measured D is **~28 % smaller** than v_F e²√(πn) computed from the Pauli-blocking E_F. After heavy chemical doping (2E_F ≈ 1.12 eV), D = 1.8×10⁴ versus a theoretical 2.8×10⁴, i.e. **~35 % smaller**, with Γ = 91 cm⁻¹. The epitaxial samples have Γ ~270 cm⁻¹ (~33 meV, τ ~20 fs) — [arXiv:1111.3714](https://arxiv.org/abs/1111.3714)
- "The reduction of Drude weight varies from sample to sample and **some samples have no reduction at all**." The mechanisms are called "not clear at this stage" — [arXiv:1111.3714](https://arxiv.org/abs/1111.3714)

**Mak, Ju, Wang, Heinz, "Optical spectroscopy of graphene: from the far infrared to the ultraviolet", Solid State Commun. 152, 1341–1349 (2012), DOI 10.1016/j.ssc.2012.04.064 [V-full via author PDF]**
- Review by the Heinz and Wang groups. Near charge neutrality the Drude scattering rate is ~100 cm⁻¹ (momentum lifetime ~50 fs). D_intra is "somewhat lower" than D_inter = (e²/ħ)E_F. "The measured reduction in D_intra is, therefore, correlated to the **imperfect Pauli blocking of interband transitions at energies below 2E_F**. Modification of the free-carrier Drude conductivity can arise from **impurity and defect states** or from many-body interactions." Following Yan et al., "for samples with nearly perfect Pauli blocking… D_intra appears to be only slightly reduced" — [Stanford author PDF](https://web.stanford.edu/group/heinz/publications/Pub193.pdf)
- The review quotes a mean sheet absorbance of A = (2.28 ± 0.14) % over 0.5–1.2 eV, compared with πα = 2.293 % — [Stanford author PDF](https://web.stanford.edu/group/heinz/publications/Pub193.pdf)

**Frenzel, Lui, Shin, Kong, Gedik, "Semiconducting-to-metallic photoconductivity crossover and temperature-dependent Drude weight in graphene", Phys. Rev. Lett. 113, 056602 (2014), DOI 10.1103/PhysRevLett.113.056602 [V-abs]**
- Gated CVD graphene, optical-pump THz-probe. Photoconductivity is positive near zero density and turns negative at high density. This is explained by photoinduced changes of both the Drude weight and the scattering rate. The fluence dependence is non-monotonic, showing that the **Drude weight depends non-monotonically on (electron) temperature** — [DOI](https://doi.org/10.1103/PhysRevLett.113.056602)
- This is not a residual-absorption measurement. It matters because the equilibrium D (and hence, by the sum rule, the sub-2E_F weight) must be compared at the correct carrier temperature.

**Theoretical explanations quoted by the experiments (for context; theory is covered elsewhere)**
- Scharf, Perebeinos, Fabian, Avouris, PRB 87, 035414 (2013), DOI 10.1103/PhysRevB.87.035414 [V-abs]. In a Kubo calculation, phonon-assisted (intrinsic optical plus surface-polar) absorption "in the optical gap… can be as large as **20–25 % of the universal ac conductivity** for graphene on polar substrates (Al2O3, HfO2)" at room temperature. It increases with T, and the Drude weight decreases with T — [DOI](https://doi.org/10.1103/PhysRevB.87.035414)

### Inferences
- The experimental sub-2E_F residual is ~0.2–0.3 σ0. Li's upper estimate is 0.3 σ0. For Yan's CVD sample, 0.4 %/1.9 % ≈ 0.21 of the high-energy plateau (my ratio of their numbers). By the sum rule, Horng's 20–45 % Drude weight deficit corresponds to a comparable fraction of the blocked interband weight.
- A vacancy-only T-matrix theory has to reach these magnitudes at vacancy densities that Raman allows (see Q3), which is probably impossible for high-quality exfoliated samples. So the dilute-vacancy calculation will probably not "explain" the gated-graphene residual. It is better framed as a prediction for irradiated samples, where the vacancy contribution can be separated by its linear scaling with dose.
- Two signatures could discriminate vacancies from other mechanisms:
  - The gate dependence of the residual. Li finds it nearly V-independent, whereas the unitary-scatterer theory gives a residual that falls with V.
  - The width of the 2E_F edge. Li finds it independent of density, ~170 meV. Doping inhomogeneity (puddles, Yan's Raman maps) is a competing, non-intrinsic broadening that any comparison must subtract or convolve.
- Any comparison with experiment must also include substrate phonons (SiO2 features near 1000–1200 cm⁻¹ appear in both Li and Horng) and the 20–25 % phonon-assisted background predicted on polar substrates.

### Gaps
- I did not read the full text of Wang et al. (Science 2008), so I have no residual magnitude from that paper.
- The published Nature Physics version of Li et al. differs from the arXiv v1 I read ("The published version is modified", per the arXiv comment). The 0.3 σ0 figure comes from the arXiv supplementary text.
- I found no gated IR study that measures the residual sub-2E_F absorption against a **controlled vacancy density** at fixed E_F.
- "Frenzel et al." in the request could also mean Frenzel et al., APL 102, 113111 (2013), "Observation of suppressed terahertz absorption in photoexcited graphene" (arXiv:1301.6108). That one is **UNVERIFIED**: only the arXiv title was seen in search results.

---

## Q2. Optical, IR and THz measurements on intentionally defected graphene (ion or electron irradiation, plasma, hydrogenation, fluorination), including the van Hove/excitonic peak

### Takeaway
Experiments that measure σ(ω) in the IR or UV of graphene with a *controlled* defect density are surprisingly
scarce.
- **Hydrogenation:** one broadband transmission study (FIR to UV) shows the Dirac-band interband and M-point absorption shrinking with H coverage, and non-Drude far-IR response (localisation).
- **Heavy irradiation** (MeV Cu ions, long e-beam exposure) lowers the optical absorption and the THz conductivity, but only qualitative dose dependence is available.
- **Vacancies:** the best-quantified dose dependence for (probably vacancy-type) defects is in DC transport (Chen et al. 2009: 1/µ ∝ dose, midgap-state scattering).
- **UV peak:** the pristine reference is well established (4.62 eV peak, Fano exciton ~0.6 eV below the GW band-to-band M-point energy, Γ ≈ 0.78 eV). I found no experiment that follows its damping or shift against a Raman-calibrated vacancy density. The only defect-dependent UV data are for hydrogenation (intensity decrease), and the controlled UV experiment available is doping (shift and broadening).

### Cited Findings

**UV/visible reference: the M-point (saddle-point) exciton peak in pristine graphene**
- Mak, Shan, Heinz, "Seeing many-body effects in single- and few-layer graphene: observation of two-dimensional saddle-point excitons", PRL 106, 046401 (2011), DOI 10.1103/PhysRevLett.106.046401, arXiv:1012.2922 [V-full]:
  - σ(E) measured over 0.2–5.3 eV. It is universal (σ0) over 0.5–1.5 eV, rises by **~80 % at 3.0 eV**, and peaks at **E_exp = 4.62 eV**. The peak is asymmetric and "red-shifted by nearly 600 meV" from the GW band-to-band M-point energy E_GW = 5.20 eV — [arXiv:1012.2922](https://arxiv.org/abs/1012.2922)
  - Fano fit (continuum = GW σ convolved with **250 meV** broadening): q = −1, E_res = 5.02 eV (180 meV below E_GW), **Γ = 780 meV**, "an exciton lifetime of only ~0.5 fs". The resonance depends only weakly on layer number — [arXiv:1012.2922](https://arxiv.org/abs/1012.2922)
- Kravets, Grigorenko, Nair, Blake, Anissimova, Novoselov, Geim, "Spectroscopic ellipsometry of graphene and an exciton-shifted van Hove peak in absorption", PRB 81, 155413 (2010), DOI 10.1103/PhysRevB.81.155413, arXiv:1003.2618 [V-full]:
  - Exfoliated flakes on SiO2/Si and on amorphous quartz, measured by spectroscopic ellipsometry. The absorption is 2.3 % in the visible and shows a "pronounced peak reaching ~10 % in ultraviolet" at **4.6 eV**, "downshifted by 0.5 eV probably due to excitonic effects" — [arXiv:1003.2618](https://arxiv.org/abs/1003.2618)
  - There is sample-to-sample variation. "Most (~75%) of our samples very closely followed the lower curve", which deviates little from 2.3 % in the visible. "The sample dependence of the onset of the ultraviolet peak remains unclear and **could be due to contamination or rippling**." Yet "dozens of the measured samples exhibited the absorption peak precisely at the same energy" — [arXiv:1003.2618](https://arxiv.org/abs/1003.2618)
- Chae, Utikal, Weisenburger, Giessen, v. Klitzing, Lippitz, Smet, "Excitonic Fano resonance in free-standing graphene", Nano Lett. 11, 1379–1382 (2011), DOI 10.1021/nl200040q [V-abs]. On suspended monolayers and bilayers, measured by transmission over 1.5–5.5 eV, the UV is "dominated by an asymmetric Fano resonance". The Fano model "quantitatively describes the experimental data all the way down to the infrared" — [DOI](https://doi.org/10.1021/nl200040q)
- Mak, da Jornada, He, Deslippe, Petrone, Hone, Shan, Louie, Heinz, "Tuning many-body interactions in graphene: the effects of doping on excitons and carrier lifetimes", PRL 112, 207401 (2014), DOI 10.1103/PhysRevLett.112.207401 [V-abs]. Electron and hole doping cause "shift, broadening, and modification in shape of the saddle-point exciton resonance". These reflect carrier screening and changes in quasiparticle lifetimes, and are reproduced by GW-BSE — [DOI](https://doi.org/10.1103/PhysRevLett.112.207401)

**Hydrogenated graphene (sp³ adsorbates, which are resonant scatterers like vacancies), FIR to UV**
- Lee, Leconte, Kim, Cho, Lyo, Choi, "Optical spectroscopy study on the effect of hydrogen adsorption on graphene", Carbon 103, 109–114 (2016), DOI 10.1016/j.carbon.2016.03.008 [V-abs, full abstract read]:
  - Transmission spectroscopy over far-IR to UV. "For low hydrogen concentration, the absorption intensities of the interband transitions occurring in the Dirac band (mid-IR and visible) and the **M-point van Hove singularity (UV) decrease with increasing hydrogen coverage**." The change is "quantified successfully using the effective medium theory".
  - At the highest coverage, a band gap of >6 eV opens.
  - "The optical conductivity in the Far-IR regime is behaving in a **non-Drude type manner** along with the hydrogenation, implying H-induced localization of the free Dirac π electrons."
  - Source: [UOS repository record](https://pure.uos.ac.kr/en/publications/optical-spectroscopy-study-on-the-effect-of-hydrogen-adsorption-o/)
- Nair, Ren, Jalil, Riaz, Kravets, Britnell, Blake, Schedin, Mayorov, Yuan, Katsnelson, Cheng et al., "Fluorographene: a two-dimensional counterpart of Teflon", Small 6, 2877–2884 (2010), DOI 10.1002/smll.201001555 [V-abs]. This is the fully functionalised limit (CF): an insulator (>10¹² Ω) with an **optical gap of 3 eV** — [DOI](https://doi.org/10.1002/smll.201001555)

**Ion or electron irradiation: optical and THz data**
- Jahanzaib, Jalil, Aisida, Zhao, Dee, Sorokin, Ahmad, Ul-Hamid, "Cu ions irradiation-induced defects in graphene and their effects on optical properties", Radiat. Phys. Chem. 193, 110008 (2022), DOI 10.1016/j.radphyschem.2022.110008 [V-abs, partial]. 8 MeV Cu²⁺ at 1×10¹⁵ to 1×10¹⁶ ions/cm², characterised by SEM, TEM and Raman. With increasing fluence the metallic character changes to semiconducting. At the highest fluence the lattice is converted "into the massively disordered layer with **decreased optical absorption**". No quantitative σ(ω) was available in the abstract — [KFUPM record](https://pure.kfupm.edu.sa/en/publications/cu-ions-irradiation-induced-defects-in-graphene-and-their-effects/)
- Feng, Hu, Zhang, Gong, Zhou, Zhong, Liu, Wu, Zhao, Zhang, Liu, "Decrease in terahertz conductivity of graphene under electron beam irradiations", J. Infrared Millim. Terahertz Waves 40, 297–305 (2019), DOI 10.1007/s10762-018-0559-2 [V-meta]. Per the repository snippet: "After a long-time electron beam irradiation, the graphene conductivity decreases, indicating that defects and damages are created in the graphene." I could not access the abstract or full text (publisher-elided and 403 errors) — [UOW record (search snippet)](https://ro.uow.edu.au/aiimpapers/3476)

**Defect dose dependence in the DC limit (the ω → 0 end of σ(ω)), in the same samples as Raman D-band data**
- Chen, Cullen, Jang, Fuhrer, Williams, "Defect scattering in graphene", PRL 102, 236805 (2009), DOI 10.1103/PhysRevLett.102.236805, arXiv:0903.2602 [V-full]:
  - Exfoliated graphene on SiO2, irradiated in UHV by **500 eV Ne⁺ and He⁺**. Each ion produces "one atomic-scale defect, most likely a carbon vacancy with a trapped rare-gas molecule". A dose of 10¹² cm⁻² is ~1 Ne⁺ per 4×10³ C. A D band appears after irradiation — [arXiv:0903.2602](https://arxiv.org/abs/0903.2602)
  - σ(n) is linear in n, and **1/µ increases linearly with dose**: 7.9×10⁻¹⁶ V s (Ne⁺) and 9.3×10⁻¹⁶ V s (He⁺). With the midgap-state (resonant) formula this gives a defect radius R = 2.3 Å (Ne⁺) and 2.9 Å (He⁺), "a reasonable value for single-carbon vacancies". The mobility drop is **4× larger than for the same concentration of singly charged impurities**, and the minimum conductivity falls below 4e²/πh. Defected samples become insulating at low T — [arXiv:0903.2602](https://arxiv.org/abs/0903.2602)
  - Raman caveat from 2009: the Tuinstra–Koenig-type formula gave "L_a ~ 60 nm, larger than the expected defect spacing of 10 nm", closer to the transport mean free path of ~50 nm. Old nanographite calibrations are unsuitable for point defects (see Q3) — [arXiv:0903.2602](https://arxiv.org/abs/0903.2602)
- Shlimak et al., "Raman scattering and electrical resistance of highly disordered graphene", arXiv:1410.3425 [arXiv full text seen; journal reference **UNVERIFIED**]. CVD monolayer samples irradiated by C⁺ up to 10¹⁵ cm⁻², with Raman and I–V on the same samples. At the highest disorder the Raman lines vanish and the resistance rises exponentially. The maximal resistance of a continuous film is ~πh/4e² ≈ 20 kΩ — [arXiv:1410.3425](https://arxiv.org/abs/1410.3425)
- Shlimak et al., "Irradiation-induced metal-insulator transition in monolayer graphene", arXiv:1907.10472 [arXiv seen; journal reference **UNVERIFIED**]. Ion-irradiated (C⁺, Xe⁺) CVD graphene with the disorder monitored by Raman I_D/I_G. Transport goes from metallic to variable-range hopping — [arXiv:1907.10472](https://arxiv.org/abs/1907.10472)

**Non-equilibrium optics on controlled defects (same samples, defect density varied)**
- Alencar, Silva, Malard, de Paula, "Defect-induced supercollision cooling of photoexcited carriers in graphene", Nano Lett. 14, 5621–5624 (2014), DOI 10.1021/nl502163d [V-abs]. Defects were "optically generated in a controlled manner". The pump–probe transient-absorption **decay time decreases as defect density increases**, as supercollision cooling predicts. I did not access how the density was calibrated (see Gaps) — [DOI](https://doi.org/10.1021/nl502163d)

### Inferences
- Lee et al. 2016 is the closest experimental analogue to "σ(ω) of graphene with dilute resonant scatterers across IR–UV". They find that the M-point peak and the Dirac interband absorption *decrease* with coverage. A dilute-vacancy T-matrix calculation should check whether it reproduces this loss of spectral weight, not only broadening.
- Their non-Drude far-IR response (localisation) signals that a Boltzmann/Drude-like dilute-limit treatment will fail at higher defect densities. Chen et al. also see insulating behaviour at low T.
- Chen et al. 2009 gives a quantitative, same-sample dose dependence for vacancy-type defects (1/µ ∝ n_d, R ≈ 2–3 Å, midgap scattering). A T-matrix σ(ω → 0) should reproduce it and serves as the natural cross-check at ω → 0: the same T-matrix fixes the Drude scattering rate and the Pauli-window absorption.
- Doping itself shifts and broadens the 4.6 eV exciton (Mak et al. 2014). Irradiation also dopes graphene (charged/adsorbate defects), so any experimental UV change after irradiation mixes doping, defects and possibly contamination (Kravets' caveat). An ab initio DFT (non-BSE) calculation will place the vHs peak near the band-to-band energy, ~0.5–0.6 eV above experiment. Changes are best compared as relative changes (ratios to pristine) rather than absolute positions.

### Gaps
- I found **no experiment reporting IR σ(ω) inside the Pauli window, or the 2E_F edge width, as a function of a Raman-calibrated vacancy density** in gated or chemically doped graphene. This appears to be an open experimental gap (my assessment after the searches; it could exist in less-indexed literature).
- I found no experiment tracking the **4.6 eV peak damping or shift versus controlled vacancy density**. Only hydrogenation (Lee 2016, intensity decrease) and doping (Mak 2014) are documented.
- Full quantitative data for Feng et al. 2019 (THz vs e-beam dose) and Jahanzaib et al. 2022 (optical vs Cu fluence) were not accessible.
- Plasma-treated graphene optics (e.g., oxygen-plasma photoluminescence) and graphene oxide UV absorption were not verified in this pass and are not reported.

---

## Q3. Links to Raman defect metrics (I_D/I_G vs L_D, I_D/I_D′) and same-sample Raman ↔ optical/THz correlations

### Takeaway
The Raman calibration for point defects is quantitative only in the low-density regime (L_D ≥ 10 nm, "stage 1").
- **Cançado et al. 2011:** n_D (cm⁻²) = (1.8 ± 0.5)×10²² λ_L⁻⁴ (I_D/I_G), with λ_L in nm. I_D/I_G is non-monotonic, with a maximum at L_D ~3 nm.
- **Eckmann et al. 2012:** I_D/I_D′ ≈ 7 for vacancy-like defects, versus ≈ 13 for sp³ and ≈ 3.5 for boundaries.
- **Bruna et al. 2014:** I_D/I_G also depends on doping.
- **Same-sample correlations:** the only verified Raman ↔ optical-conductivity correlation is a THz-TDS study of polycrystalline CVD graphene. The THz mean free path (~11 nm) matches the Raman inter-defect distance (~9 ± 3 nm). There is no equivalent for vacancies and IR/UV σ(ω).

### Cited Findings

**Lucchese, Stavale, Martins Ferreira, Vilani, Moutinho, Capaz, Achete, Jorio, "Quantifying ion-induced defects and Raman relaxation length in graphene", Carbon 48, 1592–1597 (2010), DOI 10.1016/j.carbon.2009.12.057 [V-meta: Crossref, OpenAlex. Abstract publisher-elided; content taken from Cançado 2011 and a companion abstract]**
- Introduces the local-activation model for Ar⁺-bombarded graphene. Each impact creates a structurally disordered region of radius r_S, surrounded by an "activated" ring out to r_A where the sp² lattice survives but D-band scattering occurs. r_A − r_S is "the Raman relaxation length". Cançado et al. use **r_S = 1 nm "as determined in Ref. 27 [Lucchese]"** and r_A = 3.1 nm, "in excellent agreement with the values obtained in Refs. 27, 28, 33" — [Cançado 2011, arXiv:1105.0175](https://arxiv.org/abs/1105.0175)
- The same group's companion paper studies "low energy (90 eV) Ar⁺ ion bombardment in graphene samples as a function of the number of layers": Jorio, Lucchese, Stavale, Martins Ferreira, Moutinho, Capaz, Achete, J. Phys.: Condens. Matter 22, 334204 (2010), DOI 10.1088/0953-8984/22/33/334204 [V-meta] — [BVS record (abstract)](https://busqueda.bvsalud.org/portal/resource/pt/mdl-21386494). The ion energy used in Lucchese 2010 itself is not independently verified here (probably the same procedure).
- The fit constants often quoted for Lucchese 2010 (C_A ≈ 4.2, C_S ≈ 0.87 at 514 nm) are **UNVERIFIED**, because I could not access the text.

**Cançado, Jorio, Martins Ferreira, Stavale, Achete, Capaz, Moutinho, Lombardo, Kulmala, Ferrari, "Quantifying defects in graphene via Raman spectroscopy at different excitation energies", Nano Lett. 11, 3190–3196 (2011), DOI 10.1021/nl201432g, arXiv:1105.0175 [V-full; arXiv journal-ref "Nano Letters 11, 3190 (2011)"; pages 3190–3196 confirmed via OpenAlex]**
- Samples: exfoliated single-layer graphene, Ar⁺-bombarded "as for the procedure outlined in Ref. 27 [Lucchese]". Doses run from **10¹¹ Ar⁺/cm² (one defect per 4×10⁴ C) to 10¹⁵ Ar⁺/cm² (one defect per four C)**, giving L_D = 24, 14, 13, 7, 5 and 2 nm. Excitations: E_L = 1.58, 1.96 and 2.41 eV (785, 633, 514.5 nm) — [arXiv:1105.0175](https://arxiv.org/abs/1105.0175)
- Key formulae, valid for **L_D ≥ 10 nm** and visible excitation — [arXiv:1105.0175](https://arxiv.org/abs/1105.0175):
  - L_D² (nm²) = (4.3 ± 1.3)×10³ E_L⁻⁴ (I_D/I_G)⁻¹ = (1.8 ± 0.5)×10⁻⁹ λ_L⁴ (I_D/I_G)⁻¹
  - n_D (cm⁻²) = 10¹⁴/(πL_D²) = (7.3 ± 2.2)×10⁹ E_L⁴ (I_D/I_G) = **(1.8 ± 0.5)×10²² λ_L⁻⁴ (I_D/I_G)**
  - Fit constants: C_A = (160 ± 48) E_L⁻⁴, r_A = 3.1 nm, r_S = 1 nm, C_S = 0.
- "For all excitations, the D to G intensity ratio **reaches a maximum for an inter-defect distance ~3 nm**. Thus, a given ratio could correspond to two different defect densities… The analysis of the G peak width and its dispersion with excitation energy solves this ambiguity." Γ_G always increases with disorder — [arXiv:1105.0175](https://arxiv.org/abs/1105.0175)

**Eckmann, Felten, Mishchenko, Britnell, Krupke, Novoselov, Casiraghi, "Probing the nature of defects in graphene by Raman spectroscopy", Nano Lett. 12, 3925–3930 (2012), DOI 10.1021/nl300901a, arXiv:1207.2058 [V-abs]**
- "The intensity ratio of the D and D′ peak is maximum (**~13**) for sp³-defects, it decreases for **vacancy-like defects (~7)**, and it reaches a minimum for boundaries in graphite (~3.5)" — [DOI](https://doi.org/10.1021/nl300901a); [arXiv:1207.2058](https://arxiv.org/abs/1207.2058)
- The low-defect-density validity range of these ratios, which I believe the paper specifies, was not checked in the full text (see Gaps).

**Bruna, Ott, Ijäs, Yoon, Sassi, Ferrari, "Doping dependence of the Raman spectrum of defected graphene", ACS Nano 8, 7432–7441 (2014), DOI 10.1021/nn502676g [V-abs]**
- With polymer-electrolyte gating up to E_F = 0.7 eV (Hall-calibrated), "for a given number of defects… the intensities of the D and D′ peaks **decrease with increasing doping**". The authors attribute this to doping-dependent electron–electron scattering of the photoexcited carriers and give "a general relation between D peak intensity and defects valid for any doping level" — [DOI](https://doi.org/10.1021/nn502676g)
- As quoted by Whelan et al. 2024 below, the doping-corrected relation is L_D² (nm²) = (1.2 ± 0.3)×10³ (I_D/I_G)⁻¹ E_F^(−0.54 ± 0.04) / E_L⁴. I verified it only as quoted, not in Bruna's text — [Whelan 2024, PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850153/)

**Same-sample Raman ↔ optical/THz or transport correlations**
- Whelan, De Fazio, Pasternak et al., "Mapping nanoscale carrier confinement in polycrystalline graphene by terahertz spectroscopy", Sci. Rep. 14, 3163 (2024), DOI 10.1038/s41598-024-51548-z [V-abs + page summary]:
  - THz-TDS spectra of polycrystalline CVD graphene follow the **Drude–Smith** model, with median τ = 8.76 fs, backscattering c = −0.83 and ℓ_mfp ≈ 11.3 nm.
  - Raman gives I_D/I_G = 0.8 ± 0.2, E_F ≈ −0.4 ± 0.2 eV and **ℓ_D ≈ 9 ± 3 nm** (Bruna-type formula). The two are "consistent", and the scattering is attributed to small-angle grain boundaries.
  - [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850153/)
- Whelan et al., "Case studies of electrical characterisation of graphene by terahertz time-domain spectroscopy", 2D Materials (2021), DOI 10.1088/2053-1583/abdbcb [V-abs; volume and article number not returned by OpenAlex]. THz-TDS used as a non-contact metrology of sheet conductivity, τ, n and µ — [DOI](https://doi.org/10.1088/2053-1583/abdbcb)
- Ni, Ponomarenko, Nair, Yang, Anissimova, Grigorieva, Schedin, Blake, Shen, Hill, Novoselov, Geim, "On resonant scatterers as a factor limiting carrier mobility in graphene", Nano Lett. 10, 3868–3872 (2010), DOI 10.1021/nl101399r [V-abs]. Substrate-supported graphene shows "a previously unnoticed D peak… with intensity of **~1% with respect to the G peak**". By correlating the Raman D intensity with mobility using H adsorbates as mimics, the authors conclude that "if the intervalley scatterers responsible for the D peak are monovalent, their concentration is sufficient to account for the limited mobilities" — [DOI](https://doi.org/10.1021/nl101399r)
- Chen et al. 2009 (Q2) and Shlimak et al. (Q2) also pair Raman D-band data and transport on the same irradiated samples.
- Yan et al. 2011 (Q1) used Raman on the same samples only to cross-check E_F (G-peak shift) and doping inhomogeneity, not the defect density — [arXiv:1111.3714](https://arxiv.org/abs/1111.3714)

### Inferences
These are my calculations from the formulae above: graphene carbon density ≈ 3.82×10¹⁵ cm⁻²; Cançado's n_D formula at λ_L = 514.5 nm gives n_D ≈ 2.57×10¹¹ (I_D/I_G) cm⁻²; at 532 nm, n_D ≈ 2.25×10¹¹ (I_D/I_G).

**Vacancy fraction → Raman observables (514 nm)**

| Vacancy fraction x (per C) | n_D (cm⁻²) | L_D (nm) | I_D/I_G | Regime |
|---|---|---|---|---|
| 7×10⁻⁷ | 2.7×10⁹ | ~110 | ~0.01 | Order of the ~1 % D/G in Ni et al. 2010; formula applies |
| 4×10⁻⁵ (unitary density in Peres–Stauber–Castro Neto) | 1.5×10¹¹ | ~14 | ~0.6 | Clearly visible D peak |
| 10⁻⁴ | 3.8×10¹¹ | ~9 | ~1.5 | Edge of the formula's validity (L_D ≥ 10 nm) |
| 10⁻³ | 3.8×10¹² | ~2.9 | — | Near the I_D/I_G maximum; formula invalid |
| 1/162 (one vacancy in a 9×9 supercell, 162 C) | 2.4×10¹³ | ~1.2 | — | Stage 2; formula invalid |

- A DFT *supercell* concentration is far too high to compare with Raman stage-1 samples. In the T-matrix approach, however, the concentration is a free parameter that enters linearly (dilute limit). The calculation should be reported per unit vacancy density, e.g. Δσ(ω)/σ0 per 10¹¹ cm⁻², so it can be mapped onto I_D/I_G through Cançado's relation.
- The unitary-scatterer density Peres et al. needed (~1.5×10¹¹ cm⁻²) would give a D peak with I_D/I_G ~0.6. Good exfoliated samples typically show a much smaller D peak (Ni et al.: ~1 % of G). This supports Li's conclusion that vacancies alone do not explain the residual absorption in high-quality gated samples. Caveats: Li's own Raman D intensity was not reported in the text I read, and Ni's excitation wavelength was not checked.
- For irradiated samples, a comparison with experiment needs three Raman checks on the same spot:
  - I_D/I_D′ ≈ 7, to confirm vacancy-like rather than sp³ defects;
  - Γ_G, to break the stage-1/stage-2 ambiguity;
  - E_F (G/2D positions), to apply the Bruna doping correction and fix 2|μ| for the Pauli window.

### Gaps
- I did not verify any study that measures IR or UV σ(ω) (Pauli-window residual, 2E_F edge width, 4.6 eV peak) and Raman L_D on the *same* intentionally vacancy-defected samples. The only verified same-sample Raman ↔ optical-conductivity work is the THz study of grain boundaries (Whelan 2024).
- The Lucchese 2010 fit constants and the exact ion energy were not verified from the paper itself.
- I did not check the Eckmann et al. validity range (the defect-density regime where I_D/I_D′ ≈ 7 holds) in the full text.
- I did not check how Alencar et al. 2014 calibrated defect density (probably Raman, as in the UFMG group) — **UNVERIFIED**.

---

## Q4. Magnetic and spin experiments on vacancies (brief; the theory is covered elsewhere)

### Takeaway
Magnetometry on irradiated graphene laminates shows **spin-½ paramagnetism with 0.1–0.4 µ_B per vacancy**, with no
magnetic order down to liquid-He temperatures. STM on graphite vacancies shows a sharp resonance at E_F. Spin transport
shows that vacancies (and H adatoms) form local moments that scatter pure spin currents. So a spin-resolved T-matrix is
physically motivated, but the measured moment per vacancy is well below the idealised 1 µ_B.

### Cited Findings
- Nair, Sepioni, Tsai, Lehtinen, Keinonen, Krasheninnikov, Thomson, Geim, Grigorieva, "Spin-half paramagnetism in graphene induced by point defects", Nature Physics 8, 199–202 (2012), DOI 10.1038/nphys2183, arXiv:1111.3775 [V-full]:
  - Vacancies were made in graphene laminates with 350–400 keV protons and 20 MeV C⁴⁺ ions; vacancy numbers come from SRIM. Fluorine adatoms were also studied.
  - Brillouin fits require **J = S = 1/2**. The moment per vacancy is **0.1–0.4 µ_B**, "in reasonable agreement with the theoretical value of µ_B expected for isolated vacancies" given the SRIM uncertainties and possible reconstruction into non-magnetic divacancies and similar defects.
  - No magnetic ordering was found. The maximum was about one moment per ~1000 C. Samples above ~10²⁰ vacancies g⁻¹ disintegrated.
  - [arXiv:1111.3775](https://arxiv.org/abs/1111.3775)
- Nair, Tsai, Sepioni, Lehtinen, Keinonen, Krasheninnikov, Castro Neto, Katsnelson, Geim, Grigorieva, "Dual origin of defect magnetism in graphene and its reversible switching by molecular doping", Nat. Commun. 4, 2010 (2013), DOI 10.1038/ncomms3010 [V-meta; abstract not retrieved, so content **UNVERIFIED**] — [DOI](https://doi.org/10.1038/ncomms3010)
- Ugeda, Brihuega, Guinea, Gómez-Rodríguez, "Missing atom as a source of carbon magnetism", PRL 104, 096804 (2010), DOI 10.1103/PhysRevLett.104.096804 [V-abs]. STM/STS of artificially created single vacancies on graphite reveals "a sharp electronic resonance at the Fermi energy around each single graphite vacancy". This is associated with local moments and "implies a dramatic reduction of the charge carriers' mobility" — [DOI](https://doi.org/10.1103/PhysRevLett.104.096804)
- McCreary, Swartz, Han, Fabian, Kawakami, "Magnetic moment formation in graphene detected by scattering of pure spin currents", PRL 109, 186604 (2012), DOI 10.1103/PhysRevLett.109.186604 [V-abs]. In UHV spin valves, H adatoms produce a dip in the non-local spin signal versus field and an exchange field seen in Hanle precession. "**Lattice vacancies also demonstrate similar behavior**, indicating that the magnetic moment formation originates from p_z-orbital defects" — [DOI](https://doi.org/10.1103/PhysRevLett.109.186604)

### Inferences
- The measured moment of 0.1–0.4 µ_B per vacancy, and the absence of ordering, suggest treating vacancies in σ(ω) calculations as independent (dilute, uncorrelated) spin-½ scatterers. Spin-resolved and spin-unpolarised T-matrices should be compared as bounds, since the experimental moment lies between 0 and 1 µ_B.
- The resonance at E_F seen by STM is the same midgap state that governs resonant scattering. This links the magnetism, the transport results of Chen et al. 2009, and the predicted Pauli-window absorption, which involves transitions to and from the midgap resonance.

### Gaps
- I did not verify the content of Nair et al. 2013 (the dual π/σ origin and switching by doping).
- I did not look for optical or magneto-optical experiments resolving the spin of vacancy-induced absorption; none were encountered.
