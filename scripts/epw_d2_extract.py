#!/usr/bin/env python
r"""
epw_d2_extract.py -- P7: <D^2_Gamma>, <D^2_K> of the 24k-24q chain by two independent routes, no new EPW run.

Definitions used here (block sums, gauge invariant):
  S_g(q, modes) = sum_{i,j in {pi, pi*}} sum_{nu in modes} |g_ij^nu(k=K, q)|^2      [eV^2]   (EPW prtgkk / ph.x 'prt' |g| column, meV)
  D^2 = |g|^2 * 2 M omega / hbar  with  g = sqrt(hbar / (2 M omega)) <k+q| dV/du_nu |k>   (u = displacement along the mass-scaled,
        unit-normalised mode eigenvector; both atoms are carbon so M = M_C)  ->  D^2 [eV^2/A^2] = |g|^2 [eV^2] * 2 M_C * (hbar omega)[eV] / hbar^2,
        the numerical factor being computed in SI below (no value from memory).
  Route 1 (direct):  S_D(Gamma) = sum over LO+TO and the 2x2 pi/pi* block at k = K, q = Gamma;  S_D(K) = A1' mode, same block, q = K.
                     Also the k+q-side "circle": q along the path with |E_{k+q} - E_D| <= 0.2 eV.
  Route 2 (inversion of gamma, Dirac cone, E_F = E_D, T -> 0):
     EPW: Im Pi_nu(q) = pi sum_k w_k sum_ij |g_ij^nu|^2 [f_ik - f_jk+q] delta(e_jk+q - e_ik - hw),  sum_k w_k = 2 (spin; EPW prints DOS per spin,
          wkf = 2/N_k on the full mesh), i.e. sum_k w_k F(k) -> 2 * A_c * int d^2k/(2pi)^2 F(k).
     Dirac cone e = +-v|k|: interband pi -> pi* only, delta(2 v k - hw) -> ring k_r = hw/(2v),  int d^2k delta(2vk - hw) = 2 pi k_r /(2v) = pi hw/(2 v^2).
     Ring-averaged interband weight of one mode = S_g(mode)/4 [pseudospin algebra: intra = inter after angle average, pi->pi* = pi*->pi].
     E2g (q = 0): both valleys contribute (K->K, K'->K'), per mode nu (LO or TO):
          gamma_nu = pi * 2 * 2 * A_c/(2pi)^2 * pi hw/(2v^2) * S_g^nu/4 = A_c hw S_g^nu /(8 v^2) = A_c hw S_g(LO+TO) /(16 v^2)
     A1' (q = K): only the K -> K' valley pair (k near K' + K = Gamma has no Dirac states), single mode:
          gamma_K = pi * 2 * 1 * A_c/(2pi)^2 * pi hw_K/(2v^2) * S_g(K)/4 = A_c hw_K S_g(K) /(16 v^2)
     Both: S_g = 16 v^2 gamma / (A_c hw), gamma = EPW HWHM per mode (10 K), v = Fermi velocity fitted on OUR bands (|e - E_D| <= 0.2 eV).
  Note the omega cancellation: gamma = A_c/(16 v^2) * hw * S_g = A_c hbar^2 /(32 M v^2) * S_D  -> gamma depends on D^2 and v only.
Spin: EPW weights (2); valley: explicit above (2 at Gamma, 1 at K); the DOS N(E_F) is not used.
Usage: python scripts/epw_d2_extract.py [--root <24k-24q dir>]
"""
import argparse, os, re, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
ap = argparse.ArgumentParser(); ap.add_argument("--root", default="/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/epw/24k-24q")
ap.add_argument("--tag", default="24k24q", help="validation_<tag>.npz and d2_extract_<tag>.npz"); ap.add_argument("--phself-tag", default="path_1200_dg0.02"); ap.add_argument("--E_D", type=float, default=-4.2389); ap.add_argument("--out-tag", default=None, help="d2_extract_<out-tag>.npz (default: --tag)"); a = ap.parse_args()
R = a.root; E_D = a.E_D

# ---------------- constants (SI), conversion factor computed, not quoted
hbar = 1.054571817e-34; eV = 1.602176634e-19; u_kg = 1.66053906660e-27; M_C = 12.011 * u_kg
D2_over_g2_per_eV = 2 * M_C * eV / hbar ** 2 * 1e-20          # [A^-2 per eV of hbar*omega]: D^2[eV^2/A^2] = g^2[eV^2] * hw[eV] * this
alat_bohr = np.linalg.norm([4.0354919061, -2.3298923383]); a_A = alat_bohr * 0.529177210903; A_c = np.sqrt(3) / 2 * a_A ** 2
print(f"[const] D^2/g^2 = {D2_over_g2_per_eV:.2f} A^-2 per eV of hbar*omega ; a = {a_A:.4f} A, A_c = {A_c:.4f} A^2")

# ---------------- Fermi velocity of OUR chain: EPW band.eig (474 k on Gamma-K-M-Gamma, absolute eV) and bands.x (181 k), around K
V = np.load(f"results/epw/validation_{a.tag}.npz", allow_pickle=True)
L_path = (2 / 3 + 1 / 3 + 1 / np.sqrt(3)) * 2 * np.pi / a_A                     # A^-1, Gamma-K-M-Gamma
sK = (2 / 3) / (2 / 3 + 1 / 3 + 1 / np.sqrt(3))
def read_plot(f):
    txt = open(f).read().split("\n"); hdr = re.search(r"nbnd=\s*(\d+),\s*nks=\s*(\d+)", txt[0]); nb, nks = int(hdr.group(1)), int(hdr.group(2)); vals = []; i = 1
    while i < len(txt) and len(vals) < nks:
        if len(txt[i].split()) == 3:
            e = []; j = i + 1
            while len(e) < nb: e += [float(x) for x in txt[j].split()]; j += 1
            vals.append(e); i = j
        else: i += 1
    return np.array(vals)
def vfit(s, E_bands, ED, label, emax):
    out = []
    for n, side in ((3, "pi"), (4, "pi*")):
        e = E_bands[:, n] - ED; m = (np.abs(e) <= emax) & (np.abs(e) > 0.02); k = (s[m] - sK) * L_path
        vv = np.polyfit(np.abs(k), np.abs(e[m]), 1)[0]; out.append(vv)
    print(f"[v_F] {label}: pi {out[0]:.3f}, pi* {out[1]:.3f} eV.A  (fit |E-E_D| in [0.02, {emax}] eV, {m.sum()} points)"); return float(np.mean(out))
qb = np.loadtxt(f"{R}/band_freq_interp/graphene_band.kpt", skiprows=1)[:, :2]
Bm = np.array([[1 / np.sqrt(3.0), -1.0], [1 / np.sqrt(3.0), 1.0]]).T; kc = (Bm @ qb.T).T; s_e = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(kc, axis=0), axis=1))]); s_e /= s_e[-1]
E_e = read_plot(f"{R}/band_freq_interp/band.eig")
v_epw = vfit(s_e, E_e, E_D, "EPW Wannier (band.eig)", 0.3)
v_dft = vfit(V["bands_s"], V["bands_E_dft"], float(V["bands_ED_dft"]), "bands.x (DFT)", 0.5)
v = v_epw

# ---------------- route 1a: EPW prtgkk at k = K (epw_g_K), q on the 474-point path
qk = np.loadtxt(f"{R}/epw_g_K/graphene_band.kpt", skiprows=1)[:, :2]
def qdist(q, ref): d = np.mod(q - np.array(ref) + 0.5, 1.0) - 0.5; return np.linalg.norm(d, axis=1)
iG = int(np.argmin(qdist(qk, (0, 0)))); iK = int(np.argmin(qdist(qk, (2 / 3, 1 / 3))))
g = {}
for ib in (4, 5):
    for jb in (4, 5):
        for nu in range(1, 7): g[(ib, jb, nu)] = np.loadtxt(f"{R}/epw_g_K/epw_g_g{ib}{jb}_{nu}.dat")
om = np.array([g[(4, 4, nu)][:, 5] for nu in range(1, 7)]).T * 1e-3                  # (nq, 6) eV
Ekq = g[(4, 4, 1)][:, 4]; Ek = g[(4, 4, 1)][:, 3]
def Sg(iq, modes): return sum(g[(ib, jb, nu)][iq, 7] ** 2 for ib in (4, 5) for jb in (4, 5) for nu in modes) * 1e-6   # eV^2
print(f"[EPW k=K] E_k(4,5) = {g[(4,4,1)][iK,3]:.4f}/{g[(5,5,1)][iK,3]:.4f} eV ; at q=Gamma (iq {iG+1}) E_k+q(4,5) = {g[(4,4,1)][iG,4]:.4f}/{g[(5,5,1)][iG,4]:.4f}, "
      f"omega = {np.round(om[iG]*1e3,2)} meV ; at q=K (iq {iK+1}) E_k+q = {g[(4,4,1)][iK,4]:.4f}/{g[(5,5,1)][iK,4]:.4f}, omega = {np.round(om[iK]*1e3,2)} meV")
oG = np.argsort(om[iG]); modesG = (int(oG[-2]) + 1, int(oG[-1]) + 1); assert abs(om[iG, modesG[0] - 1] - om[iG, modesG[1] - 1]) * 1e3 < 0.5, f"E2g at Gamma not degenerate: {om[iG]*1e3} meV"
iA1 = int(np.argmax([Sg(iK, (nu,)) for nu in range(1, 7)])) + 1     # A1' at K by CHARACTER (largest pi/pi* block sum), never by index or by a target frequency
print(f"[modes] E2g(q=Gamma) = modes {modesG} ({om[iG, modesG[0]-1]*1e3:.2f} meV, two highest), A1'(q=K) = mode {iA1} ({om[iK, iA1-1]*1e3:.2f} meV; branch 3 at degauss 0.002 Ry, 6 at 0.02) ; omega(K) all modes = {np.round(om[iK]*1e3, 2).tolist()} meV")
SgG = Sg(iG, modesG); SgK = Sg(iK, (iA1,)); wG = om[iG, [m - 1 for m in modesG]].mean(); wK = om[iK, iA1 - 1]
SDG = sum(g[(ib, jb, nu)][iG, 7] ** 2 * 1e-6 * om[iG, nu - 1] for ib in (4, 5) for jb in (4, 5) for nu in modesG) * D2_over_g2_per_eV
SDK = SgK * wK * D2_over_g2_per_eV
print(f"[route 1a EPW] q=Gamma: S_g(LO+TO, 2x2) = {SgG:.5f} eV^2 (per mode {SgG/2:.5f}; |g| rms per entry {np.sqrt(SgG/8)*1e3:.1f} meV), hw = {wG*1e3:.2f} meV -> S_D(Gamma) = {SDG:.2f} eV^2/A^2")
print(f"                q=K    : S_g(A1' mode {iA1}, 2x2) = {SgK:.5f} eV^2 (|g| rms {np.sqrt(SgK/4)*1e3:.1f} meV), hw = {wK*1e3:.2f} meV -> S_D(K) = {SDK:.2f} eV^2/A^2 ; S_D(K)/S_D(Gamma) = {SDK/SDG:.3f}")
# the k+q-side circle: q near Gamma / K with |E_{k+q} - E_D| <= 0.2 eV (bands 4/5 at k+q)
def strongest_mode(iq): return int(np.argmax([Sg(iq, (nu,)) for nu in range(1, 7)])) + 1        # A1' branch reorders away from K: take the branch carrying the coupling
for name, iq0, modes_of in (("Gamma", iG, lambda iq: modesG), ("K", iK, lambda iq: (strongest_mode(iq),))):
    near = [iq for iq in range(len(qk)) if abs(g[(4, 4, 1)][iq, 4] - E_D) <= 0.2 and abs(g[(5, 5, 1)][iq, 4] - E_D) <= 0.2 and qdist(qk[[iq]], qk[iq0])[0] < 0.05]
    vals = [sum(g[(ib, jb, nu)][iq, 7] ** 2 * 1e-6 * om[iq, nu - 1] for ib in (4, 5) for jb in (4, 5) for nu in modes_of(iq)) * D2_over_g2_per_eV for iq in near]
    print(f"                circle {name}: {len(near)} q with |E_k+q - E_D| <= 0.2 eV: S_D mean {np.mean(vals):.2f}, min {np.min(vals):.2f}, max {np.max(vals):.2f} eV^2/A^2 (at the point: {vals[near.index(iq0)]:.2f})")

# ---------------- route 1b: ph.x 'prt' DFPT at k = K (dfpt_g_K/ph.out), q = Gamma (block 1) and q = K (block 11); bands 4,5 (E = E_D)
blocks = []; cur = None
for line in open(f"{R}/dfpt_g_K/ph.out"):
    if "q coord.:" in line: cur = {"q": [float(x) for x in line.split(":")[1].split()], "rows": []}; blocks.append(cur)
    p = line.split()
    if cur is not None and len(p) == 10 and p[0].isdigit() and p[1].isdigit() and p[2].isdigit():
        try: cur["rows"].append([float(x) for x in p])
        except ValueError: pass
def dfpt_S(b, modes):
    rows = np.array(b["rows"]); sel = np.isin(rows[:, 0], (4, 5)) & np.isin(rows[:, 1], (4, 5)) & np.isin(rows[:, 2], modes)
    r = rows[sel]; Sg_ = (r[:, 7] ** 2).sum() * 1e-6; SD_ = ((r[:, 7] ** 2) * 1e-6 * r[:, 5] * 1e-3).sum() * D2_over_g2_per_eV; w_ = np.unique(np.round(r[:, 5], 3))
    return Sg_, SD_, w_, r[:, 3].min(), r[:, 3].max(), r[:, 4].min(), r[:, 4].max()
bG = blocks[0]; bK = blocks[10]
a1_d = int(max(bK["rows"], key=lambda r: r[5])[2])                                     # A1' = highest-frequency mode at K in this run (159.8 meV, mode 6)
SgGd, SDGd, wGd, *eG = dfpt_S(bG, (5, 6)); SgKd, SDKd, wKd_, *eK = dfpt_S(bK, (a1_d,))
print(f"[route 1b DFPT prt, 16x16 k] q={bG['q']}: modes 5,6 hw {wGd} meV, E_k {eG[0]:.4f}..{eG[1]:.4f}, E_k+q {eG[2]:.4f}..{eG[3]:.4f}: S_g = {SgGd:.5f} eV^2 -> S_D(Gamma) = {SDGd:.2f} eV^2/A^2")
print(f"                              q={bK['q']}: A1' = mode {a1_d} hw {wKd_} meV, E_k+q {eK[2]:.4f}..{eK[3]:.4f}: S_g = {SgKd:.5f} eV^2 -> S_D(K) = {SDKd:.2f} eV^2/A^2 ; ratio {SDKd/SDGd:.3f}")
print(f"                EPW/DFPT: S_g Gamma {SgG/SgGd:.3f}, K {SgK/SgKd:.3f} ; S_D Gamma {SDG/SDGd:.3f}, K {SDK/SDKd:.3f} ; omega ratio Gamma {wG*1e3/wGd.mean():.3f}, K {wK*1e3/wKd_.mean():.3f}")

# ---------------- route 2: inversion of EPW gamma (1200^2, sigma 0.02, 10 K), per-mode HWHM
P = np.load(f"results/epw/phself_{a.phself_tag}.npz", allow_pickle=True)
from electron_defect_interaction.electron_phonon import phself
sp = phself.special_points(P["q"]); jG, jK = sp["G"][0], sp["K"][0]; i10 = int(np.argmin(np.abs(P["T"] - 10)))
Rp = dict(T=P["T"], omega=P["omega"], gamma_epw=P["gamma_epw"]); mE2 = list(phself.modes_E2g(Rp, jG)); mA2 = phself.mode_A1p(Rp, jK)   # by frequency / by character (largest gamma at K)
print(f"[route 2 modes] phonselfen: E2g(Gamma) = modes {mE2[0]+1}+{mE2[1]+1}, A1'(K) = mode {mA2+1} ({P['omega'][jK, mA2]:.2f} meV)")
gG = P["gamma_epw"][i10, jG, mE2].mean() * 1e-3; gK = P["gamma_epw"][i10, jK, mA2] * 1e-3; wG2 = P["omega"][jG, mE2].mean() * 1e-3; wK2 = P["omega"][jK, mA2] * 1e-3   # RAW EPW gamma___ (route 2 tests the raw output)
for vv, lab in ((v_epw, "v_F EPW"), (v_dft, "v_F bands.x")):
    SgG_inv = 16 * vv ** 2 * gG / (A_c * wG2); SgK_inv = 16 * vv ** 2 * gK / (A_c * wK2)
    print(f"[route 2, {lab} = {vv:.3f} eV.A] gamma_HWHM(E2g, per mode) = {gG*1e3:.4f} meV, hw {wG2*1e3:.2f} -> S_g(LO+TO) = {SgG_inv:.5f} eV^2, S_D(Gamma) = {SgG_inv*wG2*D2_over_g2_per_eV:.2f} eV^2/A^2 ; "
          f"gamma(A1') = {gK*1e3:.4f} meV, hw {wK2*1e3:.2f} -> S_g(K) = {SgK_inv:.5f} eV^2, S_D(K) = {SgK_inv*wK2*D2_over_g2_per_eV:.2f} eV^2/A^2")
print(f"[route 2 vs route 1a] S_g Gamma: inverted/direct = {16*v**2*gG/(A_c*wG2)/SgG:.3f} ; K: {16*v**2*gK/(A_c*wK2)/SgK:.3f}   (1 = EPW gamma consistent with EPW g under the Dirac-cone formula, v = v_F EPW)")
# forward: gamma predicted from route-1 sums
print(f"[forward] gamma_HWHM predicted from S_g(EPW): E2g per mode {A_c*wG*SgG/(16*v**2)*1e3:.4f} meV (EPW phonselfen {gG*1e3:.4f}), A1' {A_c*wK*SgK/(16*v**2)*1e3:.4f} meV (EPW {gK*1e3:.4f})")
# ---------------- comparison table (per-entry mean d^2 = S_D / n_entries: 8 at Gamma (2x2 block x LO+TO), 4 at K (2x2 block x A1'))
# Identification of Piscanec's <D^2>_F with the per-entry mean rests on two independent facts (to be confirmed in the PDFs):
#   (i) TB identity <D^2_Gamma>_F = (9/4)(dt/da)^2 = 45.6 eV^2/A^2 for dt/da = 4.5 eV/A (the value quoted with Piscanec's number), and (9/4) beta^2 is
#       exactly S_D(Gamma)/8 in the EPW displacement convention (per mode Tr(DD+) = 9 beta^2, two modes, 8 entries);
#   (ii) the Dirac-cone HWHM formula above with d^2 = 45.6 and v = 5.5 eV.A gives FWHM = 11 cm^-1, the literature E2g value.
LIT = {"Piscanec 2004 LDA": (45.6, 92.05), "Lazzeri 2008 GW": (62.8, 193.0)}     # eV^2/A^2, quoted from memory -> verify in Zotero
d2 = {"route 1a EPW (24x24 chain)": (SDG / 8, SDK / 4), "route 1b DFPT prt (16x16 k)": (SDGd / 8, SDKd / 4),
      "route 2 inversion of EPW gamma (Dirac cone)": (16 * v ** 2 * gG / (A_c * wG2) * wG2 * D2_over_g2_per_eV / 8, 16 * v ** 2 * gK / (A_c * wK2) * wK2 * D2_over_g2_per_eV / 4)}
print("\n[table] <D^2>_Gamma, <D^2>_K (eV^2/A^2), per-entry mean convention; ratio to Piscanec LDA in parentheses")
for name, (dg_, dk_) in list(d2.items()) + list(LIT.items()):
    print(f"   {name:46s} {dg_:8.1f} ({dg_/45.6:4.2f})   {dk_:8.1f} ({dk_/92.05:4.2f})   K/Gamma {dk_/dg_:4.2f}")
gG_lit = A_c * 45.6 / (2 * 5.5 ** 2 * D2_over_g2_per_eV); print(f"[check] Dirac-cone HWHM from Piscanec d^2 = 45.6, v = 5.5 eV.A: {gG_lit*1e3:.3f} meV = {gG_lit*1e3*8.06554:.2f} cm^-1 -> FWHM {2*gG_lit*1e3*8.06554:.1f} cm^-1 (literature ~10-11)")
np.savez(f"results/epw/d2_extract_{a.out_tag or a.tag}.npz", d2_table=np.array([(k, v1, v2) for k, (v1, v2) in list(d2.items()) + list(LIT.items())], dtype=object), units="S_g eV^2; S_D eV^2/A^2; v eV.A; omega eV; gamma eV (HWHM)", E_D=E_D, v_F_epw=v_epw, v_F_dft=v_dft, A_c=A_c, D2_over_g2_per_eV=D2_over_g2_per_eV,
         Sg_G_epw=SgG, Sg_K_epw=SgK, SD_G_epw=SDG, SD_K_epw=SDK, omega_G_epw=wG, omega_K_epw=wK, Sg_G_dfpt=SgGd, Sg_K_dfpt=SgKd, SD_G_dfpt=SDGd, SD_K_dfpt=SDKd,
         omega_G_dfpt=wGd.mean() * 1e-3, omega_K_dfpt=wKd_.mean() * 1e-3, gamma_G_hwhm_10K=gG, gamma_K_hwhm_10K=gK, Sg_G_inv=16 * v ** 2 * gG / (A_c * wG2), Sg_K_inv=16 * v ** 2 * gK / (A_c * wK2))
print(f"saved results/epw/d2_extract_{a.out_tag or a.tag}.npz")
