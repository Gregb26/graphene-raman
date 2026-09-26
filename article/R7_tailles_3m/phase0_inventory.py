"""R7 phase 0: vacancy-position rule (5..12), nbnd rule, window coverage from XML, timing fit."""
import os, re, sys, glob
import numpy as np
import xml.etree.ElementTree as ET
HA = 27.211386245988
GQ = os.environ["PROJECTS"] + "/graphene/qe"
SC = GQ + "/defects/super_cell"
SCR = os.path.expanduser("~/links/scratch/qe_tmp")

def read_scf_in(path):
    txt = open(path).read()
    nat = int(re.search(r"nat\s*=\s*(\d+)", txt).group(1))
    nbnd = int(re.search(r"nbnd\s*=\s*(\d+)", txt).group(1))
    cell = re.search(r"CELL_PARAMETERS\s+bohr\s*\n((?:.*\n){3})", txt).group(1)
    cell = np.array([[float(v) for v in l.split()] for l in cell.strip().splitlines()])
    pos = re.search(r"ATOMIC_POSITIONS\s+crystal\s*\n((?:\s*C\s+.*\n?)+)", txt).group(1)
    pos = np.array([[float(v) for v in l.split()[1:4]] for l in pos.strip().splitlines()])
    assert len(pos) == nat, (path, len(pos), nat)
    kp = re.search(r"K_POINTS\s+(\w+)", txt).group(1)
    return dict(nat=nat, nbnd=nbnd, cell=cell, pos=pos, kpoints=kp, txt=txt)

def wall(path):
    m = re.search(r"PWSCF\s+:\s+([\d.]+)s CPU\s+([\dhms .]+) WALL", open(path).read())
    s = m.group(2).strip()
    t = 0.0
    for val, unit in re.findall(r"([\d.]+)\s*([hms])", s):
        t += float(val) * dict(h=3600, m=60, s=1)[unit]
    return t

def niter(path):
    m = re.search(r"convergence has been achieved in\s+(\d+) iterations", open(path).read())
    return int(m.group(1)) if m else -1

def xml_eigs(save):
    root = ET.parse(save + "/data-file-schema.xml").getroot()
    bs = root.find(".//band_structure")
    ef = float(bs.find("fermi_energy").text) * HA
    nel = float(bs.find("nelec").text)
    nb = int(bs.find("nbnd").text)
    ks = bs.findall("ks_energies")
    assert len(ks) == 1
    e = np.array([float(v) for v in ks[0].find("eigenvalues").text.split()]) * HA
    occ = np.array([float(v) for v in ks[0].find("occupations").text.split()])
    npw = int(ks[0].find("npw").text)
    return dict(ef=ef, nel=nel, nbnd=nb, e=e, occ=occ, npw=npw)

print("=== A. Règle de position de la lacune (5..12) ===")
print(f"{'N':>3} {'nat_p':>5} {'nat_d':>5} {'idx1':>5} {'(i,j)':>8} {'SR':>2} {'frac':>22} {'centre-décalage':>18} {'nbnd_p':>6} {'occ_p':>5} {'+':>3} {'nbnd_d':>6} {'occ_d':>5} {'+':>3}")
rows = {}
for N in range(5, 13):
    t = f"{N}x{N}"
    p = read_scf_in(f"{SC}/{t}/pristine/scf.in"); d = read_scf_in(f"{SC}/{t}/defective/scf.in")
    assert np.allclose(p["cell"], d["cell"])
    # atom of p missing in d (minimal image)
    missing = []
    for ia, x in enumerate(p["pos"]):
        dx = d["pos"] - x; dx -= np.rint(dx)
        if np.min(np.linalg.norm(dx, axis=1)) > 1e-6:
            missing.append(ia)
    assert len(missing) == 1, (t, missing)
    ia = missing[0]; x = p["pos"][ia]
    # check that d = p minus that atom, same order
    keep = np.delete(p["pos"], ia, axis=0)
    same_order = np.allclose(keep, d["pos"], atol=1e-9)
    # sublattice: A if N*x - 1/3 integer, B if N*x - 2/3 integer
    fa = N * x[:2] - 1.0 / 3; fb = N * x[:2] - 2.0 / 3
    if np.allclose(fa, np.rint(fa), atol=1e-6): sr = "A"; ij = np.rint(fa).astype(int)
    elif np.allclose(fb, np.rint(fb), atol=1e-6): sr = "B"; ij = np.rint(fb).astype(int)
    else: sr = "?"; ij = (-1, -1)
    off = x[:2] - 0.5
    occ_p = int(p["nat"] * 4 // 2); occ_d = int(d["nat"] * 4 // 2)
    print(f"{N:>3} {p['nat']:>5} {d['nat']:>5} {ia+1:>5} {str(tuple(ij)):>8} {sr:>2} ({x[0]:.6f},{x[1]:.6f}) ({off[0]:+.4f},{off[1]:+.4f}) {p['nbnd']:>6} {occ_p:>5} {p['nbnd']-occ_p:>+3} {d['nbnd']:>6} {occ_d:>5} {d['nbnd']-occ_d:>+3}  ordre_conservé={same_order} kp={p['kpoints']}/{d['kpoints']}")
    rows[N] = dict(ia=ia, ij=tuple(ij), sr=sr, x=x, p=p, d=d)

# check generator hypothesis: positions = ((i+1/3)/N,(j+1/3)/N) A then ((i+2/3)/N,(j+2/3)/N) B, i outer, j inner
print("\n=== B. Hypothèse de construction des positions de la parfaite ===")
for N in range(5, 13):
    gen = []
    for i in range(N):
        for j in range(N):
            gen.append([(i + 1/3) / N, (j + 1/3) / N, 0.0]); gen.append([(i + 2/3) / N, (j + 2/3) / N, 0.0])
    gen = np.array(gen)
    ok = np.allclose(gen, rows[N]["p"]["pos"], atol=2e-9)
    # cell check: a = |a1| in bohr, a1 = (a*sqrt(3)/2 * N?, ...)
    cell = rows[N]["p"]["cell"]; a1 = np.linalg.norm(cell[0]) / N; ang = np.degrees(np.arccos(cell[0] @ cell[1] / (np.linalg.norm(cell[0]) * np.linalg.norm(cell[1]))))
    print(f"N={N:>2}: ordre (i,j,A,B) i extérieur j intérieur : {ok} ; a = {a1:.10f} bohr ({a1*0.529177210903:.6f} Å), angle {ang:.4f}°, c = {cell[2,2]:.4f} bohr ; vacance A-nearest-centre i=j={int(np.floor(N/2))} -> {'oui' if rows[N]['ij']==(int(np.floor(N/2)),)*2 and rows[N]['sr']=='A' else 'non'} ; vacance = (i,j)={rows[N]['ij']} {rows[N]['sr']}")

print("\n=== C. Proposition 15x15 / 18x18 (A, même règle) ===")
for N in (15, 18):
    cands = [(i, abs((i + 1/3) / N - 0.5)) for i in range(N)]
    i_near = min(cands, key=lambda c: c[1])[0]
    i_floor = N // 2
    print(f"N={N}: A le plus proche du centre : i=j={i_near} -> ({(i_near+1/3)/N:.10f}, idem) ; floor(N/2)={i_floor} ; indice 1-based dans la parfaite = {2*(i_near*N+i_near)+1} ; nat_p={2*N*N}, nat_d={2*N*N-1} ; occ_p={2*N*N*4//2}, occ_d={(2*N*N-1)*4//2} -> nbnd(+30) p={2*N*N*2+30}, d={(2*N*N-1)*2+30}")

print("\n=== D. XML : E_F, quadruplet, couverture de la fenêtre par nbnd (parfaite et lacune) ===")
print(f"{'cas':>16} {'nbnd':>5} {'npw':>8} {'E_F':>9} {'E_D(quad)':>9} {'quad spread':>11} {'e[nbnd-1]-E_D':>13} {'e[nbnd-1]-E_F':>13} {'n(e<=E_D+1)':>11} {'n(e<=E_D+1.5)':>13} {'marge(n bandes au-dessus de E_D+1)':>10}")
for N in range(5, 13):
    for s in "pd":
        t = f"{N}x{N}"; save = f"{SCR}/defect_{t}_{s}/defect_{t}_{s}.save"
        if not os.path.isdir(save): print(t, s, "absent"); continue
        z = xml_eigs(save); e = z["e"]; ef = z["ef"]
        # E_D from pristine quadruplet : 4 states nearest E_F (for 3m) -> use pristine of the same N
        zp = xml_eigs(f"{SCR}/defect_{t}_p/defect_{t}_p.save"); ep = np.sort(zp["e"]); efp = zp["ef"]
        nocc = int(round(zp["nel"] / 2))
        quad = ep[nocc - 2:nocc + 2]; E_D = quad.mean(); spread = quad.max() - quad.min()
        n1 = int(np.sum(e <= E_D + 1.0)); n15 = int(np.sum(e <= E_D + 1.5))
        print(f"{t+'_'+s:>16} {z['nbnd']:>5} {z['npw']:>8} {ef:>9.4f} {E_D:>9.4f} {spread:>11.4f} {e[-1]-E_D:>13.3f} {e[-1]-ef:>13.3f} {n1:>11} {n15:>13} {z['nbnd']-n1:>10}")

print("\n=== E. Temps de mur et ajustement t ∝ N^p ===")
Ns = np.arange(5, 13); td = []; tp = []; itd = []; itp = []
for N in Ns:
    t = f"{N}x{N}"
    td.append(wall(f"{SC}/{t}/defective/scf.out")); tp.append(wall(f"{SC}/{t}/pristine/scf.out"))
    itd.append(niter(f"{SC}/{t}/defective/scf.out")); itp.append(niter(f"{SC}/{t}/pristine/scf.out"))
td = np.array(td); tp = np.array(tp)
for N, a, b, c, d_ in zip(Ns, td, tp, itd, itp): print(f"N={N:>2}: d {a:7.1f} s ({c} it, {a/c:5.1f} s/it)   p {b:7.1f} s ({d_} it, {b/d_:5.1f} s/it)")
for lab, t, it in (("d", td, np.array(itd)), ("p", tp, np.array(itp)), ("d+p", td + tp, None)):
    for sub, mask in (("5-12", Ns >= 5), ("8-12", Ns >= 8)):
        pfit = np.polyfit(np.log(Ns[mask]), np.log(t[mask]), 1)
        pred = lambda n: np.exp(np.polyval(pfit, np.log(n)))
        print(f"  {lab:>3} N∈{sub}: p = {pfit[0]:.2f} ; t(15) = {pred(15)/60:.0f} min, t(18) = {pred(18)/60:.0f} min")
    if it is not None:
        for sub, mask in (("5-12", Ns >= 5), ("8-12", Ns >= 8)):
            pfit = np.polyfit(np.log(Ns[mask]), np.log(t[mask] / it[mask]), 1)
            pred = lambda n: np.exp(np.polyval(pfit, np.log(n)))
            print(f"  {lab:>3} par itération N∈{sub}: p = {pfit[0]:.2f} ; s/it(15) = {pred(15):.0f}, s/it(18) = {pred(18):.0f} -> ×{it[-1]} it (comme 12x12) = {pred(15)*it[-1]/60:.0f} / {pred(18)*it[-1]/60:.0f} min")

print("\n=== F. Mémoire QE estimée par rang (scf.out) et extrapolation N^4 ===")
for N in Ns:
    t = f"{N}x{N}"
    for s in ("defective", "pristine"):
        txt = open(f"{SC}/{t}/{s}/scf.out").read()
        m = re.search(r"Estimated max dynamical RAM per process >\s+([\d.]+) (\w+)", txt); m2 = re.search(r"Estimated total dynamical RAM >\s+([\d.]+) (\w+)", txt)
        fft = re.search(r"FFT dimensions: \(\s*(\d+),\s*(\d+),\s*(\d+)\)", txt); g = re.search(r"Dense  grid:\s+(\d+) G-vectors", txt)
        print(f"N={N:>2} {s[:1]}: par rang {m.group(1)} {m.group(2)}, total {m2.group(1)} {m2.group(2)}, FFT {fft.groups()}, G dense {g.group(1)}")
m12 = 946.65
for N in (15, 18): print(f"N={N}: par rang (×(N/12)^4 de 946.65 MB) = {m12*(N/12)**4/1024:.2f} GB ; total = {m12*64*(N/12)**4/1024:.0f} GB ; G dense ≈ {5487066*(N/12)**2/1e6:.2f} M ; FFT ({30*N},{30*N},192) ; npw(Γ, demi-sphère?) ≈ {343*(N/12)**2:.0f} k")
print("\n=== G. Taille des fichiers wfc/charge (12x12) et extrapolation ===")
for N in (9, 12):
    for s in "dp":
        t = f"{N}x{N}"; save = f"{SCR}/defect_{t}_{s}/defect_{t}_{s}.save"
        w = os.path.getsize(save + "/wfc1.hdf5"); c = os.path.getsize(save + "/charge-density.hdf5"); z = xml_eigs(save)
        print(f"{t}_{s}: wfc1.hdf5 {w/1e9:.3f} GB (npw {z['npw']} × nbnd {z['nbnd']} × 16 B = {z['npw']*z['nbnd']*16/1e9:.3f} GB), charge-density {c/1e6:.1f} MB")
w12 = 6636821124; c12 = 153641944; v12 = 428008254; pp12 = 327647000
for N in (15, 18):
    f = (N / 12) ** 2
    print(f"N={N}: wfc ≈ {w12*f*f/1e9:.1f} GB/cellule, charge-density ≈ {c12*f/1e6:.0f} MB, Vks (filplot) ≈ {v12*f/1e6:.0f} MB, pp.out (cube sur stdout) ≈ {pp12*f/1e6:.0f} MB ; 2 cellules : wfc {2*w12*f*f/1e9:.0f} GB")
