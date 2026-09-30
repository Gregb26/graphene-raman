#!/usr/bin/env python
"""
em1_check.py -- EM1 verifications a-d on the seedname_tb.dat produced by `restart = plot`
(graphene unit cell, 27x27x1 coarse grid, reference wannierization of chapter 4).

Reads (in the working directory):
    wannier_tb.dat            new _tb.dat (H block then r block)
    wannier_hr.dat            new _hr.dat
    ref/wannier_hr.dat.ref    reference _hr.dat (chapter 4)
    ref/wannier_tb.dat.ref    reference _tb.dat (chapter 4, same run)
    ref/wannier.wout.ref      reference .wout (final centres)
    wannier_wsvec.dat         new _wsvec.dat
Prints raw numbers only (no judgement). Usage: python em1_check.py [workdir]
"""
import os
import re
import sys

import numpy as np

wd = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.chdir(wd)


def read_tb(path):
    """Return A (3x3, rows = a_i, Ang), R (nrpts,3), ndegen, H (nrpts,nw,nw), r (nrpts,nw,nw,3)
    with H[ir,j,i] = <0j|H|Ri> and r[ir,j,i,:] = <0j|r|Ri> exactly as written (NOT divided by ndegen)."""
    with open(path) as f:
        header = f.readline()
        A = np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
        nw = int(f.readline()); nr = int(f.readline())
        nd = []
        while len(nd) < nr:
            nd += [int(x) for x in f.readline().split()]
        nd = np.array(nd[:nr])
        R = np.zeros((nr, 3), int); H = np.zeros((nr, nw, nw), complex)
        for ir in range(nr):
            assert f.readline().strip() == ""
            R[ir] = [int(x) for x in f.readline().split()]
            for _ in range(nw * nw):
                j, i, re_, im_ = f.readline().split()
                H[ir, int(j) - 1, int(i) - 1] = float(re_) + 1j * float(im_)
        Rr = np.zeros((nr, 3), int); r = np.zeros((nr, nw, nw, 3), complex)
        for ir in range(nr):
            assert f.readline().strip() == ""
            Rr[ir] = [int(x) for x in f.readline().split()]
            for _ in range(nw * nw):
                t = f.readline().split()
                j, i = int(t[0]) - 1, int(t[1]) - 1
                v = [float(x) for x in t[2:8]]
                r[ir, j, i, :] = [v[0] + 1j * v[1], v[2] + 1j * v[3], v[4] + 1j * v[5]]
        rest = f.read().strip()
    assert (R == Rr).all(), "R lists of H and r blocks differ"
    assert rest == "", "trailing data after the r block"
    return header, A, R, nd, H, r


def read_hr(path):
    with open(path) as f:
        header = f.readline()
        nw = int(f.readline()); nr = int(f.readline())
        nd = []
        while len(nd) < nr:
            nd += [int(x) for x in f.readline().split()]
        nd = np.array(nd[:nr])
        R = np.zeros((nr, 3), int); H = np.zeros((nr, nw, nw), complex)
        for ir in range(nr):
            for k in range(nw * nw):
                t = f.readline().split()
                if k == 0:
                    R[ir] = [int(x) for x in t[:3]]
                else:
                    assert [int(x) for x in t[:3]] == list(R[ir])
                H[ir, int(t[3]) - 1, int(t[4]) - 1] = float(t[5]) + 1j * float(t[6])
        assert f.read().strip() == ""
    return header, R, nd, H


def read_wout_centres(path):
    txt = open(path).read()
    fin = txt[txt.rfind("Final State"):]
    c = re.findall(r"WF centre and spread\s+\d+\s+\(\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\)\s+([-\d.]+)", fin)
    c = np.array(c, float)
    return c[:, :3], c[:, 3]


def read_wsvec(path):
    """dict (R tuple, iw, jw) -> list of shift T vectors ; header string."""
    d = {}
    with open(path) as f:
        header = f.readline().strip()
        while True:
            line = f.readline()
            if not line:
                break
            t = [int(x) for x in line.split()]
            R, iw, jw = tuple(t[:3]), t[3], t[4]
            n = int(f.readline())
            d[(R, iw, jw)] = [tuple(int(x) for x in f.readline().split()) for _ in range(n)]
    return header, d


hdr_tb, A, R, nd, H, r = read_tb("wannier_tb.dat")
hdr_tbref, Aref, Rref_tb, ndref_tb, Href_tb, rref = read_tb("ref/wannier_tb.dat.ref")
hdr_hr, R_hr, nd_hr, H_hr = read_hr("wannier_hr.dat")
hdr_hrref, Rref, ndref, Href = read_hr("ref/wannier_hr.dat.ref")
nr, nw = H.shape[0], H.shape[1]
print(f"new _tb.dat : {hdr_tb.strip()!r}  nw={nw} nrpts={nr}")
print(f"ref _hr.dat : {hdr_hrref.strip()!r}  nw={Href.shape[1]} nrpts={Href.shape[0]}")
print(f"lattice (Ang, rows a_i):\n{A}")

print("\n=== (a) H(R): new _hr.dat and H block of new _tb.dat vs reference _hr.dat ===")
print("same R list (order incl.)     new_hr vs ref_hr:", (R_hr == Rref).all(), "| tb vs ref_hr:", (R == Rref).all(),
      "| tb vs ref_tb:", (R == Rref_tb).all())
print("same ndegen                   new_hr vs ref_hr:", (nd_hr == ndref).all(), "| tb vs ref_hr:", (nd == ndref).all())
print("ndegen values: {1: %d, 2: %d, 3: %d}, sum = %d, N1*N2*N3 = 729" % ((nd == 1).sum(), (nd == 2).sum(), (nd == 3).sum(), (1.0 / nd).sum().round(6)))
print("max |H_new_hr - H_ref_hr|     = %.3e eV   (both files written with F12.6)" % np.abs(H_hr - Href).max())
print("max |H_tb - H_ref_hr|         = %.3e eV   (tb: E15.8 vs hr: F12.6 -> rounding of the hr file)" % np.abs(H - Href).max())
print("max |H_tb - H_ref_tb|         = %.3e eV   (new _tb.dat vs reference _tb.dat)" % np.abs(H - Href_tb).max())
print("max |r_tb - r_ref_tb|         = %.3e Ang  (new r block vs reference _tb.dat r block)" % np.abs(r - rref).max())
i0 = int(np.where((R == 0).all(axis=1))[0][0])
print("on-site energies <0i|H|0i> (eV):", np.round(H[i0].diagonal().real, 6))

print("\n=== (b) diagonal <0i|r|0i> vs final centres of the .wout ===")
cen, spr = read_wout_centres("ref/wannier.wout.ref")
rd = np.array([r[i0, i, i, :] for i in range(nw)])
print("max |Im <0i|r|0i>| = %.3e Ang" % np.abs(rd.imag).max())
print("  i   <0i|r|0i> (Ang)                         wout centre (Ang)                        diff (Ang)")
for i in range(nw):
    print("  %d   %s   %s   %s" % (i + 1, np.array2string(rd[i].real, precision=6, floatmode="fixed"),
                                   np.array2string(cen[i], precision=6, floatmode="fixed"),
                                   np.array2string(rd[i].real - cen[i], precision=1)))
dmax = np.abs(rd.real - cen).max(axis=0)
print("max |diff| per component (x,y,z) = %s Ang   (wout prints 6 decimals)" % np.array2string(dmax, precision=2))
# atoms and bond midpoints for reference (C1 at (1/3,1/3), C2 at (2/3,2/3))
tau1 = (A[0] + A[1]) / 3.0; tau2 = 2 * (A[0] + A[1]) / 3.0
mids = [(tau1 + tau2) / 2, (tau1 + tau2 - A[0]) / 2, (tau1 + tau2 - A[1]) / 2]
print("atoms: C1 =", np.round(tau1, 6), " C2 =", np.round(tau2, 6), " Ang")
print("bond midpoints C1-C2 (3 images):", [np.round(m, 6).tolist() for m in mids])

print("\n=== (c) hermiticity: max |r_ij(R) - conj(r_ji(-R))|, all (i,j,R) with (i,R) != (j,0) ===")
idx = {tuple(v): k for k, v in enumerate(R)}
missing = [tuple(v) for v in R if tuple(-v) not in idx]
print("R vectors whose -R is absent from the list:", len(missing))
herm = np.zeros(3); hermH = 0.0; cnt = 0
for ir, Rv in enumerate(R):
    jr = idx.get(tuple(-Rv))
    if jr is None:
        continue
    for i in range(nw):
        for j in range(nw):
            if i == j and (Rv == 0).all():
                continue
            cnt += 1
            herm = np.maximum(herm, np.abs(r[ir, i, j, :] - np.conj(r[jr, j, i, :])))
            hermH = max(hermH, abs(H[ir, i, j] - np.conj(H[jr, j, i])))
print("pairs checked:", cnt)
print("max |r_ij(R) - conj(r_ji(-R))| per component (x,y,z) = %s Ang" % np.array2string(herm, precision=3))
print("max |H_ij(R) - conj(H_ji(-R))|                       = %.3e eV" % hermH)
# also imaginary part of the diagonal at R=0 (should vanish for hermitian r)
print("max |Im r_ii(0)| = %.3e Ang ; max |Im H_ii(0)| = %.3e eV" % (np.abs(rd.imag).max(), np.abs(H[i0].diagonal().imag).max()))

print("\n=== (d) decay: max_{i!=j} |r_ij(R)| (any component) and max |H_ij(R)| vs |R| (Ang) ===")
Rc = R @ A                      # Cartesian R
dist = np.linalg.norm(Rc, axis=1)
off = ~np.eye(nw, dtype=bool)
rmax = np.array([np.abs(r[ir][off]).max() for ir in range(nr)])
rdiag = np.array([np.abs(np.array([r[ir, i, i, :] for i in range(nw)])).max() for ir in range(nr)])
hmax = np.array([np.abs(H[ir][off]).max() for ir in range(nr)])
hall = np.array([np.abs(H[ir]).max() for ir in range(nr)])
shells = np.unique(np.round(dist, 4))
print(" |R| (Ang)   n_R   max|r_ij| i!=j (Ang)   max|r_ii| (Ang)   max|H_ij| i!=j (eV)   max|H_ij| all (eV)")
for s in shells[:14]:
    m = np.isclose(np.round(dist, 4), s)
    print(" %8.4f   %3d   %12.3e         %12.3e       %12.3e          %12.3e" % (s, m.sum(), rmax[m].max(), rdiag[m].max(), hmax[m].max(), hall[m].max()))
print(" ... (%d shells in total, |R|max = %.4f Ang)" % (len(shells), dist.max()))
for s in shells[-3:]:
    m = np.isclose(np.round(dist, 4), s)
    print(" %8.4f   %3d   %12.3e         %12.3e       %12.3e          %12.3e" % (s, m.sum(), rmax[m].max(), rdiag[m].max(), hmax[m].max(), hall[m].max()))
# largest off-diagonal r at R=0 and at first neighbours
print("R=0 |r_ij| i!=j max = %.4e Ang ; |r_ij| at R=0 (Ang, abs, x/y/z max over comp):" % rmax[i0])
print(np.array2string(np.abs(r[i0]).max(axis=2), precision=4, suppress_small=True))

print("\n=== (e) wsvec ===")
hdr_ws, ws = read_wsvec("wannier_wsvec.dat")
print("header:", hdr_ws)
nent = len(ws); nmulti = sum(1 for v in ws.values() if len(v) > 1)
nshift = sum(1 for v in ws.values() if any(t != (0, 0, 0) for t in v))
print("entries (R,i,j):", nent, "= nrpts*nw*nw =", nr * nw * nw, "| entries with ndeg>1:", nmulti,
      "| entries with a non-zero shift T:", nshift)
# magnitude of r and H on the shifted entries
big = 0.0; bigH = 0.0
for (Rv, iw, jw), v in ws.items():
    if any(t != (0, 0, 0) for t in v):
        ir = idx[Rv]
        big = max(big, np.abs(r[ir, iw - 1, jw - 1, :]).max()); bigH = max(bigH, abs(H[ir, iw - 1, jw - 1]))
print("max |r_ij(R)| on entries with a non-zero T shift = %.3e Ang ; max |H_ij(R)| there = %.3e eV" % (big, bigH))
