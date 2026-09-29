#!/usr/bin/env python
"""
r8_driver.py -- R8 (reproduction of Kaasbjerg, PRB 101, 045433 (2020), Figs. 13 and 14) pilot.

Sub-commands
    extract   step 7a: vector extraction of Fig. 13 (top panel, vacancies V_A) and reading of Fig. 14 (left column, V_A)
              from the article PDF with `mutool trace`; axis calibration from the tick marks, residuals, curve metrics.
              Runs on the login node (seconds). No computation of our DOS.
    prep      M2 dense 9x9 and 12x12 (v2 gate, md5, gauge gate), Wannier rotation, recentering, V_loc for R_cut 2/3/4, box, E_D,
              neighbours of the vacancy -> cache/prep_<S>.npz
    g0        every g0 of the campaign (local_green_batch, 5 WF, one path; pi block = index selection) -> cache/g0*_*.npy
    gate      step 1: gates P1 (Born), P2 (Gamma vs resonance_9x9.npz), P3 (linear term = Lloyd), P4 (regression on R9 B / R.3 / C.3);
              information norms. Exit code 3 if a gate fails.
    dos       steps 2 and 2 bis (DOS 9x9, pristine, eta_unique; LDOS of the three neighbours)
    spec      step 3 (A_k on Gamma-K-M)
    sens      step 5 (DOS c_i = 1 %: N_k^int, R_cut, 12x12, eta_G)
    fig       step 7b (superposition on Kaasbjerg Fig. 13), 2 bis table, figures, tables

Outputs: out/ (json, npz, csv, md), fig/ ; nothing in results/ (TEST, decision 7 of 2026-09-28). Energies in eV relative to E_D
(Wannier Dirac point). DOS in states / eV / unit cell / spin (pi block). Concentrations c_i = defects per unit cell.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.environ["PROJECTS"], "graphene-raman")
PDF = os.path.join(PROJ, "article", "R8_kaasbjerg", "ref", "kaasbjerg_2020_prb101_045433.pdf")
OUT = os.path.join(HERE, "out")
FIG = os.path.join(HERE, "fig")

# ----------------------------------------------------------------------------------------------------------------------
# mutool trace reader
# ----------------------------------------------------------------------------------------------------------------------


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def trace_page(pdf, page):
    """Run `mutool trace` on one page; return the list of paths with device coordinates (pt, origin top-left, y down).

    Each path: dict(kind, color (tuple of floats), lw, dash, subs=[(n, 2) arrays], ncurve, nsub, npts). Curves (curveto) are
    replaced by their end point: enough for bounding boxes, centres and the polylines of plotted data (lineto only)."""
    xml = subprocess.run(["mutool", "trace", pdf, str(page)], check=True, capture_output=True).stdout
    root = ET.fromstring(xml)
    paths = []
    for el in root.iter():
        if el.tag not in ("stroke_path", "fill_path"):
            continue
        a, b, c, d, e, f = map(float, el.get("transform").split())
        M, t = np.array([[a, b], [c, d]]), np.array([e, f])
        subs, cur, ncurve = [], [], 0
        for c_ in el:
            if c_.tag == "moveto":
                if cur:
                    subs.append(cur)
                cur = [(float(c_.get("x")), float(c_.get("y")))]
            elif c_.tag == "lineto":
                cur.append((float(c_.get("x")), float(c_.get("y"))))
            elif c_.tag == "curveto":
                cur.append((float(c_.get("x3")), float(c_.get("y3"))))
                ncurve += 1
        if cur:
            subs.append(cur)
        if not subs:
            continue
        subs = [np.array(s) @ M + t for s in subs]
        paths.append(dict(kind=el.tag, color=tuple(float(v) for v in el.get("color", "0").split()), lw=el.get("linewidth"),
                          dash=el.get("dash"), subs=subs, ncurve=ncurve, nsub=len(subs), npts=sum(len(s) for s in subs)))
    return paths


def bbox(p):
    P = np.vstack(p["subs"])
    return P[:, 0].min(), P[:, 0].max(), P[:, 1].min(), P[:, 1].max()


def inside(p, box, tol=0.02):
    x0, x1, y0, y1 = box
    a, b, c, d = bbox(p)
    return a >= x0 - tol and b <= x1 + tol and c >= y0 - tol and d <= y1 + tol


def is_color(p, rgb, tol=2e-3):
    c = p["color"]
    c = (c[0],) * 3 if len(c) == 1 else c
    return all(abs(u - v) < tol for u, v in zip(c, rgb))


INK = (0.13725, 0.12157, 0.12549)            # axes / text colour of the article's figures


def ticks(paths, frame, side, length=(0.9, 1.5)):
    """Tick marks: short ink strokes (0.216 pt) perpendicular to the given frame edge, starting on it. Returns sorted positions."""
    x0, x1, y0, y1 = frame
    pos = []
    for p in paths:
        if p["kind"] != "stroke_path" or not is_color(p, INK) or p["npts"] != 2:
            continue
        (xa, ya), (xb, yb) = p["subs"][0][0], p["subs"][0][-1]
        # one end of the tick on the frame edge, the tick pointing inwards or outwards
        if side == "bottom" and abs(xa - xb) < 1e-3 and length[0] <= abs(ya - yb) <= length[1] \
                and min(abs(ya - y1), abs(yb - y1)) < 0.05 and x0 - 0.05 <= xa <= x1 + 0.05:
            pos.append(xa)
        if side == "left" and abs(ya - yb) < 1e-3 and length[0] <= abs(xa - xb) <= length[1] \
                and min(abs(xa - x0), abs(xb - x0)) < 0.05 and y0 - 0.05 <= ya <= y1 + 0.05:
            pos.append(ya)
    return np.array(sorted(set(np.round(pos, 3))))


# glyph classes of the tick/legend digits (fonts drawn as paths): signature (nsub, npts, ncurve), minus = flat 4-point glyph > 4 pt wide
def glyph_class(p):
    a, b, c, d = bbox(p)
    w = b - a
    sig = (p["nsub"], p["npts"], p["ncurve"])
    if p["ncurve"] == 0 and p["npts"] == 4 and w > 4:
        return "-"
    if (p["ncurve"] == 0 and p["npts"] == 4) or sig == (1, 10, 6):
        return "."
    if sig in ((2, 18, 16), (2, 27, 25)):
        return "0"
    if sig in ((1, 12, 0), (1, 14, 4)):
        return "1"
    if sig == (1, 22, 16):
        return "2"
    if sig == (1, 23, 16):
        return "5"
    return "?"


def label_near(paths, x, y, side):
    """Glyph string of the tick label next to a tick (below for the x axis, left for the y axis), ordered left to right."""
    out = []
    for p in paths:
        if p["kind"] != "fill_path" or not is_color(p, INK) or p["npts"] < 4:
            continue
        a, b, c, d = bbox(p)
        if (b - a) > 7 or (d - c) > 8:
            continue
        cx, cy = 0.5 * (a + b), 0.5 * (c + d)
        if side == "bottom" and abs(cx - x) < 12 and 0.5 < cy - y < 10:
            out.append((a, glyph_class(p)))
        if side == "left" and 0.5 < x - cx < 20 and abs(cy - y) < 3.5:     # labels <= 18.7 pt from the frame, axis titles >= 22 pt
            out.append((a, glyph_class(p)))
    return [g for _, g in sorted(out)]


def matches(glyphs, label):
    """True if the glyph classes spell the label string."""
    return "".join(glyphs) == label.replace("−", "-")


def calibrate(pos, values):
    """Least-squares affine map value = a*pos + b from the tick positions; returns (a, b, max |residual| in axis units)."""
    pos, values = np.asarray(pos, float), np.asarray(values, float)
    A = np.vstack([pos, np.ones_like(pos)]).T
    (a, b), *_ = np.linalg.lstsq(A, values, rcond=None)
    return float(a), float(b), float(np.max(np.abs(A @ [a, b] - values)))


# ----------------------------------------------------------------------------------------------------------------------
# curve metrics
# ----------------------------------------------------------------------------------------------------------------------


def half_width(x, y, i):
    """Full width at half maximum of the peak of y at index i (linear interpolation of the two half-height crossings); NaN if a
    crossing is not found inside the array."""
    h = 0.5 * y[i]
    left = np.where(y[:i] < h)[0]
    right = np.where(y[i:] < h)[0]
    if len(left) == 0 or len(right) == 0:
        return float("nan"), float("nan"), float("nan")
    j = left[-1]
    xl = x[j] + (h - y[j]) * (x[j + 1] - x[j]) / (y[j + 1] - y[j])
    k = i + right[0]
    xr = x[k - 1] + (h - y[k - 1]) * (x[k] - x[k - 1]) / (y[k] - y[k - 1])
    return float(xr - xl), float(xl), float(xr)


def local_maxima(x, y):
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))[0] + 1
    return [(float(x[j]), float(y[j])) for j in i]


# ----------------------------------------------------------------------------------------------------------------------
# step 7a
# ----------------------------------------------------------------------------------------------------------------------

# Labels read on the rendered figure (pdftoppm), checked below against the glyph shapes of the vector file.
FIG13_X_LABELS = ["−1.0", "−0.5", "0.0", "0.5", "1.0"]           # Energy (eV), left to right
FIG13_Y_LABELS = ["0.00", "0.05", "0.10", "0.15", "0.20"]        # DOS (eV^-1), bottom to top
FIG13_LEGEND = ["pristine", "0.1", "1.0"]                          # legend lines, top to bottom (c_i in %)
FIG14_Y_LABELS = ["−1.0", "−0.5", "0.0", "0.5", "1.0"]           # Energy (eV), bottom to top
BLUE, GREEN, WHITE, RED = (0.0, 0.0, 1.0), (0.0, 0.50196, 0.0), (1.0, 1.0, 1.0), (1.0, 0.0, 0.0)


def fig13(paths, log):
    """Fig. 13, top panel (V_A). Frame = the white axes background whose top edge is above y = 150 pt (the upper panel)."""
    frames = [p for p in paths if p["kind"] == "fill_path" and is_color(p, WHITE) and p["npts"] == 4 and bbox(p)[0] > 300]
    frames = sorted(frames, key=lambda p: bbox(p)[2])
    frame = bbox(frames[0])                                         # upper panel
    log(f"[fig13] frame (pt) x = [{frame[0]:.3f}, {frame[1]:.3f}], y = [{frame[2]:.3f}, {frame[3]:.3f}]")
    xt = ticks(paths, frame, "bottom")
    yt = ticks(paths, frame, "left")[::-1]                          # bottom to top (device y decreases upwards)
    assert len(xt) == len(FIG13_X_LABELS) and len(yt) == len(FIG13_Y_LABELS), (xt, yt)
    chk = {}
    for x, lab in zip(xt, FIG13_X_LABELS):
        g = label_near(paths, x, frame[3], "bottom"); chk[f"x {lab}"] = ("".join(g), matches(g, lab))
    for y, lab in zip(yt, FIG13_Y_LABELS):
        g = label_near(paths, frame[0], y, "left"); chk[f"y {lab}"] = ("".join(g), matches(g, lab))
    ax, bx, rx = calibrate(xt, [float(s.replace("−", "-")) for s in FIG13_X_LABELS])
    ay, by, ry = calibrate(yt, [float(s) for s in FIG13_Y_LABELS])
    log(f"[fig13] x ticks {xt.tolist()} -> {ax:.6f} eV/pt, max residual {rx:.2e} eV ; y ticks {yt.tolist()} -> {ay:.6e} /eV/pt, "
        f"max residual {ry:.2e} /eV")
    log(f"[fig13] tick labels vs glyph shapes: {chk}")
    # curves: long polylines inside the frame; legend samples: 2-point strokes of the same style inside the frame
    box = frame
    style = lambda p: (p["color"], p["dash"])
    curves = [p for p in paths if p["kind"] == "stroke_path" and p["npts"] > 50 and inside(p, box)]
    legend = sorted([p for p in paths if p["kind"] == "stroke_path" and p["npts"] == 2 and inside(p, box) and float(p["lw"]) > 0.5],
                    key=lambda p: bbox(p)[2])
    assert len(curves) == 3 and len(legend) == 3, (len(curves), len(legend))
    lab_of = {style(l): FIG13_LEGEND[i] for i, l in enumerate(legend)}
    leg_chk = {}
    for l in legend[1:]:
        a, b, c, d = bbox(l)
        g = [glyph_class(p) for p in sorted([p for p in paths if p["kind"] == "fill_path" and is_color(p, INK) and inside(p, box)
                                             and abs(0.5 * (bbox(p)[2] + bbox(p)[3]) - c) < 4.5 and 445 < bbox(p)[0] < 462],
                                            key=lambda p: bbox(p)[0])]
        leg_chk[lab_of[style(l)]] = "".join(g)
    log(f"[fig13] legend (top to bottom): {[ (lab_of[style(l)], l['color'], l['dash']) for l in legend]} ; digit glyphs after '=' : {leg_chk}")
    data = {}
    for p in curves:
        P = np.vstack(p["subs"])
        e, rho = ax * P[:, 0] + bx, ay * P[:, 1] + by
        assert np.all(np.diff(e) > 0), "curve not monotonic in energy"
        data[lab_of[style(p)]] = (e, rho)
        log(f"[fig13] curve {lab_of[style(p)]:>8s}: {len(e)} points, e in [{e.min():+.4f}, {e.max():+.4f}] eV, color {p['color']}, dash {p['dash']}")
    cal = dict(frame_pt=frame, x_ticks_pt=xt.tolist(), x_labels=FIG13_X_LABELS, x_scale_eV_per_pt=ax, x_offset=bx, x_resid_eV=rx,
               y_ticks_pt=yt.tolist(), y_labels=FIG13_Y_LABELS, y_scale_per_pt=ay, y_offset=by, y_resid=ry, label_check=chk,
               legend_check=leg_chk, legend=FIG13_LEGEND)
    return data, cal


def fig13_metrics(data, log):
    e = np.round(np.arange(-1.2, 1.2 + 5e-4, 1e-3), 6)
    r = {k: np.interp(e, *v, left=np.nan, right=np.nan) for k, v in data.items()}
    m = {}
    p = r["pristine"]
    ip = int(np.nanargmin(p))
    m["pristine"] = dict(min_pos_eV=float(e[ip]), min_val=float(p[ip]), rho_at_0=float(np.interp(0.0, e, p)),
                         rho_at_m1=float(np.interp(-1.0, e, p)), rho_at_p1=float(np.interp(1.0, e, p)))
    for side, (lo, hi) in (("neg", (-0.5, -0.1)), ("pos", (0.1, 0.5))):
        s = (e >= lo) & (e <= hi)
        alpha, beta = np.polyfit(np.abs(e[s]), p[s], 1)
        m["pristine"][f"slope_{side}_{lo:+.1f}_{hi:+.1f}"] = float(alpha)
        m["pristine"][f"intercept_{side}"] = float(beta)
    for c in ("0.1", "1.0"):
        y = r[c]; d = y - p
        s = (e >= -1.0) & (e <= 0.0)
        es, ds = e[s], d[s]
        i = int(np.argmax(ds))
        fw, xl, xr = half_width(es, ds, i)
        ic = int(np.nanargmin(y))
        m[c] = dict(drho_max_pos_eV=float(es[i]), drho_max=float(ds[i]), drho_fwhm_eV=fw, drho_half_left=xl, drho_half_right=xr,
                    drho_max_at_window_edge=bool(i in (0, len(es) - 1)),
                    rho_local_maxima_m1_0=local_maxima(es, y[s]), rho_at_0=float(np.interp(0.0, e, y)), rho_min_pos_eV=float(e[ic]),
                    rho_min=float(y[ic]), rms_drho_m1_p1=float(np.sqrt(np.nanmean(d[(e >= -1) & (e <= 1)] ** 2))))
    for k, v in m.items():
        log(f"[fig13 metrics] {k}: {v}")
    return e, r, m


def fig14(paths, log):
    """Fig. 14, left column (V_A): upper panel c_i = 0.1 %, lower panel c_i = 1 %. Frames = white axes backgrounds with x < 170 pt."""
    frames = sorted([p for p in paths if p["kind"] == "fill_path" and is_color(p, WHITE) and p["npts"] == 4 and p["ncurve"] == 0
                     and bbox(p)[1] < 170 and (bbox(p)[3] - bbox(p)[2]) > 80], key=lambda p: bbox(p)[2])
    res = {}
    for name, fr in zip(("0.1", "1.0"), frames[:2]):
        frame = bbox(fr)
        yt = ticks(paths, frame, "left")[::-1]
        xt = ticks(paths, frame, "bottom")
        assert len(yt) == len(FIG14_Y_LABELS), yt
        chk = {}
        for y, lab in zip(yt, FIG14_Y_LABELS):
            g = label_near(paths, frame[0], y, "left"); chk[lab] = ("".join(g), matches(g, lab))
        ay, by, ry = calibrate(yt, [float(s.replace("−", "-")) for s in FIG14_Y_LABELS])
        E = lambda y: ay * y + by
        # colour-map cells (4-point grey polygons, no curves) : column width and row height
        cells = [bbox(p) for p in paths if p["kind"] == "fill_path" and p["ncurve"] == 0 and p["npts"] == 4 and inside(p, frame)
                 and (bbox(p)[3] - bbox(p)[2]) < 0.5]
        ch = np.median([c[3] - c[2] for c in cells])
        cw = np.median([c[1] - c[0] for c in cells if (c[1] - c[0]) < 2.5 and (c[1] - c[0]) > 1.5])
        xK = float(xt[0]) if len(xt) == 1 else float("nan")
        ncol = (frame[1] - frame[0]) / cw
        # white dots (circles: white fill paths with curves)
        centre = lambda p: (0.5 * (bbox(p)[0] + bbox(p)[1]), 0.5 * (bbox(p)[2] + bbox(p)[3]))
        dots = np.array([centre(p) for p in paths if p["kind"] == "fill_path" and is_color(p, WHITE) and p["ncurve"] > 0
                         and frame[0] <= centre(p)[0] <= frame[1] and frame[2] <= centre(p)[1] <= frame[3]])   # centre in the frame
        # K column grey profile
        colK = sorted([(bbox(p), p["color"]) for p in paths if p["kind"] == "fill_path" and p["ncurve"] == 0 and p["npts"] == 4
                       and inside(p, frame) and (bbox(p)[3] - bbox(p)[2]) < 0.5 and bbox(p)[0] <= xK <= bbox(p)[1]], key=lambda t: t[0][2])
        prof = np.array([[E(0.5 * (b[2] + b[3])), np.mean(c)] for b, c in colK])
        sat = prof[prof[:, 1] <= prof[:, 1].min() + 1e-6, 0]
        red = [p for p in paths if p["kind"] == "stroke_path" and is_color(p, RED) and bbox(p)[0] >= frame[0] - 5 and bbox(p)[1] <= frame[1] + 5
               and bbox(p)[2] >= frame[2] - 10 and bbox(p)[3] <= frame[3] + 10]
        r = dict(frame_pt=frame, y_ticks_pt=yt.tolist(), y_labels=FIG14_Y_LABELS, y_scale_eV_per_pt=ay, y_offset=by, y_resid_eV=ry,
                 label_check=chk, x_ticks_pt=xt.tolist(), xK_pt=xK, cell_width_pt=float(cw), cell_height_pt=float(ch),
                 cell_height_eV=float(abs(ay) * ch), n_columns=float(ncol), n_dots=int(len(dots)),
                 K_column_saturated_eV=[float(sat.min()), float(sat.max())] if len(sat) else None)
        log(f"[fig14 {name} %] frame {tuple(round(v, 3) for v in frame)} ; y ticks {yt.tolist()} -> {ay:.6f} eV/pt, max residual {ry:.2e} eV ; "
            f"labels {chk} ; x ticks {xt.tolist()} ; cells {cw:.4f} x {ch:.4f} pt ({abs(ay) * ch * 1e3:.2f} meV), {ncol:.2f} columns ; "
            f"{len(dots)} white dots ; K column darkest (saturated) cells over [{sat.min():+.4f}, {sat.max():+.4f}] eV")
        if len(dots):
            col = np.rint((dots[:, 0] - xK) / cw).astype(int)             # column index relative to K (<0: towards Gamma, >0: towards M)
            eps = E(dots[:, 1])
            atK = np.sort(eps[col == 0])
            r["dots"] = [dict(col=int(c_), x_pt=float(x), y_pt=float(y), eps_eV=float(e_)) for c_, (x, y), e_ in
                         sorted(zip(col, dots, eps), key=lambda t: (t[0], t[2]))]
            r["K_dots_eV"] = atK.tolist()
            r["K_gap_eV"] = float(atK[-1] - atK[0]) if len(atK) == 2 else None
            # two branches: in every column holding two dots, lower = min, upper = max; a jump = |change| > 0.15 eV between adjacent
            # columns (regular spacing of the dots away from the resonance: <= 0.11 eV per column)
            per = {}
            for c_, e_ in zip(col, eps):
                per.setdefault(int(c_), []).append(float(e_))
            two = {c_: sorted(v) for c_, v in per.items() if len(v) == 2}
            jumps = {}
            for br, j in (("lower", 0), ("upper", 1)):
                cs = sorted(two)
                jumps[br] = [dict(cols=[a_, b_], eps_eV=[two[a_][j], two[b_][j]]) for a_, b_ in zip(cs[:-1], cs[1:])
                             if b_ == a_ + 1 and abs(two[b_][j] - two[a_][j]) > 0.15]
            r["branch_jumps"] = jumps
            r["dots_per_column"] = {str(k): v for k, v in sorted(per.items())}
            steps = {side: np.abs(np.diff([two[c_][0] for c_ in sorted(two) if (c_ < 0 if side == "Gamma" else c_ > 0)])) for side in ("Gamma", "M")}
            r["lower_branch_step_eV"] = {k: [float(v.min()), float(np.median(v)), float(v.max())] for k, v in steps.items()}
            log(f"[fig14 {name} %] K column dots: {np.round(atK, 4).tolist()} eV -> gap {r['K_gap_eV']} eV ; branch jumps > 0.15 eV between adjacent "
                f"columns: {jumps} ; lower-branch step between adjacent columns (min, median, max): {r['lower_branch_step_eV']}")
        rb = []
        for p in red:
            for s in p["subs"]:
                # energy of the dashed branches at the frame top: horizontal distance from K on each side
                rb.append(dict(first_pt=s[0].tolist(), last_pt=s[-1].tolist()))
        r["red_dashed_subpaths"] = rb
        r["K_profile"] = prof.tolist()
        res[name] = r
    return res


def cmd_extract(a):
    os.makedirs(os.path.join(OUT, "7a"), exist_ok=True); os.makedirs(FIG, exist_ok=True)
    logf = open(os.path.join(OUT, "7a", "extract_log.txt"), "w")

    def log(s):
        print(s, flush=True); logf.write(s + "\n")
    ver = subprocess.run(["mutool", "-v"], capture_output=True, text=True)
    log(f"[7a] PDF {PDF} md5 {md5(PDF)} ; {(ver.stdout + ver.stderr).strip()} ; python {sys.version.split()[0]}")
    p14 = trace_page(PDF, 14)
    p15 = trace_page(PDF, 15)
    log(f"[7a] page 14: {len(p14)} paths ; page 15: {len(p15)} paths")
    data, cal13 = fig13(p14, log)
    e, r, m13 = fig13_metrics(data, log)
    for k, (x, y) in data.items():
        np.savetxt(os.path.join(OUT, "7a", f"fig13_VA_{'pristine' if k == 'pristine' else 'c' + k.replace('.', 'p') + 'pc'}.csv"),
                   np.c_[x, y], delimiter=",", header="energy_eV,DOS_per_eV (axis units of the article: 'DOS (eV^-1)')", comments="")
    np.savetxt(os.path.join(OUT, "7a", "fig13_VA_interp_1meV.csv"), np.c_[e, r["pristine"], r["0.1"], r["1.0"]], delimiter=",",
               header="energy_eV,pristine,c0p1pc,c1p0pc (linear interpolation of the extracted polylines)", comments="")
    f14 = fig14(p15, log)
    if "dots" in f14.get("1.0", {}):
        np.savetxt(os.path.join(OUT, "7a", "fig14_VA_c1p0pc_dots.csv"),
                   np.array([[d["col"], d["x_pt"], d["y_pt"], d["eps_eV"]] for d in f14["1.0"]["dots"]]), delimiter=",",
                   header="column_rel_K,x_pt,y_pt,energy_eV", comments="")
    json.dump(dict(pdf=PDF, pdf_md5=md5(PDF), fig13_calibration=cal13, fig13_metrics=m13,
                   fig14={k: {kk: vv for kk, vv in v.items() if kk != "K_profile"} for k, v in f14.items()}),
              open(os.path.join(OUT, "7a", "7a_results.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(OUT, "7a", "fig14_K_profiles.npz"), **{f"c{k.replace('.', 'p')}": np.array(v["K_profile"]) for k, v in f14.items()})
    control_figure(data, cal13, f14, log)
    logf.close()


def control_figure(data, cal13, f14, log):
    """Extracted points drawn over the rendered page (pdftoppm, 400 dpi): visual check of the extraction for Greg."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.join(PROJ, "figures", "memoire.mplstyle"))
    sys.path.insert(0, os.path.join(PROJ, "scripts"))
    from _palette import CYCLE
    dpi = 400
    tmp = os.path.join(OUT, "7a", "render")
    for pg in (14, 15):
        subprocess.run(["pdftoppm", "-r", str(dpi), "-f", str(pg), "-l", str(pg), "-png", "-singlefile", PDF, f"{tmp}_p{pg}"], check=True)
    s = dpi / 72.0
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.3))
    im = plt.imread(f"{tmp}_p14.png")
    x0, x1, y0, y1 = cal13["frame_pt"]
    a1.imshow(im[int((y0 - 12) * s):int((y1 + 24) * s), int((x0 - 38) * s):int((x1 + 6) * s)],
              extent=((x0 - 38), (x1 + 6), (y1 + 24), (y0 - 12)))
    for (k, (x, y)), col in zip(sorted(data.items()), CYCLE):
        xp = (x - cal13["x_offset"]) / cal13["x_scale_eV_per_pt"]; yp = (y - cal13["y_offset"]) / cal13["y_scale_per_pt"]
        a1.plot(xp, yp, ".", ms=1.2, color=col, label=("parfait" if k == "pristine" else rf"$c_i$ = {k.replace('.', ',')}\,\%"))
    a1.set_axis_off(); a1.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=6, markerscale=4); a1.set_title("(a) Fig. 13 (haut) : points extraits", fontsize=8)
    im = plt.imread(f"{tmp}_p15.png")
    r = f14["1.0"]; x0, x1, y0, y1 = r["frame_pt"]
    a2.imshow(im[int((y0 - 8) * s):int((y1 + 16) * s), int((x0 - 30) * s):int((x1 + 6) * s)], extent=((x0 - 30), (x1 + 6), (y1 + 16), (y0 - 8)))
    if "dots" in r:
        d = np.array([[q["x_pt"], q["y_pt"]] for q in r["dots"]])
        a2.plot(d[:, 0], d[:, 1], "o", ms=2.2, mfc="none", mec=CYCLE[1], mew=0.5, label="points blancs extraits")
    a2.axvline(r["xK_pt"], color=CYCLE[2], lw=0.4)
    a2.set_axis_off(); a2.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), fontsize=6); a2.set_title(r"(b) Fig. 14 ($V_A$, $c_i$ = 1\,\%)", fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"7a_controle.{ext}"), dpi=300)
    for pg in (14, 15):
        os.remove(f"{tmp}_p{pg}.png")
    log(f"[7a] control figure: {os.path.join(FIG, '7a_controle.pdf')}")


# ======================================================================================================================
# Computations (steps 1, 2, 2 bis, 3, 5, 7b) -- GO of 2026-09-28
# ======================================================================================================================
import types

BOHR = 0.529177210903
SIZES = ("9x9", "12x12")
NW = 5                                   # production wannierisation: 0-2 sp2 on C1, 3 p_z on C1, 4 p_z on C2
PI = (3, 4)
WIN = 1.2                                # eV, half window of the DOS and A_k around E_D
C_LIST = (0.001, 0.01)                   # c_i, defects per unit cell
ETA_DOS, ETA_PRIS, ETA_AK = 0.050, 0.015, 0.025
NK_ALL, NK_MAIN = (300, 600, 900), 900
GRIDS_DOS, GRIDS_PRIS, GRID_SENS, GRID_RES = (300, 600), (300, 1200), 600, 240
NPATH = 601
VARIANTS = ("tel_quel", "aligne")
R9W = os.path.join(os.path.dirname(HERE), "R9_controles")
A1_JSON = os.path.join(R9W, "a", "a1_results.json")
CONV_C = ("c_i = defauts par maille unitaire (Kaasbjerg, PRB 101, 045433, eq. 22-24) ; c_at = c_i/2 ; n_i = c_i/A_cell "
          "(A_cell = 5,266 A^2 : c_i = 1 % <=> 1,90e13 cm^-2)")
DOS_UNITS = "etats / eV / maille / spin (bloc pi : 2 bandes)"

# smoke mode (R8_SMOKE=1): tiny grids, separate cache and output directories, comparisons with R9 / production skipped
SMOKE = os.environ.get("R8_SMOKE") == "1"
_SMALL = {240: 24, 300: 24, 600: 30, 900: 36, 1200: 48}
CACHE = os.path.join(HERE, "cache_smoke" if SMOKE else "cache")
OUTC = os.path.join(HERE, "out_smoke" if SMOKE else "out")


def NK(n):
    return _SMALL[n] if SMOKE else n


def ensure(*p):
    d = os.path.join(*p)
    os.makedirs(d, exist_ok=True)
    return d


def rlog(s):
    print(s, flush=True)
    with open(os.path.join(ensure(OUTC), "r8_log.txt"), "a") as f:
        f.write(f"[{time.strftime('%F %T')}] {s}\n")


def _jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(type(o))


def save_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=_jsonable)


def git_head():
    return subprocess.run(["git", "-C", PROJ, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()


def sha12(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()[:12]


def workers():
    return max(1, int(os.environ.get("R8_WORKERS") or os.environ.get("SLURM_CPUS_PER_TASK") or 1))


_MODS = {}


def mods():
    """Package imports, done on first use (step 7a runs without h5py / the scipy-stack modules)."""
    if not _MODS:
        from electron_defect_interaction.config import load_production, dense_paths
        from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
        from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
        from electron_defect_interaction.wannier import wannier_interpolation as wi
        from electron_defect_interaction.defects.many_body import local_tmatrix as lt
        from electron_defect_interaction.defects.many_body import disorder_average as da
        from electron_defect_interaction.defects.many_body.pole_criterion import local_t_cache
        from electron_defect_interaction.utils.lattice import build_k_path
        _MODS.update(load_production=load_production, dense_paths=dense_paths, qe_io=qe_io, matrix_io=matrix_io,
                     wannier_provenance=wannier_provenance, read_w90_mat=read_w90_mat, read_w90_HR=read_w90_HR, wi=wi, lt=lt, da=da,
                     local_t_cache=local_t_cache, build_k_path=build_k_path)
    return types.SimpleNamespace(**_MODS)


def cfg():
    return mods().load_production(verbose=False)


def provenance(S=None, **kw):
    m = mods(); c = cfg()
    p = dict(campaign="R8", status="TEST", smoke=SMOKE, date=time.strftime("%F %T"), head=git_head(), convention_c=CONV_C,
             dos_units=DOS_UNITS, sha256_disorder_average=sha12(m.da.__file__), sha256_local_tmatrix=sha12(m.lt.__file__),
             sha256_driver=sha12(os.path.abspath(__file__)), a1_json=A1_JSON, a1_md5=md5(A1_JSON), R_cut_prod=c["R_cut"],
             eta_t_prod=c["eta_eV"], energies="eV relatives a E_D (point de Dirac de Wannier, minimum du gap pi/pi* sur 90^2)")
    if S:
        P = load_prep(S)
        p.update(size=S, M2_md5=str(P["M2_md5"]), wannier_run_id=str(P["run_id"]), sublattice=str(P["sublattice"]),
                 vacancy_wf=int(P["vac_wf"]), E_D=float(P["E_D"]), C_N_eV=C_N(S))
    p.update(kw)
    return p


# ------------------------------------------------------------------------------------------------ prep: M2 -> V_loc
def load_wannier(dp, k):
    """U, U_dis reordered on the k of the dense .save, H(R); gauge gate (refuses if the manifest does not pass)."""
    m = mods()
    paths = m.wannier_provenance.load_wannier_checked(dp["manifest"])
    U, kU = m.read_w90_mat(paths["u"]); U = U[m.wi._match_kpoint_order(kU, k)]
    Ud, kUd = m.read_w90_mat(paths["u_dis"]); Ud = Ud[m.wi._match_kpoint_order(kUd, k)]
    Hwr, Rw, nd = m.read_w90_HR(paths["tb"])
    return U, Ud, Hwr, Rw, nd, json.load(open(dp["manifest"]))["run_id"]


def hwk_chunked(Hwr, Rw, nd, k, chunk=100_000, want_u=False):
    """H_W(k) (and eigenpairs) by chunks of k: same formula as Hwr_to_Hwk (the (nk, nrpts) phase array is the memory limit)."""
    m = mods(); nw = Hwr.shape[1]
    H = np.empty((len(k), nw, nw), complex)
    E = np.empty((len(k), nw)) if want_u else None
    U = np.empty((len(k), nw, nw), complex) if want_u else None
    for s in range(0, len(k), chunk):
        h, e, u = m.lt.Hwr_to_Hwk(Hwr, Rw, k[s:s + chunk], ndegen=nd)
        H[s:s + chunk] = h
        if want_u:
            E[s:s + chunk] = e; U[s:s + chunk] = u
    return (H, E, U) if want_u else H


def md5_listed(fname):
    for line in open(os.path.join(PROJ, "results", "M2", "MD5SUMS_2026-09-25.txt")):
        h, n = line.split()
        if n == fname:
            return h
    return None


def prep_size(S):
    m = mods(); c = cfg(); dp = m.dense_paths(c, S); n = int(S.split("x")[0]); t0 = time.time()
    side = json.load(open(dp["mfile"].replace(".npy", ".json")))
    h_file, h_list = md5(dp["mfile"]), md5_listed(os.path.basename(dp["mfile"]))
    if not (h_file == h_list == side.get("md5")):
        raise SystemExit(f"[prep] {S} : md5 de M2 {h_file} != MD5SUMS {h_list} / sidecar {side.get('md5')} -- refus")
    M = m.matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=m.matrix_io.UNIT_CELL, units=m.matrix_io.EV,
                                   require_normalization=m.matrix_io.M_NORM_V2)
    k = m.qe_io.get_k_red(dp["uc"]); MP = m.wi._infer_mp_grid(k)
    U, Ud, Hwr, Rw, nd, run_id = load_wannier(dp, k)
    Mwk = m.wi.Mbk_to_Mwk(M, U, Ud); del M
    Mwr, R = m.wi.Mwk_to_Mwr(Mwk, k, MP)
    Rn, R_d = m.lt.recenter_mwr(Mwr, R, MP); m.lt.mwr_locality(Mwr, Rn)
    A = m.qe_io.get_A_volume(dp["uc"])[0] * BOHR                       # columns a_i (Angstrom)
    Bc = 2 * np.pi * np.linalg.inv(A).T                                  # columns b_i (1/Angstrom)
    x_red = np.asarray(m.qe_io.get_x_red(dp["uc"]), float)
    K = np.array(c["K_red"], float)
    iKc = int(np.argmin(np.linalg.norm(np.mod(k - K + 0.5, 1.0) - 0.5, axis=1)))
    half_Mwk_K = float(0.5 * np.trace(Mwk[3:5, iKc, 3:5, iKc]).real)   # V^dagger M2 V at (K, K), pi block
    kg = np.vstack([np.zeros(3), k[iKc], k[37], k[250]])                # grid points only (labels, audit P-c4)
    Mf = m.wi.Mwr_to_Mwk_pairs(Mwr, Rn, kg, kg)
    Mfull_pi = np.array([Mf[3:5, i, 3:5, i] for i in range(len(kg))]); del Mwk, Mf
    out = {}
    D = np.asarray(MP[:2], int)
    for rc in (2, 3, 4):
        sel = np.linalg.norm(Rn, axis=1) <= rc + 1e-9
        Rloc, iloc = Rn[sel], np.where(sel)[0]
        V, herm = m.lt.extract_V_loc(Mwr, Rn, Rloc)
        cbox = np.mod(np.asarray(R, int)[iloc][:, :2], D)                # box of the N x N supercell (as R9 box_geometry)
        inbox = np.all(cbox < n, axis=1)
        edge = np.where(inbox, np.min(np.minimum(cbox, n - 1 - cbox), axis=1), -1)
        out.update({f"Rloc_rc{rc}": Rloc, f"V_rc{rc}": V, f"herm_rc{rc}": herm, f"inbox_rc{rc}": inbox, f"edge_rc{rc}": edge,
                    f"Mloc_rc{rc}": np.ascontiguousarray(Mwr[:, iloc][:, :, :, iloc])})
    del Mwr
    _, E_ref, _ = m.lt.Hwr_to_Hwk(Hwr, Rw, m.lt.mp_grid(90, 90, 1), ndegen=nd)
    gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
    # vacancy orbital (largest on-site p_z element at R = 0) against the sidecar, and its three nearest p_z neighbours
    V3, R3 = out["V_rc3"], out["Rloc_rc3"]; L0 = int(np.where((R3 == 0).all(axis=1))[0][0])
    vac_wf = PI[int(np.argmax([abs(V3[L0 * NW + w, L0 * NW + w]) for w in PI]))]
    if vac_wf != {"A": 3, "B": 4}[side["vacancy_sublattice"]]:
        raise SystemExit(f"[prep] {S} : orbitale de la lacune {vac_wf} incompatible avec le sidecar ({side['vacancy_sublattice']})")
    other = 7 - vac_wf
    tau_v, tau_o = x_red[vac_wf - 3], x_red[other - 3]                    # C1 -> WF 3, C2 -> WF 4 (projections C1:sp2;pz C2:pz)
    cands = [(i, j, 0) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    dist = np.array([np.linalg.norm(A @ (np.array(r, float) + tau_o - tau_v)) for r in cands])
    order = np.argsort(dist)
    nn_cells = np.array([cands[i] for i in order[:3]]); nn_d = dist[order[:3]]
    assert np.ptp(nn_d) < 1e-6 and dist[order[3]] > nn_d.max() + 0.1, dist
    if S == "9x9":
        assert {tuple(r) for r in nn_cells} == {(0, 0, 0), (-1, 0, 0), (0, -1, 0)}, nn_cells   # NN_CELLS of R4 / R9
    idx_nn = np.array([int(np.where((R3 == r).all(axis=1))[0][0]) * NW + other for r in nn_cells])
    np.savez(os.path.join(ensure(CACHE), f"prep_{S}.npz"), **out, Hwr=Hwr, Rw=Rw, nd=nd, E_D=E_D, R_d=R_d, MP=np.array(MP), A=A, Bc=Bc,
             x_red=x_red, kg=kg, Mfull_pi=Mfull_pi, half_Mwk_K=half_Mwk_K, vac_wf=vac_wf, nn_cells=nn_cells, nn_dist=nn_d,
             idx_nn=idx_nn, M2_md5=h_file, run_id=run_id, sublattice=side["vacancy_sublattice"], n=n)
    summ = dict(S=S, M2_md5=h_file, run_id=run_id, sublattice=side["vacancy_sublattice"], vac_wf=vac_wf, R_d=R_d, MP=MP, E_D=E_D,
                nn_cells=nn_cells, nn_dist_A=nn_d, half_Mwk_K_eV=half_Mwk_K,
                clusters={rc: dict(nL=len(out[f"Rloc_rc{rc}"]), herm=out[f"herm_rc{rc}"], out_of_box=int((~out[f"inbox_rc{rc}"]).sum()),
                                   edge0=int((out[f"edge_rc{rc}"] == 0).sum())) for rc in (2, 3, 4)}, elapsed_s=time.time() - t0)
    save_json(os.path.join(ensure(OUTC, "prep"), f"prep_{S}.json"), summ)
    rlog(f"[prep] {S} : M2 md5 {h_file} (= MD5SUMS = sidecar), Wannier {run_id}, R_d {R_d.tolist()}, E_D {E_D:.6f} eV, lacune WF {vac_wf} "
         f"({side['vacancy_sublattice']}), voisins {nn_cells.tolist()} ({nn_d[0]:.4f} A), amas {[summ['clusters'][r]['nL'] for r in (2, 3, 4)]}, "
         f"hors boite {[summ['clusters'][r]['out_of_box'] for r in (2, 3, 4)]}, derniere maille de la boite {[summ['clusters'][r]['edge0'] for r in (2, 3, 4)]} ; "
         f"{time.time() - t0:.0f} s")


def cmd_prep(a):
    for S in SIZES:
        prep_size(S)


_PREP = {}


def load_prep(S):
    if S not in _PREP:
        Z = np.load(os.path.join(CACHE, f"prep_{S}.npz"))
        _PREP[S] = {k: Z[k] for k in Z.files}
    return _PREP[S]


def C_N(S):
    """C_N = mean of the plateau (i) of R9 A.1 (key C_i_eV of the file, never the cache)."""
    return float(json.load(open(A1_JSON))["sizes"][S]["C_i_eV"])


def V_of(P, rc, var, S):
    """tel_quel: V_loc of M2 ; aligne: V_loc - C_N P_box (approximation (i) of R9: identity on the 5 orbitals of the box cells)."""
    V = P[f"V_rc{rc}"]
    if var == "tel_quel":
        return V
    return V - C_N(S) * np.diag(np.repeat(P[f"inbox_rc{rc}"].astype(float), NW))


def Mloc_of(P, rc, var, S):
    """M_W restricted to the cluster, (NW, nL, NW, nL), same variant as V_of (built from M_W directly, for gate P1)."""
    M = P[f"Mloc_rc{rc}"].copy()
    if var == "aligne":
        w = np.arange(NW)
        for L in np.where(P[f"inbox_rc{rc}"])[0]:
            M[w, L, w, L] -= C_N(S)
    return M


def pi_index(nL):
    return np.array([L * NW + w for L in range(nL) for w in PI])


def hk_pi(P, k):
    return np.ascontiguousarray(hwk_chunked(P["Hwr"], P["Rw"], P["nd"], k)[:, 3:5, 3:5])


# ------------------------------------------------------------------------------------------------ g0 (one path: local_green_batch)
def egrid_res(E_D):
    """Energy grid of resonance_metrics.py (« res »): step eta/ne_per_eta = 2.5 meV on E_D -+ (e_window + eta)."""
    c = cfg(); eta = c["eta_eV"]; de = eta / c["ne_per_eta"]; ew = c["e_window_eV"]
    return np.arange(E_D - ew - eta, E_D + ew + eta + de, de), de


def grid_index(egrid, E_D, grid):
    if grid == "res":
        return np.arange(len(egrid))
    w = np.where(np.abs(egrid - E_D) <= WIN + 1e-9)[0]
    return w if grid == "win" else w[::4]


def g0_file(S, rc, nk, eta, grid, deriv=False):
    return os.path.join(CACHE, f"g0{'p' if deriv else ''}_{S}_rc{rc}_nk{NK(nk)}_eta{int(round(eta * 1e3))}_{grid}.npy")


def g0_meta(S, rc, nk, eta, grid, deriv):
    P = load_prep(S); E_D = float(P["E_D"]); egrid, de = egrid_res(E_D); idx = grid_index(egrid, E_D, grid)
    return dict(S=S, rc=rc, nk=NK(nk), eta=eta, grid=grid, deriv=deriv, E_D=E_D, e0=float(egrid[idx][0]), nE=len(idx), de=de,
                nL=len(P[f"Rloc_rc{rc}"]))


def g0_make(S, rc, nk, eta, grid, deriv=False):
    m = mods(); f = g0_file(S, rc, nk, eta, grid, deriv); meta = g0_meta(S, rc, nk, eta, grid, deriv)
    if os.path.exists(f) and os.path.exists(f.replace(".npy", ".json")):
        old = json.load(open(f.replace(".npy", ".json")))
        if all(old.get(k) == v for k, v in meta.items()):
            rlog(f"[g0] {os.path.basename(f)} : deja calcule"); return
    P = load_prep(S); E_D = float(P["E_D"]); egrid, _ = egrid_res(E_D); idx = grid_index(egrid, E_D, grid); t0 = time.time()
    k_int = m.lt.mp_grid(NK(nk), NK(nk), 1)
    Hk = hwk_chunked(P["Hwr"], P["Rw"], P["nd"], k_int)
    g = m.lt.local_green_batch(Hk, k_int, P[f"Rloc_rc{rc}"], egrid[idx], eta, deriv=deriv); del Hk
    np.save(f, g)
    save_json(f.replace(".npy", ".json"), dict(meta, elapsed_s=time.time() - t0, head=git_head()))
    rlog(f"[g0] {os.path.basename(f)} : {len(idx)} energies, {NK(nk)}^2, {(time.time() - t0) / 60:.1f} min")


def g0_get(S, rc, nk, eta, grid, deriv=False):
    f = g0_file(S, rc, nk, eta, grid, deriv)
    old = json.load(open(f.replace(".npy", ".json")))
    meta = g0_meta(S, rc, nk, eta, grid, deriv)
    if any(old.get(k) != v for k, v in meta.items()):
        raise SystemExit(f"[g0] {os.path.basename(f)} : metadonnees {old} != attendu {meta}")
    return np.load(f)


def g0_for(S, rc, nk, eta, grid):
    """g0 (5 WF) of the cluster R_cut on the requested grid, taken from the stored file (sub-grid and, for R_cut 2, sub-block of R_cut 3)."""
    P = load_prep(S); E_D = float(P["E_D"]); egrid, _ = egrid_res(E_D); want = grid_index(egrid, E_D, grid)
    if rc == 2:
        g, eg = g0_for(S, 3, nk, eta, grid)
        R3, R2 = P["Rloc_rc3"], P["Rloc_rc2"]
        pos = [int(np.where((R3 == r).all(axis=1))[0][0]) for r in R2]
        idx = np.array([L * NW + w for L in pos for w in range(NW)])
        return g[:, idx][:, :, idx], eg
    for stored in ("res", "win"):
        if os.path.exists(g0_file(S, rc, nk, eta, stored)):
            g = g0_get(S, rc, nk, eta, stored); have = grid_index(egrid, E_D, stored)
            pos = np.searchsorted(have, want)
            assert np.array_equal(have[pos], want), "grille demandee absente du g0 stocke"
            return g[pos], egrid[want]
    raise SystemExit(f"[g0] aucun g0 pour {S} R_cut {rc} {nk} eta {eta}")


def g0_tasks():
    eta_p = cfg()["eta_eV"]
    tasks = [(S, 3, nk, eta_p, "res", False) for S in SIZES for nk in NK_ALL]
    tasks += [("9x9", 4, NK_MAIN, eta_p, "win", False), ("9x9", 3, 600, ETA_DOS, "win", False), ("9x9", 3, 600, ETA_AK, "win", False),
              ("9x9", 3, 300, eta_p, "win4", True)]
    return tasks


def cmd_g0(a):
    for t in g0_tasks():
        g0_make(*t)


def t_pi(S, rc, nk, eta_t, var, grid="win"):
    """pi-block t = V_pi [1 - g0_pi V_pi]^-1 (local_t_cache) and its reduction tau(D) (Q1)."""
    m = mods(); P = load_prep(S)
    g, eg = g0_for(S, rc, nk, eta_t, grid)
    nL = len(P[f"Rloc_rc{rc}"]); ip = pi_index(nL); V = V_of(P, rc, var, S)
    tp = m.local_t_cache(V[np.ix_(ip, ip)], np.ascontiguousarray(g[:, ip][:, :, ip])); del g
    Du, tau = m.da.tbar_reduce(tp, P[f"Rloc_rc{rc}"], 2)
    return dict(P=P, eg=eg, E_D=float(P["E_D"]), tau=tau, Du=Du, t=tp)


def dos_par(Hk, tau, Du, k, cs, eg, eta, linear=False):
    """dos_average (Q4) over energy blocks in Python threads (1 BLAS thread each; LAPACK and BLAS release the GIL)."""
    m = mods(); W = workers(); nE = len(eg)
    parts = [p for p in np.array_split(np.arange(nE), min(nE, 4 * W)) if len(p)]

    def run(ix):
        return ix, m.da.dos_average(Hk, tau[ix], Du, k, cs, eg[ix], eta, linear=linear)
    with ThreadPoolExecutor(W) as ex:
        res = list(ex.map(run, parts))
    scalar = np.ndim(cs) == 0
    rho = np.zeros(nE) if scalar else np.zeros((len(np.atleast_1d(cs)), nE))
    rho0 = np.zeros(nE); lin = np.zeros(nE) if linear else None
    for ix, r in res:
        rho[..., ix] = r["rho"]; rho0[ix] = r["rho0"]
        if linear:
            lin[ix] = r["drho_lin"]
    out = dict(rho=rho, rho0=rho0)
    if linear:
        out["drho_lin"] = lin
    return out


def dmetrics(eg, rho, rho0):
    """Maximum of rho - rho0 on [-1, 0] eV (position, height, FWHM), rho(E_D), local maxima of rho on [-1, 0], minimum of rho."""
    d = rho - rho0
    msk = (eg >= -1.0 - 1e-9) & (eg <= 0.0 + 1e-9)
    w = np.where(msk)[0]; i_rel = int(np.argmax(d[msk])); i = int(w[i_rel])
    fw, xl, xr = half_width(eg, d, i)
    ic = int(np.argmin(rho))
    return dict(drho_max_pos_eV=float(eg[i]), drho_max=float(d[i]), drho_fwhm_eV=fw, drho_half_left=xl, drho_half_right=xr,
                drho_max_at_window_edge=bool(i_rel in (0, len(w) - 1)), rho_at_ED=float(np.interp(0.0, eg, rho)),
                rho_local_maxima_m1_0=local_maxima(eg[msk], rho[msk]), rho_min_pos_eV=float(eg[ic]), rho_min=float(rho[ic]))


def lor(x, eta):
    return (eta / np.pi) / (x * x + eta * eta)


# ------------------------------------------------------------------------------------------------ gate (step 1)
def res_setup(S):
    """States of the production output grid (240^2) and the « res » energy grid, as resonance_metrics.py / R9."""
    P = load_prep(S); c = cfg(); E_D = float(P["E_D"]); egrid, de = egrid_res(E_D)
    k_out = mods().lt.mp_grid(NK(GRID_RES), NK(GRID_RES), 1)
    _, E_out, U_out = hwk_chunked(P["Hwr"], P["Rw"], P["nd"], k_out, want_u=True)
    sel = np.abs(E_out - E_D) <= c["e_window_eV"]
    js = np.rint((E_out - egrid[0]) / de).astype(int)
    return dict(P=P, E_D=E_D, egrid=egrid, de=de, eta=c["eta_eV"], k_out=k_out, E_out=E_out, U_out=U_out, sel=sel, js=js,
                K=np.array(c["K_red"], float))


def res_observables(st, V, g0):
    """Curve Gamma_T (Lorentzian average of the on-shell rates of the 240^2 states, step 2.5 meV) and -Im Tbar(K) = -Im Tr/2 of the
    pi block, from t = V[1 - g0 V]^-1 (5 WF); on-shell elements by Q1 + Q2 in the band basis (-2 Im Tbar^{nn}_k(eps_j))."""
    m = mods(); P = st["P"]; egrid, E_D = st["egrid"], st["E_D"]
    t = m.local_t_cache(V, g0)
    Du, tau = m.da.tbar_reduce(t, P["Rloc_rc3"], NW); del t
    E_out, sel, js = st["E_out"], st["sel"], st["js"]
    G = np.full(E_out.shape, np.nan)
    for j in np.unique(js[sel]):
        ik, nb = np.nonzero(sel & (js == j))
        uk, inv = np.unique(ik, return_inverse=True)
        Tn = m.da.tbar_k(tau[j], Du, st["k_out"][uk], U=st["U_out"][uk])
        G[ik, nb] = -2.0 * Tn[inv, nb, nb].imag
    Es, Gs = E_out[sel], G[sel]
    curve = np.empty(len(egrid))
    for s in range(0, len(egrid), 64):
        w = lor(egrid[s:s + 64, None] - Es[None, :], st["eta"]); curve[s:s + 64] = (w * Gs[None, :]).sum(1) / w.sum(1)
    TK = m.da.tbar_k(tau, Du, st["K"][None])[:, 0]
    tr = 0.5 * np.trace(TK[:, 3:5, 3:5], axis1=1, axis2=2)
    return dict(peak_GT_fine=float(egrid[int(np.argmax(curve))] - E_D), max_GT_fine=float(curve.max()),
                peak_ImTbar=float(egrid[int(np.argmax(-tr.imag))] - E_D), max_mImTbar=float((-tr.imag).max()),
                ReTbar_at_ED=float(np.interp(E_D, egrid, tr.real)), curve=curve, tr=tr, G=G)


def ldos_nn(P, V, g0, egrid, E_D):
    """LDOS of the three p_z neighbours of the vacancy (cluster_ldos, G = g0 + g0 T g0 on the cluster), sum of the three, and its
    maximum on [-1, 0] eV."""
    rho, rho0 = mods().lt.cluster_ldos(g0, V, P["idx_nn"])
    tot = rho.sum(1); eg = egrid - E_D; msk = (eg >= -1.0 - 1e-9) & (eg <= 0.0 + 1e-9); i = int(np.argmax(tot[msk]))
    return dict(at_eV=float(eg[msk][i]), value=float(tot[msk][i]), per_orbital_at_max=rho[msk][i].tolist()), rho, rho0


def rot120(Bc, k):
    c_, s_ = np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)
    Rz = np.array([[c_, -s_, 0], [s_, c_, 0], [0, 0, 1.0]])
    return np.linalg.solve(Bc, Rz @ (Bc @ np.asarray(k, float)))


def cmd_gate(a):
    m = mods(); c = cfg(); eta_p = c["eta_eV"]; d = ensure(OUTC, "gate"); t0 = time.time()
    res = dict(prov=provenance()); ok = {}
    rng = np.random.default_rng(1)
    K = np.array(c["K_red"], float); Kp = np.array([1 / 3, 2 / 3, 0.0])
    # ---- P1 : Born, Q1 + Q2 with t <- V_loc against Mwr_to_Mwk_pairs on M_W restricted to the cluster
    kP1 = np.vstack([np.zeros(3), K, [0.5, 0, 0], K + [0.01, 0, 0], np.c_[rng.random((4, 2)), np.zeros(4)]])
    p1 = {}
    for S in SIZES:
        P = load_prep(S); R3 = P["Rloc_rc3"]
        for var in VARIANTS:
            Du, tau = m.da.tbar_reduce(V_of(P, 3, var, S), R3, NW, wfs=PI)
            TB = m.da.tbar_k(tau, Du, kP1)
            Mr = m.wi.Mwr_to_Mwk_pairs(Mloc_of(P, 3, var, S), R3, kP1, kP1)
            ref = np.array([Mr[3:5, i, 3:5, i] for i in range(len(kP1))])
            err = float(np.abs(TB - ref).max() / np.abs(ref).max())
            p1[f"{S} {var}"] = dict(rel_err=err, ok=bool(err <= 1e-10), half_tr_K_eV=float(0.5 * np.trace(TB[1]).real))
        Du, tau = m.da.tbar_reduce(P["V_rc3"], R3, NW, wfs=PI); TBg = m.da.tbar_k(tau, Du, P["kg"]); Mf = P["Mfull_pi"]
        p1[f"{S} information"] = dict(k_grid=P["kg"], rel_diff_Rcut3_vs_full=[float(np.abs(TBg[i] - Mf[i]).max() / np.abs(Mf[i]).max())
                                                                             for i in range(len(Mf))],
                                      half_tr_K_VdagM2V=float(P["half_Mwk_K"]), half_tr_K_Mw_full=float(0.5 * np.trace(Mf[1]).real),
                                      half_tr_K_Rcut3=float(0.5 * np.trace(TBg[1]).real))
    ok["P1"] = all(v["ok"] for v in p1.values() if "ok" in v); res["P1"] = p1
    rlog(f"[gate] P1 Born : {({k: v['rel_err'] for k, v in p1.items() if 'rel_err' in v})} -> {'OK' if ok['P1'] else 'ECHEC'}")
    # ---- P2 : -2 Im Tbar^nn_k(eps_j) (5 WF, band basis, Q1 + Q2) against G_T of results/M2/resonance_9x9.npz (N_k^int 300)
    S = "9x9"; P = load_prep(S); E_D = float(P["E_D"]); egrid, de = egrid_res(E_D)
    if not SMOKE:
        Z = np.load(os.path.join(PROJ, "results", "M2", "resonance_9x9.npz"), allow_pickle=True)
        eg_diff = float(np.abs(Z["egrid"] - egrid).max()) if Z["egrid"].shape == egrid.shape else float("inf")
        E_out, GT = Z["E_out"], Z["G_T"]; k_out = m.lt.mp_grid(GRID_RES, GRID_RES, 1)
        gapmin = np.min(np.abs(E_out[:, :, None] - E_out[:, None, :]) + 1e9 * np.eye(NW)[None], axis=2)
        cand = np.argwhere(np.isfinite(GT) & (np.abs(E_out - E_D) <= WIN) & (gapmin > 1e-3))
        pick = cand[rng.choice(len(cand), 18, replace=False)]
        near = cand[(E_out[cand[:, 0], cand[:, 1]] - E_D > -0.3) & (E_out[cand[:, 0], cand[:, 1]] - E_D < -0.1)]
        pick = np.vstack([pick, near[np.argmax(GT[near[:, 0], near[:, 1]])]])
        js = np.rint((E_out[pick[:, 0], pick[:, 1]] - egrid[0]) / de).astype(int); uj = np.unique(js)
        t = m.local_t_cache(V_of(P, 3, "tel_quel", S), g0_get(S, 3, 300, eta_p, "res")[uj])
        Du5, tau5 = m.da.tbar_reduce(t, P["Rloc_rc3"], NW); del t
        ours = []
        for (ik, nb), j in zip(pick, js):
            _, _, U1 = m.lt.Hwr_to_Hwk(P["Hwr"], P["Rw"], k_out[ik][None], ndegen=P["nd"])
            Tn = m.da.tbar_k(tau5[int(np.searchsorted(uj, j))], Du5, k_out[ik][None], U=U1)
            ours.append(-2.0 * Tn[0, nb, nb].imag)
        ours = np.array(ours); ref = GT[pick[:, 0], pick[:, 1]]
        err = float(np.abs(ours - ref).max() / np.abs(ref).max())
        ok["P2"] = bool(err <= 1e-10 and eg_diff < 1e-12)
        res["P2"] = dict(rel_err=err, egrid_max_diff=eg_diff, n_states=len(pick), states=pick, eps_rel=E_out[pick[:, 0], pick[:, 1]] - E_D,
                         Gamma_ours_eV=ours, Gamma_prod_eV=ref, ok=ok["P2"])
        rlog(f"[gate] P2 Gamma : {len(pick)} etats, ecart relatif max {err:.2e}, grille « res » identique a {eg_diff:.1e} -> {'OK' if ok['P2'] else 'ECHEC'}")
    # ---- P3 : linear term of Q4 = Lloyd (1/pi) Im Tr[t dg0/deps], pi block, 300^2 = internal grid, eta_G = eta_t
    w4 = grid_index(egrid, E_D, "win4"); R3 = P["Rloc_rc3"]; ip = pi_index(len(R3))
    g = g0_get(S, 3, 300, eta_p, "res")[w4]
    gp = g0_get(S, 3, 300, eta_p, "win4", deriv=True)
    Vp = V_of(P, 3, "tel_quel", S)[np.ix_(ip, ip)]
    tp = m.local_t_cache(Vp, np.ascontiguousarray(g[:, ip][:, :, ip])); del g
    Du, tau = m.da.tbar_reduce(tp, R3, 2)
    k_int = m.lt.mp_grid(NK(300), NK(300), 1)
    lin = dos_par(hk_pi(P, k_int), tau, Du, k_int, 0.0, egrid[w4], eta_p, linear=True)["drho_lin"]
    gpp = gp[:, ip][:, :, ip]
    lloyd = np.array([np.trace(tp[j] @ gpp[j]).imag / np.pi for j in range(len(w4))])
    err = float(np.abs(lin - lloyd).max() / np.abs(lloyd).max())
    ok["P3"] = bool(err <= 1e-10); res["P3"] = dict(rel_err=err, n_energies=len(w4), ok=ok["P3"])
    rlog(f"[gate] P3 terme lineaire = Lloyd : {len(w4)} energies, ecart relatif max {err:.2e} -> {'OK' if ok['P3'] else 'ECHEC'}")
    # ---- P4 : regression on R9 (B, R.3: peaks of the Gamma_T curve and of -Im Tbar(K); C.3: LDOS of the three neighbours)
    B9 = {} if SMOKE else {"tel_quel": json.load(open(os.path.join(R9W, "b", "b_results.json"))),
                           "aligne": json.load(open(os.path.join(R9W, "b", "b_results_plateau.json")))}
    C9 = {} if SMOKE else json.load(open(os.path.join(R9W, "c", "c_results.json")))["C3"]
    key9 = {"tel_quel": "brut", "aligne": "aligne_plateau"}
    p4, p4ok = {}, True
    for S in SIZES:
        st = res_setup(S); P = st["P"]
        for nk in NK_ALL:
            g0 = g0_get(S, 3, nk, eta_p, "res")
            for var in VARIANTS:
                V = V_of(P, 3, var, S)
                o = res_observables(st, V, g0)
                row = {k: v for k, v in o.items() if k not in ("curve", "tr", "G")}
                if S == "9x9" and nk in (600, 900) and var == "tel_quel":
                    row["ldos_nn"] = ldos_nn(P, V, g0, st["egrid"], st["E_D"])[0]
                if not SMOKE:
                    r9 = B9[var][S][key9[var]][str(nk)]
                    cmp_ = dict(peak_GT=abs(row["peak_GT_fine"] - r9["peak_GT_fine"]) < 1e-6,
                                peak_ImTbar=abs(row["peak_ImTbar"] - r9["peak_ImTbar"]) < 1e-6,
                                max_GT=abs(row["max_GT_fine"] / r9["max_GT_fine"] - 1) < 1e-8,
                                max_ImTbar=abs(row["max_mImTbar"] / r9["max_mImTbar"] - 1) < 1e-8)
                    row["R9"] = {k: r9[k] for k in ("peak_GT_fine", "max_GT_fine", "peak_ImTbar", "max_mImTbar")}
                    row["rel_max_GT"] = row["max_GT_fine"] / r9["max_GT_fine"] - 1
                    row["rel_max_ImTbar"] = row["max_mImTbar"] / r9["max_mImTbar"] - 1
                    if "ldos_nn" in row:
                        c3 = C9[str(nk)]["brut"]
                        cmp_["ldos_at"] = abs(row["ldos_nn"]["at_eV"] - c3["at_eV"]) < 1e-6
                        cmp_["ldos_value"] = abs(row["ldos_nn"]["value"] / c3["value"] - 1) < 1e-8
                        row["R9_C3"] = c3
                    row["checks"] = cmp_; p4ok &= all(cmp_.values())
                p4[f"{S} {var} {nk}"] = row
                np.savez(os.path.join(d, f"res_{S}_{var}_nk{nk}.npz"), egrid_rel=st["egrid"] - st["E_D"], curve=o["curve"], tr=o["tr"],
                         prov=json.dumps(provenance(S, variant=var, N_k_int=NK(nk), R_cut=3, eta_t=eta_p)))
                rlog(f"[gate] P4 {S} {var} {nk}^2 : pic Gamma_T {row['peak_GT_fine']:+.4f} (max {row['max_GT_fine']:.4f}), pic -Im Tbar(K) "
                     f"{row['peak_ImTbar']:+.4f} (max {row['max_mImTbar']:.4f})" + (f" ; LDOS voisins {row['ldos_nn']['at_eV']:+.4f} "
                     f"({row['ldos_nn']['value']:.4f})" if "ldos_nn" in row else "") + ("" if SMOKE else f" ; R9 : {row['checks']}"))
            del g0
    ok["P4"] = bool(p4ok) if not SMOKE else True; res["P4"] = p4
    # ---- information (not gates): off-diagonal / C3 / K' / hermiticity of Tbar, t_pi vs t_5WF[pi, pi]
    info = {}
    for S in SIZES:
        for var in VARIANTS:
            ch = t_pi(S, 3, NK_MAIN, eta_p, var); P = ch["P"]; eg = ch["eg"] - ch["E_D"]
            TK = m.da.tbar_k(ch["tau"], ch["Du"], K[None])[:, 0]; TKp = m.da.tbar_k(ch["tau"], ch["Du"], Kp[None])[:, 0]
            fro = np.linalg.norm(TK, axis=(1, 2))
            k1 = K + np.array([0.01, 0.0, 0.0]); k2 = rot120(P["Bc"], k1); kk = np.vstack([k1, k2])   # C3 image of K + delta
            Hp = hk_pi(P, kk); _, Up = np.linalg.eigh(Hp)
            Tn = m.da.tbar_k(ch["tau"], ch["Du"], kk, U=Up)
            d12 = np.abs(np.diagonal(Tn[:, 0], axis1=1, axis2=2) - np.diagonal(Tn[:, 1], axis1=1, axis2=2)).max(1) / \
                np.abs(np.diagonal(Tn[:, 0], axis1=1, axis2=2)).max(1)
            iD = int(np.argmin(np.abs(eg))); ipk = int(np.argmax(-0.5 * np.trace(TK, axis1=1, axis2=2).imag))
            q = dict(AB_over_fro=np.abs(TK[:, 0, 1]) / fro, BA_over_fro=np.abs(TK[:, 1, 0]) / fro,
                     K_minus_Kp_over_fro=np.linalg.norm(TK - TKp, axis=(1, 2)) / fro,
                     herm_over_fro=np.linalg.norm(TK - np.conj(np.transpose(TK, (0, 2, 1))), axis=(1, 2)) / fro,
                     C3_bandiag_rel=d12)
            info[f"{S} {var}"] = {k_: dict(max=float(v.max()), at_E_D=float(v[iD]), at_peak_ImTbarK=float(v[ipk])) for k_, v in q.items()}
            info[f"{S} {var}"]["k_C3"] = kk; info[f"{S} {var}"]["peak_ImTbarK_pi_eV"] = float(eg[ipk])
            if S == "9x9" and var == "tel_quel":
                g, _ = g0_for(S, 3, NK_MAIN, eta_p, "win"); tf = m.local_t_cache(P["V_rc3"], g); del g
                ip = pi_index(len(P["Rloc_rc3"])); tfp = tf[:, ip][:, :, ip]
                info["t_pi_vs_t5WF_pipi_rel"] = float(np.abs(ch["t"] - tfp).max() / np.abs(tfp).max()); del tf, tfp
    res["information"] = info
    rlog(f"[gate] information : { {k: {kk: v[kk]['max'] for kk in v if isinstance(v[kk], dict)} if isinstance(v, dict) else v for k, v in info.items()} }")
    verdict = all(ok.values())
    res.update(ok=ok, verdict="PASS" if verdict else "ECHEC", elapsed_s=time.time() - t0)
    save_json(os.path.join(d, "gate_results.json"), res)
    rlog(f"[gate] verdict {res['verdict']} {ok} ({(time.time() - t0) / 60:.1f} min)")
    if not verdict:
        raise SystemExit(3)


# ------------------------------------------------------------------------------------------------ dos (steps 2, 2 bis)
def cmd_dos(a):
    m = mods(); c = cfg(); eta_p = c["eta_eV"]; d = ensure(OUTC, "dos"); t0 = time.time()
    S = "9x9"; P = load_prep(S); E_D = float(P["E_D"]); egrid, _ = egrid_res(E_D); eg = egrid[grid_index(egrid, E_D, "win")]; egr = eg - E_D
    grids = sorted(set(GRIDS_DOS) | set(GRIDS_PRIS))
    kgrid = {N: m.lt.mp_grid(NK(N), NK(N), 1) for N in grids}
    Hk = {N: hk_pi(P, kgrid[N]) for N in grids}
    curves = dict(eg_rel=egr); res = dict(prov=provenance(S), metrics={}, ldos={})
    tau0 = np.zeros((len(eg), 1, 2, 2), complex); Du0 = np.zeros((1, 3), int)
    for N in GRIDS_PRIS:
        curves[f"pristine_eta15_{N}"] = dos_par(Hk[N], tau0, Du0, kgrid[N], 0.0, eg, ETA_PRIS)["rho0"]
        rlog(f"[dos] pristine eta 15 meV, {NK(N)}^2 : rho(E_D) {np.interp(0, egr, curves[f'pristine_eta15_{N}']):.5f}, "
             f"integrale sur +-{WIN} eV {np.trapz(curves[f'pristine_eta15_{N}'], egr):.5f}")
    runs = [(var, NK_MAIN, eta_p, var) for var in VARIANTS] + [("eta_unique", 600, ETA_DOS, "tel_quel")]
    for name, nk, eta_t, var in runs:
        ch = t_pi(S, 3, nk, eta_t, var)
        for N in GRIDS_DOS:
            r = dos_par(Hk[N], ch["tau"], ch["Du"], kgrid[N], list(C_LIST), ch["eg"], ETA_DOS)
            key = f"{name}_{N}"; curves[f"{key}_rho0"] = r["rho0"]
            for ic, cv in enumerate(C_LIST):
                rho = r["rho"][ic]; curves[f"{key}_c{cv}"] = rho
                curves[f"{key}_c{cv}_protocole"] = curves["pristine_eta15_300"] + (rho - r["rho0"])
                mt = dmetrics(egr, rho, r["rho0"]); mt.update(variant=name, N_k_int=NK(nk), eta_t=eta_t, eta_G=ETA_DOS, grid=NK(N), c_i=cv,
                                                            C_N_eV=(C_N(S) if var == "aligne" else 0.0))
                res["metrics"][f"{key}_c{cv}"] = mt
                rlog(f"[dos] {name} {NK(N)}^2 c_i {cv:.3%} : max rho-rho0 {mt['drho_max_pos_eV']:+.4f} eV ({mt['drho_max']:.5f}), "
                     f"FWHM {mt['drho_fwhm_eV']:.4f}, rho(E_D) {mt['rho_at_ED']:.5f}")
    # 2 bis : LDOS of the three neighbours (cluster_ldos), N_k^int 900, both sizes and variants (« res » grid)
    for S2 in SIZES:
        P2 = load_prep(S2); g0, egr2 = g0_for(S2, 3, NK_MAIN, eta_p, "res")
        for var in VARIANTS:
            r, rho, rho0 = ldos_nn(P2, V_of(P2, 3, var, S2), g0, egr2, float(P2["E_D"]))
            res["ldos"][f"{S2} {var}"] = dict(r, nn_cells=P2["nn_cells"], neighbour_wf=int(7 - P2["vac_wf"]))
            curves[f"ldos_{S2}_{var}"] = rho.sum(1); curves[f"ldos0_{S2}_{var}"] = rho0.sum(1); curves[f"ldos_eg_rel_{S2}"] = egr2 - float(P2["E_D"])
            rlog(f"[dos] 2 bis LDOS voisins {S2} {var} : max {r['at_eV']:+.4f} eV ({r['value']:.4f} etats/eV)")
        del g0
    np.savez(os.path.join(d, "dos_9x9.npz"), **curves, prov=json.dumps(res["prov"]))
    res["elapsed_s"] = time.time() - t0
    save_json(os.path.join(d, "dos_results.json"), res)
    rlog(f"[dos] termine ({(time.time() - t0) / 60:.1f} min)")


# ------------------------------------------------------------------------------------------------ spec (step 3)
def branch_table(pk, s, egmax=WIN):
    """Two highest maxima of A_k at each k (by height), sorted by energy -> lower / upper branch; jumps > 0.15 eV between neighbours."""
    low, up = np.full(len(pk), np.nan), np.full(len(pk), np.nan)
    for i, p in enumerate(pk):
        p = p[np.abs(p[:, 0]) <= egmax]
        if len(p) >= 2:
            top = p[np.argsort(-p[:, 1])[:2]]; e = np.sort(top[:, 0]); low[i], up[i] = e
    jumps = {}
    for name, b in (("lower", low), ("upper", up)):
        jumps[name] = [dict(i=[int(i), int(i + 1)], s=[float(s[i]), float(s[i + 1])], eps_eV=[float(b[i]), float(b[i + 1])])
                       for i in range(len(b) - 1) if np.isfinite(b[i]) and np.isfinite(b[i + 1]) and abs(b[i + 1] - b[i]) > 0.15]
    return low, up, jumps


def cmd_spec(a):
    m = mods(); c = cfg(); eta_p = c["eta_eV"]; d = ensure(OUTC, "spec"); t0 = time.time()
    S = "9x9"; P = load_prep(S); E_D = float(P["E_D"]); K = np.array(c["K_red"], float)
    kp, labels, idx, s = m.build_k_path([("Γ", [0, 0, 0]), ("K", K), ("M", [0.5, 0, 0])], NPATH, P["Bc"])
    Hp = hk_pi(P, kp); Ep = np.linalg.eigvalsh(Hp) - E_D
    arrays = dict(k=kp, s=s, idx=np.array(idx), E_pristine=Ep); res = dict(prov=provenance(S, path="Gamma-K-M", n_k=NPATH), runs={})
    for name, nk, eta_t, var in [(v, NK_MAIN, eta_p, v) for v in VARIANTS] + [("eta_unique", 600, ETA_AK, "tel_quel")]:
        ch = t_pi(S, 3, nk, eta_t, var); egr = ch["eg"] - E_D
        Tb = m.da.tbar_k(ch["tau"], ch["Du"], kp)
        for cv in C_LIST:
            A = m.da.spectral_path(Hp, Tb, cv, ch["eg"], ETA_AK)
            pk = m.da.spectral_maxima(A, egr, prominence=1e-3 * A.max())
            iK = idx[1]; aK = pk[iK]; nearK = aK[np.abs(aK[:, 0]) <= 0.5]
            top = nearK[np.argsort(-nearK[:, 1])[:2]] if len(nearK) >= 2 else nearK
            top = top[np.argsort(top[:, 0])]
            low, up, jumps = branch_table(pk, s)
            nearK_j = {br: [j for j in v if max(abs(x - s[iK]) for x in j["s"]) <= 0.25] for br, v in jumps.items()}   # |k - K| <= 0.25 1/A
            key = f"{name}_c{cv}"
            res["runs"][key] = dict(variant=name, N_k_int=NK(nk), eta_t=eta_t, eta_G=ETA_AK, c_i=cv, maxima_at_K=aK, K_two_highest=top,
                                    K_gap_eV=float(top[1, 0] - top[0, 0]) if len(top) == 2 else None, jumps=jumps, jumps_near_K=nearK_j,
                                    maxima=[p for p in pk])
            arrays[f"A_{key}"] = A.astype(np.float32); arrays[f"low_{key}"] = low; arrays[f"up_{key}"] = up
            rlog(f"[spec] {key} : maxima a K dans +-0,5 eV {np.round(nearK, 4).tolist()} ; ecart des deux plus hauts "
                 f"{res['runs'][key]['K_gap_eV']} eV ; sauts > 0,15 eV a |k-K| <= 0,25 1/A : {nearK_j} ({sum(len(v) for v in jumps.values())} sur tout le chemin)")
        arrays["eg_rel"] = egr
    np.savez(os.path.join(d, "spectral_GKM_9x9.npz"), **arrays, prov=json.dumps(res["prov"]))
    res["elapsed_s"] = time.time() - t0
    save_json(os.path.join(d, "spec_results.json"), res)
    rlog(f"[spec] termine ({(time.time() - t0) / 60:.1f} min)")


# ------------------------------------------------------------------------------------------------ sens (step 5)
def cmd_sens(a):
    m = mods(); c = cfg(); eta_p = c["eta_eV"]; d = ensure(OUTC, "sens"); t0 = time.time()
    kg = m.lt.mp_grid(NK(GRID_SENS), NK(GRID_SENS), 1); Hk = {}; rows = []; curves = {}

    def one(group, S, rc, nk, var, etaG):
        ch = t_pi(S, rc, nk, eta_p, var); P = ch["P"]
        if S not in Hk:
            Hk[S] = hk_pi(P, kg)
        r = dos_par(Hk[S], ch["tau"], ch["Du"], kg, [0.01], ch["eg"], etaG)
        egr = ch["eg"] - ch["E_D"]; mt = dmetrics(egr, r["rho"][0], r["rho0"])
        mt.update(group=group, size=S, R_cut=rc, N_k_int=NK(nk), variant=var, eta_t=eta_p, eta_G=etaG, grid=NK(GRID_SENS), c_i=0.01,
                  C_N_eV=(C_N(S) if var == "aligne" else 0.0), cells_at_box_edge=int((P[f"edge_rc{rc}"] == 0).sum()),
                  cells_out_of_box=int((~P[f"inbox_rc{rc}"]).sum()), E_D=ch["E_D"])
        rows.append(mt); key = f"{group}_{S}_rc{rc}_nk{nk}_{var}_eta{int(etaG * 1e3)}"
        curves[key] = r["rho"][0] - r["rho0"]; curves[f"eg_{S}"] = egr
        rlog(f"[sens] {key} : max rho-rho0 {mt['drho_max_pos_eV']:+.4f} eV ({mt['drho_max']:.5f}), FWHM {mt['drho_fwhm_eV']:.4f}")
    for var in VARIANTS:
        for nk in NK_ALL:
            one("N_k_int", "9x9", 3, nk, var, ETA_DOS)
        for rc in (2, 4):
            one("R_cut", "9x9", rc, NK_MAIN, var, ETA_DOS)
        one("taille", "12x12", 3, NK_MAIN, var, ETA_DOS)
        one("eta_G", "9x9", 3, NK_MAIN, var, ETA_AK)
    np.savez(os.path.join(d, "sens_dos.npz"), **curves, prov=json.dumps(provenance()))
    save_json(os.path.join(d, "sens_results.json"), dict(prov=provenance(), rows=rows, elapsed_s=time.time() - t0))
    rlog(f"[sens] termine ({(time.time() - t0) / 60:.1f} min)")


# ------------------------------------------------------------------------------------------------ fig (7b, 2 bis table, figures)
def cmd_fig(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.join(PROJ, "figures", "memoire.mplstyle"))
    sys.path.insert(0, os.path.join(PROJ, "scripts"))
    from _palette import NAVY, ORANGE, GREEN, SKY, REF, MUTED, CMAP_SEQ
    d = ensure(OUTC, "fig"); fdir = ensure(HERE, "fig_smoke" if SMOKE else "fig")
    D = np.load(os.path.join(OUTC, "dos", "dos_9x9.npz")); Dj = json.load(open(os.path.join(OUTC, "dos", "dos_results.json")))
    Sp = np.load(os.path.join(OUTC, "spec", "spectral_GKM_9x9.npz")); Sj = json.load(open(os.path.join(OUTC, "spec", "spec_results.json")))
    Se = json.load(open(os.path.join(OUTC, "sens", "sens_results.json"))); Gj = json.load(open(os.path.join(OUTC, "gate", "gate_results.json")))
    Kx = np.loadtxt(os.path.join(OUTC.replace("out_smoke", "out"), "7a", "fig13_VA_interp_1meV.csv"), delimiter=",", skiprows=1)
    eK, rK = Kx[:, 0], {"pristine": Kx[:, 1], 0.001: Kx[:, 2], 0.01: Kx[:, 3]}
    K7 = json.load(open(os.path.join(OUTC.replace("out_smoke", "out"), "7a", "7a_results.json")))["fig13_metrics"]
    egr = D["eg_rel"]; col = {"tel_quel": NAVY, "aligne": ORANGE, "eta_unique": GREEN}
    lab = {"tel_quel": "tel quel", "aligne": r"aligné ($C_N$ plateau)", "eta_unique": r"$\eta_t = \eta_G$"}
    Nmain = NK(600)
    # ---- 7b : our curves at the protocol of the article vs the extracted curves
    sup = {}
    mK = (eK >= -1.0) & (eK <= 1.0)
    for name in ("tel_quel", "aligne", "eta_unique"):
        for N in GRIDS_DOS:
            for cv in C_LIST:
                key = f"{name}_{N}_c{cv}"
                ours = np.interp(eK, egr, D[f"{key}_protocole"])
                row = dict(pos_ours=Dj["metrics"][key]["drho_max_pos_eV"], pos_K=K7["0.1" if cv == 0.001 else "1.0"]["drho_max_pos_eV"])
                row["pos_diff_eV"] = row["pos_ours"] - row["pos_K"]
                for f_ in (1, 2):
                    row[f"rms_x{f_}"] = float(np.sqrt(np.mean((f_ * ours[mK] - rK[cv][mK]) ** 2)))
                sup[key] = row
    # ---- 2 bis table
    tab2 = {}
    for S in SIZES:
        for var in VARIANTS:
            g = Gj["P4"][f"{S} {var} {NK_MAIN}"]; L = Dj["ldos"][f"{S} {var}"]
            r = dict(ldos_nn_max_eV=L["at_eV"], ldos_nn_value=L["value"], peak_ImTbarK_eV=g["peak_ImTbar"], peak_GammaT_curve_eV=g["peak_GT_fine"])
            if S == "9x9":
                for cv in C_LIST:
                    r[f"drho_max_c{cv}_600"] = Dj["metrics"][f"{var}_{GRID_SENS}_c{cv}"]["drho_max_pos_eV"]
                    r[f"drho_max_c{cv}_300"] = Dj["metrics"][f"{var}_300_c{cv}"]["drho_max_pos_eV"]
            else:
                rr = [x for x in Se["rows"] if x["group"] == "taille" and x["variant"] == var][0]
                r["drho_max_c0.01_600"] = rr["drho_max_pos_eV"]
            tab2[f"{S} {var}"] = r
    save_json(os.path.join(d, "fig_results.json"), dict(prov=provenance(), superposition_7b=sup, table_2bis=tab2))
    # ---- figure dos_c
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.5), sharex=True)
    ax[0].plot(egr, D["pristine_eta15_300"], color=REF, lw=0.8, label=r"parfait ($\eta$ 15 meV)")
    for cv, cc in ((0.001, SKY), (0.01, NAVY)):
        ax[0].plot(egr, D[f"tel_quel_{GRID_SENS}_c{cv}_protocole"], color=cc, lw=0.9, label=rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%")
    for i, cv in enumerate(C_LIST):
        for name in ("tel_quel", "aligne", "eta_unique"):
            ax[i + 1].plot(egr, D[f"{name}_{GRID_SENS}_c{cv}"] - D[f"{name}_{GRID_SENS}_rho0"], color=col[name], lw=0.9, label=lab[name])
            ax[i + 1].plot(egr, D[f"{name}_300_c{cv}"] - D[f"{name}_300_rho0"], color=col[name], lw=0.5, ls=":")
        ax[i + 1].set_title(rf"({'bc'[i]}) $\rho - \rho_0$, $c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%", fontsize=9)
    ax[0].set_title(r"(a) $\rho(\varepsilon)$, tel quel", fontsize=9)
    ax[0].set_ylabel(r"DOS (états/eV/maille/spin)"); ax[1].legend(fontsize=6); ax[0].legend(fontsize=6)
    for x in ax:
        x.set_xlabel(r"Énergie $\varepsilon - E_D$ (eV)"); x.set_xlim(-1.2, 1.2)
    fig.tight_layout(); [fig.savefig(os.path.join(fdir, f"dos_c.{e}"), dpi=300) for e in ("pdf", "png")]; plt.close(fig)
    # ---- figure spectral_GKM
    names = ("tel_quel", "aligne", "eta_unique")
    fig, ax = plt.subplots(2, 3, figsize=(6.5, 4.4), sharex=True, sharey=True)
    s, idx, ege = Sp["s"], Sp["idx"], Sp["eg_rel"]
    for i, cv in enumerate(C_LIST):
        for j, name in enumerate(names):
            A = Sp[f"A_{name}_c{cv}"]; x = ax[i, j]
            x.imshow(A, origin="lower", aspect="auto", extent=(s[0], s[-1], ege[0], ege[-1]), cmap=CMAP_SEQ, vmin=0, vmax=np.percentile(A, 99.5))
            for br in ("low", "up"):                                 # plotted only where A >= 1 % of the map maximum (json keeps all)
                e_ = Sp[f"{br}_{name}_c{cv}"]; ok_ = np.isfinite(e_)
                h_ = np.zeros_like(e_); h_[ok_] = A[np.abs(ege[:, None] - e_[None, ok_]).argmin(0), np.where(ok_)[0]]
                vis = ok_ & (h_ >= 0.01 * A.max()); x.plot(s[vis], e_[vis], ".", ms=0.8, color=ORANGE)
            x.plot(s, Sp["E_pristine"], color=REF, lw=0.4, ls="--")
            x.set_xticks(s[idx]); x.set_xticklabels([r"$\Gamma$", r"$K$", r"$M$"]); x.set_ylim(-1.2, 1.2)
            x.set_title(f"({'abcdef'[3 * i + j]}) {lab[name]}, " + rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%", fontsize=8)
        ax[i, 0].set_ylabel(r"$\varepsilon - E_D$ (eV)")
    fig.tight_layout(); [fig.savefig(os.path.join(fdir, f"spectral_GKM.{e}"), dpi=300) for e in ("pdf", "png")]; plt.close(fig)
    # ---- figure sensibilites
    fig, ax = plt.subplots(1, 4, figsize=(6.5, 2.3), sharey=True)
    groups = (("N_k_int", "N_k_int", r"$N_k^\mathrm{int}$"), ("R_cut", "R_cut", r"$R_\mathrm{cut}$"), ("taille", "size", "taille"),
              ("eta_G", "eta_G", r"$\eta_G$ (eV)"))
    for x, (g, fld, xl) in zip(ax, groups):
        for var in VARIANTS:
            rows = [r for r in Se["rows"] if (r["group"] == g or (g == "R_cut" and r["group"] == "N_k_int" and r["N_k_int"] == NK(NK_MAIN))
                                               or (g in ("taille", "eta_G") and r["group"] == "N_k_int" and r["N_k_int"] == NK(NK_MAIN)))
                    and r["variant"] == var]
            xs = [r[fld] if fld != "size" else (9 if r["size"] == "9x9" else 12) for r in rows]
            o = np.argsort(xs); x.plot(np.array(xs)[o], np.array([r["drho_max_pos_eV"] for r in rows])[o], "o-", ms=3, lw=0.8, color=col[var],
                                       label=lab[var])
        x.axhline(K7["1.0"]["drho_max_pos_eV"], color=MUTED, lw=0.6, ls="--")
        x.set_xlabel(xl)
    ax[0].set_ylabel(r"max de $\rho - \rho_0$ (eV)"); ax[0].legend(fontsize=6)
    fig.tight_layout(); [fig.savefig(os.path.join(fdir, f"sensibilites.{e}"), dpi=300) for e in ("pdf", "png")]; plt.close(fig)
    # ---- figure superposition
    fig, ax = plt.subplots(2, 2, figsize=(6.5, 4.2), sharex=True)
    for i, f_ in enumerate((1, 2)):
        for j, cv in enumerate(C_LIST):
            x = ax[i, j]
            x.plot(eK, rK["pristine"], color=MUTED, lw=0.6, ls="--", label="Kaasbjerg, parfait")
            x.plot(eK, rK[cv], color=REF, lw=1.2, label="Kaasbjerg")
            for name in ("tel_quel", "aligne", "eta_unique"):
                x.plot(egr, f_ * D[f"{name}_{GRID_SENS}_c{cv}_protocole"], color=col[name], lw=0.8, label=lab[name])
            x.set_xlim(-1.2, 1.2)
            x.set_title(f"({'abcd'[2 * i + j]}) " + rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%, " + ("par spin" if f_ == 1 else r"$\times 2$ (spin compté)"),
                        fontsize=8)
        ax[i, 0].set_ylabel(r"DOS (eV$^{-1}$)")
    for x in ax[1]:
        x.set_xlabel(r"Énergie $\varepsilon - E_D$ (eV)")
    ax[0, 0].legend(fontsize=6)
    fig.tight_layout(); [fig.savefig(os.path.join(fdir, f"superposition.{e}"), dpi=300) for e in ("pdf", "png")]; plt.close(fig)
    rlog(f"[fig] figures dos_c, spectral_GKM, sensibilites, superposition -> {fdir} ; 7b : { {k: (round(v['pos_diff_eV'], 4), round(v['rms_x1'], 5), round(v['rms_x2'], 5)) for k, v in sup.items()} }")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name, h in (("extract", "step 7a (login node)"), ("prep", "M2 -> V_loc caches"), ("g0", "all g0 (16 BLAS threads)"),
                    ("gate", "step 1"), ("dos", "steps 2, 2 bis"), ("spec", "step 3"), ("sens", "step 5"), ("fig", "7b, figures, tables")):
        sp.add_parser(name, help=h)
    a = ap.parse_args()
    {"extract": cmd_extract, "prep": cmd_prep, "g0": cmd_g0, "gate": cmd_gate, "dos": cmd_dos, "spec": cmd_spec, "sens": cmd_sens,
     "fig": cmd_fig}[a.cmd](a)


if __name__ == "__main__":
    main()
