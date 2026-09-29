#!/usr/bin/env python
"""
r8_driver.py -- R8 (reproduction of Kaasbjerg, PRB 101, 045433 (2020), Figs. 13 and 14) pilot.

Sub-commands
    extract   step 7a: vector extraction of Fig. 13 (top panel, vacancies V_A) and reading of Fig. 14 (left column, V_A)
              from the article PDF with `mutool trace`; axis calibration from the tick marks, residuals, curve metrics.
              Runs on the login node (seconds). No computation of our DOS.

Outputs: out/7a/ (csv, json), fig/7a_controle.{pdf,png}. All energies in eV relative to the zero of the article's axes.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("extract", help="step 7a (login node)")
    a = ap.parse_args()
    {"extract": cmd_extract}[a.cmd](a)


if __name__ == "__main__":
    main()
