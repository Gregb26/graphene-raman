#!/usr/bin/env python
"""
assemble_M2.py -- R6 (2026-09-25): reassemble a production M with the corrected normalization of the local part.

    M2 = N_cells * M^L(v1) + M^NL          (both parts in unit-cell Bloch norm)

v1 = any M^L file produced before commit ee051ec (kernels compute_ML_R* with psi normalized over the SUPERCELL: the local
part is N_cells = Omega_sc/Omega_uc times too small relative to M^NL, R5-A.2). Nothing is recomputed and no input is ever
modified. Outputs are new files (refused if they exist, unless --force), each with a JSON sidecar carrying
bloch_norm=unit_cell, units=hartree, M_normalization="v2 ...", N_cells, assembled_from, date, md5:

    --out-l    N_cells * M^L(v1)   (or a copy of --ml-v2, an M^L already computed with the corrected kernel)
    --out-nl   M^NL: copy of --nl, symlink to --nl (--nl-link), or M_ed(v1) - M^L(v1) (--nl-from-diff; v1 assembly inverted
               exactly, the June/September files being plain sums)
    --out-m    M2 = out_l + out_nl

Blockwise over the bra-band axis (memmap) so 6.7 GB files never fill the memory. Checks, re-read from disk after writing:
max|out_m - (out_l + out_nl)| (must be 0), max|out_l - N_cells*ml| (must be 0), max|out_nl - source| (0), hermiticity
max|M2 - M2^dag| / max|M2|, max|.| of each part. One JSON summary line per call (prefix [assemble_M2]).
"""
import os
import sys
import json
import time
import hashlib
import argparse
import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(HERE)
try:
    import electron_defect_interaction  # noqa: F401
except ImportError:
    sys.path.insert(0, os.path.join(PROJ, "src"))
from electron_defect_interaction.io import matrix_io

M_NORM_V2 = "v2 : L et NL en norme unit_cell, 2026-09-25"


def md5sum(path, chunk=64 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sidecar(npy_path, M_shape, part, base, extra):
    meta = {k: v for k, v in (base or {}).items() if k not in ("shape", "note", "part", "N_kd", "N_cells", "M_normalization", "md5", "assembled_from", "date")}
    meta.update({"bloch_norm": matrix_io.UNIT_CELL, "units": matrix_io.HARTREE, "shape": [int(x) for x in M_shape], "part": part,
                 "M_normalization": M_NORM_V2, "campaign": "R6", "date": datetime.date.today().isoformat()})
    meta.update(extra)
    with open(os.path.splitext(npy_path)[0] + ".json", "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)


def refuse_existing(paths, force):
    for p in paths:
        if p and (os.path.lexists(p) or os.path.exists(os.path.splitext(p)[0] + ".json")):
            if not force:
                raise SystemExit(f"assemble_M2: {p} (or its sidecar) exists; refusing to overwrite (use --force)")
            for q in (p, os.path.splitext(p)[0] + ".json"):
                if os.path.lexists(q):
                    os.remove(q)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--ml", help="M^L v1 file (supercell-normalized local part): scaled by N_cells")
    g.add_argument("--ml-v2", help="M^L already in unit-cell norm (corrected kernel, commit ee051ec): copied, not scaled")
    ap.add_argument("--nl", help="M^NL file (unchanged by R6)")
    ap.add_argument("--nl-from-diff", metavar="M_ED_V1", help="M^NL = M_ed(v1) - M^L(v1); M^L(v1) is --ml, or --diff-ml, or --ml-v2 / N_cells")
    ap.add_argument("--diff-ml", help="explicit M^L(v1) for --nl-from-diff (e.g. a June M^L known to be wrong but summed into M_ed)")
    ap.add_argument("--nl-link", action="store_true", help="--out-nl is a symbolic link to --nl (sidecar written as a real file)")
    ap.add_argument("--n-cells", type=int, required=True, help="N_cells = Omega_sc / Omega_uc = N^2 (NOT N_kd = D^2 for the dense M)")
    ap.add_argument("--out-l", required=True); ap.add_argument("--out-nl"); ap.add_argument("--out-m")
    ap.add_argument("--block", type=int, default=0, help="bra bands per block (0 = auto, ~1.5 GB per array)")
    ap.add_argument("--label", default=""); ap.add_argument("--force", action="store_true"); ap.add_argument("--no-md5", action="store_true")
    ap.add_argument("--summary-file", help="append the JSON summary line to this file")
    a = ap.parse_args()
    t0 = time.time()
    N = float(a.n_cells)

    ml_path = a.ml or a.ml_v2
    ml = np.load(ml_path, mmap_mode="r")
    if ml.ndim != 4:
        raise SystemExit(f"{ml_path}: expected (nb, nk, nb, nk), got {ml.shape}")
    nb, nk = ml.shape[0], ml.shape[1]
    base = matrix_io.read_manifest(ml_path) or {}
    if base.get("bloch_norm") not in (None, matrix_io.UNIT_CELL) or base.get("units") not in (None, matrix_io.HARTREE):
        raise SystemExit(f"{ml_path}: unexpected sidecar {base}")
    if a.ml_v2 and not str(base.get("M_normalization", "")).startswith("v2"):
        raise SystemExit(f"{a.ml_v2}: --ml-v2 requires a sidecar tagged M_normalization=v2 (corrected kernel), got {base.get('M_normalization')!r}")
    if a.ml and str(base.get("M_normalization", "")).startswith("v2"):
        raise SystemExit(f"{a.ml}: already v2; use --ml-v2 (would double-apply N_cells)")

    nl = m_ed = diff_ml = None
    if a.out_nl or a.out_m:
        if a.nl:
            nl = np.load(a.nl, mmap_mode="r"); assert nl.shape == ml.shape, (nl.shape, ml.shape)
            nl_meta = matrix_io.read_manifest(a.nl)
            if nl_meta is not None and (nl_meta.get("bloch_norm") != matrix_io.UNIT_CELL or nl_meta.get("units") != matrix_io.HARTREE):
                raise SystemExit(f"{a.nl}: unexpected sidecar {nl_meta}")
        elif a.nl_from_diff:
            m_ed = np.load(a.nl_from_diff, mmap_mode="r"); assert m_ed.shape == ml.shape, (m_ed.shape, ml.shape)
            if a.diff_ml:
                diff_ml = np.load(a.diff_ml, mmap_mode="r"); assert diff_ml.shape == ml.shape
        else:
            raise SystemExit("--out-nl/--out-m need --nl or --nl-from-diff")
    if a.nl_link and not (a.nl and a.out_nl):
        raise SystemExit("--nl-link needs --nl and --out-nl")

    refuse_existing([a.out_l, a.out_nl, a.out_m], a.force)
    for p in (a.out_l, a.out_nl, a.out_m):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    per_band = nk * nb * nk * 16
    blk = a.block or max(1, min(nb, int(1.5e9 // per_band)))
    shape = tuple(int(x) for x in ml.shape)
    out_l = np.lib.format.open_memmap(a.out_l, mode="w+", dtype=np.complex128, shape=shape)
    out_nl = None if (not a.out_nl or a.nl_link) else np.lib.format.open_memmap(a.out_nl, mode="w+", dtype=np.complex128, shape=shape)
    out_m = np.lib.format.open_memmap(a.out_m, mode="w+", dtype=np.complex128, shape=shape) if a.out_m else None
    mx = dict(L=0.0, NL=0.0, M=0.0)

    def nl_block(b0, b1):
        if nl is not None:
            return np.asarray(nl[b0:b1])
        L1 = np.asarray(diff_ml[b0:b1]) if diff_ml is not None else (np.asarray(ml[b0:b1]) if a.ml else np.asarray(ml[b0:b1]) / N)
        return np.asarray(m_ed[b0:b1]) - L1

    for b0 in range(0, nb, blk):
        b1 = min(nb, b0 + blk)
        L2 = np.asarray(ml[b0:b1]) * N if a.ml else np.asarray(ml[b0:b1])
        out_l[b0:b1] = L2; mx["L"] = max(mx["L"], float(np.abs(L2).max()))
        if out_nl is not None or out_m is not None:
            NLb = nl_block(b0, b1); mx["NL"] = max(mx["NL"], float(np.abs(NLb).max()))
            if out_nl is not None:
                out_nl[b0:b1] = NLb
            if out_m is not None:
                Mb = L2 + NLb; out_m[b0:b1] = Mb; mx["M"] = max(mx["M"], float(np.abs(Mb).max()))
            del NLb
        del L2
    for arr in (out_l, out_nl, out_m):
        if arr is not None:
            arr.flush()
    del out_l, out_nl, out_m
    if a.nl_link:
        os.symlink(os.path.relpath(os.path.abspath(a.nl), os.path.dirname(os.path.abspath(a.out_nl))), a.out_nl)
    t_write = time.time() - t0

    # ---- verification, re-read from disk
    OL = np.load(a.out_l, mmap_mode="r"); ONL = np.load(a.out_nl, mmap_mode="r") if a.out_nl else None; OM = np.load(a.out_m, mmap_mode="r") if a.out_m else None
    chk = dict(max_dL_vs_source=0.0, max_dNL_vs_source=0.0, max_dM_vs_sum=0.0, herm_abs=0.0)
    for b0 in range(0, nb, blk):
        b1 = min(nb, b0 + blk)
        Lsrc = np.asarray(ml[b0:b1]) * N if a.ml else np.asarray(ml[b0:b1])
        chk["max_dL_vs_source"] = max(chk["max_dL_vs_source"], float(np.abs(np.asarray(OL[b0:b1]) - Lsrc).max()))
        if ONL is not None:
            chk["max_dNL_vs_source"] = max(chk["max_dNL_vs_source"], float(np.abs(np.asarray(ONL[b0:b1]) - nl_block(b0, b1)).max()))
        if OM is not None:
            Mb = np.asarray(OM[b0:b1])
            chk["max_dM_vs_sum"] = max(chk["max_dM_vs_sum"], float(np.abs(Mb - (np.asarray(OL[b0:b1]) + np.asarray(ONL[b0:b1]))).max()))
            R = Mb.reshape((b1 - b0) * nk, nb * nk)                                   # rows (bra band block)
            Cc = np.asarray(OM[:, :, b0:b1, :]).reshape(nb * nk, (b1 - b0) * nk)     # columns (ket band block)
            chk["herm_abs"] = max(chk["herm_abs"], float(np.abs(R - Cc.conj().T).max()))
            del Mb, R, Cc
        del Lsrc
    chk["herm_rel"] = chk["herm_abs"] / mx["M"] if mx["M"] > 0 else 0.0
    ok = chk["max_dL_vs_source"] == 0.0 and chk["max_dNL_vs_source"] == 0.0 and chk["max_dM_vs_sum"] == 0.0 and (OM is None or chk["herm_rel"] < 1e-12)

    md5 = {}
    if not a.no_md5:
        for key, p in (("out_l", a.out_l), ("out_nl", a.out_nl), ("out_m", a.out_m)):
            if p:
                md5[key] = md5sum(os.path.realpath(p))
    src = dict(ml=os.path.abspath(ml_path), ml_scale=(a.n_cells if a.ml else 1), ml_is_v2=bool(a.ml_v2),
               nl=(os.path.abspath(a.nl) if a.nl else None),
               nl_from_diff=(dict(m_ed=os.path.abspath(a.nl_from_diff), minus=os.path.abspath(a.diff_ml or ml_path) + ("" if (a.diff_ml or a.ml) else f" / {a.n_cells}")) if a.nl_from_diff else None),
               nl_link=bool(a.nl_link), label=a.label,
               ml_mtime=datetime.datetime.fromtimestamp(os.path.getmtime(ml_path)).isoformat(timespec="seconds"), ml_size=os.path.getsize(ml_path))
    part_base = str(base.get("part", "M_L")).replace("M_L", "")
    sidecar(a.out_l, shape, "M_L" + part_base, base, dict(N_cells=a.n_cells, assembled_from=src, md5=md5.get("out_l"), max_abs=mx["L"]))
    if a.out_nl:
        sidecar(a.out_nl, shape, "M_NL" + part_base, base, dict(N_cells=a.n_cells, assembled_from=src, md5=md5.get("out_nl"), max_abs=mx["NL"],
                                                                symlink_to=(os.readlink(a.out_nl) if a.nl_link else None)))
    if a.out_m:
        sidecar(a.out_m, shape, "M2" + part_base, base, dict(N_cells=a.n_cells, assembled_from=src, md5=md5.get("out_m"), max_abs=mx["M"],
                                                             checks=chk, N_kd=base.get("N_kd", nk)))
    summ = dict(out_l=a.out_l, out_nl=a.out_nl, out_m=a.out_m, n_cells=a.n_cells, shape=list(shape), max_abs=mx, checks=chk, md5=md5, ok=bool(ok),
                seconds=round(time.time() - t0, 1), seconds_write=round(t_write, 1), label=a.label, block=blk)
    line = "[assemble_M2] " + json.dumps(summ, ensure_ascii=False)
    print(line, flush=True)
    if a.summary_file:
        with open(a.summary_file, "a") as f:
            f.write(line + "\n")
    print(f"[assemble_M2] {a.label or os.path.basename(a.out_l)}: N_cells={a.n_cells} max|L2|={mx['L']:.4e} max|NL|={mx['NL']:.4e} max|M2|={mx['M']:.4e} Ha ; "
          f"dL={chk['max_dL_vs_source']:.1e} dNL={chk['max_dNL_vs_source']:.1e} dM={chk['max_dM_vs_sum']:.1e} herm_rel={chk['herm_rel']:.1e} -> {'OK' if ok else 'ECHEC'} ({time.time()-t0:.0f} s)", flush=True)
    if not ok:
        sys.exit(3)


if __name__ == "__main__":
    main()
