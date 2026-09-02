"""
matrix_io.py
    Save/load the electron-defect scattering matrix M together with a JSON sidecar manifest
    recording the Bloch-normalization convention. This guards against the most likely future
    bug: a script silently re-applying the 1/N_cells factor to an already-normalized M (or
    running on an un-normalized one). Consumers MUST use load_M_checked so a wrong/missing tag
    aborts the run instead of producing quietly wrong physics.

    Convention values for `bloch_norm`:
        "supercell"  -> M[n,k] built from Bloch states normalized over the supercell, i.e. the
                        matrix elements already carry the 1/N_cells factor (correct for the
                        single-defect T-matrix). This is what the physics code requires.
        "unit_cell"  -> raw matrix elements, O(1), MISSING the 1/N_cells factor.
"""
import json
import os

import numpy as np

SUPERCELL = "supercell"
UNIT_CELL = "unit_cell"


def _manifest_path(npy_path):
    return os.path.splitext(npy_path)[0] + ".json"


def save_M(npy_path, M, bloch_norm, **extra):
    """Save M to npy_path and a sidecar <stem>.json recording bloch_norm (+ any extra metadata)."""
    np.save(npy_path, M)
    meta = {"bloch_norm": bloch_norm, "shape": list(np.shape(M)), **extra}
    with open(_manifest_path(npy_path), "w") as f:
        json.dump(meta, f, indent=2)


def read_manifest(npy_path):
    p = _manifest_path(npy_path)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def load_M_checked(npy_path, require_bloch_norm=SUPERCELL):
    """
    Load M, refusing (ValueError) unless its manifest declares require_bloch_norm. Prevents both
    running on an un-normalized M and double-applying the 1/N_cells factor.
    """
    meta = read_manifest(npy_path)
    if meta is None:
        raise ValueError(
            f"{npy_path}: no manifest sidecar (.json), unknown Bloch normalization. Refusing to run; "
            f"expected bloch_norm='{require_bloch_norm}'. Regenerate/migrate M with matrix_io.save_M.")
    got = meta.get("bloch_norm")
    if got != require_bloch_norm:
        raise ValueError(
            f"{npy_path}: bloch_norm='{got}' but '{require_bloch_norm}' required. Refusing to run "
            f"(this M is likely un-normalized, or already normalized and about to be double-counted).")
    return np.load(npy_path)
