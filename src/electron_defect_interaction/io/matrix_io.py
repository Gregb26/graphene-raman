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
M_NORM_V2 = "v2"             # R6 (2026-09-25): M^L and M^NL both in unit-cell Bloch norm (M2 = N_cells M^L + M^NL)
M_NORM_V1 = "v1"             # pre-R6 files (results/M, frozen): M^L in supercell norm, factor N_cells missing
HARTREE = "hartree"          # the only unit M files are ever stored in (kernels: V_ed in Ha, psi dimensionless)
EV = "eV"
HA2EV = 27.211386245988


def _manifest_path(npy_path):
    return os.path.splitext(npy_path)[0] + ".json"


def save_M(npy_path, M, bloch_norm, units=HARTREE, **extra):
    """Save M to npy_path and a sidecar <stem>.json recording bloch_norm and units (+ any extra metadata).
    M must be in Hartree (the kernels' unit); saving in any other unit is refused."""
    if units != HARTREE:
        raise ValueError(f"save_M: M files are stored in Hartree only (got units='{units}')")
    np.save(npy_path, M)
    meta = {"bloch_norm": bloch_norm, "units": units, "shape": list(np.shape(M)), **extra}
    with open(_manifest_path(npy_path), "w") as f:
        json.dump(meta, f, indent=2)


def read_manifest(npy_path):
    p = _manifest_path(npy_path)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def load_M_checked(npy_path, require_bloch_norm=SUPERCELL, units=None, require_normalization=None):
    """
    Load M, refusing (ValueError) unless its manifest declares require_bloch_norm AND units='hartree'
    AND (when require_normalization is given, e.g. matrix_io.M_NORM_V2 from config["M_normalization"][:2])
    the sidecar key 'M_normalization' starts with that version; a sidecar without the key is a v1 file.
    `units` (mandatory) is the unit the CALLER wants back:
        units='hartree' -> raw array (M construction / checks);
        units='eV'      -> M * HA2EV, the ONLY place this conversion may happen (t-matrix consumers:
                           the Wannier Hamiltonian and all energies are eV).
    Prevents un-normalized/double-normalized M and Ha-vs-eV mixing (the 2026-09-05 bug).
    """
    if units not in (HARTREE, EV):
        raise ValueError(f"load_M_checked: units must be '{HARTREE}' or '{EV}' (explicit), got {units!r}")
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
    stored = meta.get("units")
    if stored != HARTREE:
        raise ValueError(
            f"{npy_path}: manifest units={stored!r}, expected '{HARTREE}'. Refusing to run: unknown unit "
            f"would silently mix Hartree and eV in the t-matrix. Tag the sidecar (units='hartree') only if "
            f"the file was produced by the Ha kernels.")
    if require_normalization is not None:
        ver = str(meta.get("M_normalization", M_NORM_V1))
        if not ver.startswith(require_normalization):
            raise ValueError(
                f"{npy_path}: M_normalization={ver!r} but {require_normalization!r} required. Refusing to run "
                f"(pre-R6 file: M^L in supercell norm, factor N_cells missing; use results/M2).")
    M = np.load(npy_path)
    return M * HA2EV if units == EV else M
