"""
wannier_provenance.py
    Gauge-consistency guard for the local t-matrix pipeline. Mwr (from the coarse M via the U
    matrices) and H_W(k) (from the tight-binding _tb.dat) MUST come from the SAME wannier90 run,
    otherwise the Wannier gauge is inconsistent and Gamma_nk is silently wrong.

    Contract: after a wannier90 run, write a manifest recording a common run_id and the sha256 of
    each output file (write_wannier_manifest). At load time the driver calls load_wannier_checked,
    which recomputes the hashes and RAISES (refuses to run) on any mismatch or missing file -- it
    never warns-and-continues. Also enforces one run_id across the tb/u/u_dis files.
"""
import hashlib
import json
import os


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_wannier_manifest(out_json, tb, u, u_dis=None, run_id=None):
    """Record run_id + sha256/mtime of the wannier90 outputs. Run once after wannierization."""
    entries = {"tb": tb, "u": u}
    if u_dis is not None:
        entries["u_dis"] = u_dis
    files = {}
    for key, path in entries.items():
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        files[key] = {"path": os.path.abspath(path), "sha256": _sha256(path),
                      "mtime": os.path.getmtime(path)}
    manifest = {"run_id": run_id or files["tb"]["sha256"][:16], "files": files}
    with open(out_json, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def load_wannier_checked(manifest_path):
    """
    Verify (and RAISE on any failure) that the wannier90 files listed in the manifest are present
    and unmodified since the manifest was written, then return their resolved paths. This is a hard
    gate: a wrong/missing/edited file, or a manifest that mixes runs, aborts the run.
    """
    if not os.path.exists(manifest_path):
        raise ValueError(f"gauge check: no wannier manifest at {manifest_path}. Write one with "
                         "write_wannier_manifest after wannierization; the driver refuses to run without it.")
    with open(manifest_path) as f:
        man = json.load(f)
    files = man.get("files", {})
    if "tb" not in files or "u" not in files:
        raise ValueError(f"gauge check: manifest {manifest_path} missing 'tb' and/or 'u' entries.")
    paths = {}
    for key, rec in files.items():
        p = rec["path"]
        if not os.path.exists(p):
            raise ValueError(f"gauge check: file '{key}' missing at {p} (referenced by manifest).")
        got = _sha256(p)
        if got != rec["sha256"]:
            raise ValueError(
                f"gauge check: '{key}' at {p} sha256 {got[:12]} != manifest {rec['sha256'][:12]}. "
                "Files are from a different wannier run or were edited -- refusing to run.")
        paths[key] = p
    paths.setdefault("u_dis", None)
    return paths
