import os, subprocess, hashlib, platform
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np
from graphene_raman.config import ROOT, HA2EV, wannier_dir, load_production, twied_dir
from graphene_raman.io import qe_io, pseudo_io, wannier_io
from graphene_raman.wannier.wannier_interpolation import _match_kpoint_order
import json

# paths
EX = os.path.join(ROOT, "data/export_chaine_5x5/data/graphene")
INPUTS = {
    "unit_cell":   f"{EX}/unit_cell/qe/defect_unit_cell_5x5_wann.save",
    "p": f"{EX}/supercell/qe/defect_5x5_p.save",
    "d": f"{EX}/supercell/qe/defect_5x5_d.save",
    "wannier":     wannier_dir(load_production(verbose=False), 5),
}
POT_PATH = {
    "p":          f"{INPUTS["supercell_p"]}/Vks_5x5_p",
    "d":          f"{INPUTS["supercell_d"]}/Vks_5x5_d"
}
PSEUDO_PATH = f"{INPUTS["unit_cell"]}/C.upf"

def _git(*args):
    """Run `git <args>` in the repository root and return its stdout, stripped."""
    res = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True)
    return res.stdout.strip()

def _put(g, name, data, units, axes=None):

    if np.asarray(data).nbytes > 1 << 20:
        ds = g.create_dataset(name, data=data, compression='gzip')
    else:
        ds = g.create_dataset(name, data=data)

    ds.attrs["units"] = units

    if axes is not None: 
        ds.attrs["axes"] = axes
    return ds

def md5_of_file(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def input_files(dirs):
    """Get all files in dirs in stable order"""
    files = []
    for d in dirs:
        d = Path(d)
        if not d.is_dir():
            raise FileNotFoundError(d)
        files += sorted(p for p in d.iterdir() if p.is_file())
    return files

        
def collect_provenance(input_dirs, allow_dirty=False):

    status = _git("status", "--porcelain", "src", "config", "campagnes/TWIED")
    git_clean = status == ""
    if not git_clean and not allow_dirty:
        raise RuntimeError(f"working tree not clean:\n{status}")
    git_commit = _git("rev-parse", "HEAD")

    # get date
    date = datetime.now().astimezone().isoformat(timespec="seconds")
    run_id = date
    code = "graphene-raman"
    script = os.path.relpath(__file__, ROOT)
    source_files = {os.path.relpath(p, ROOT): md5_of_file(p) for p in input_files(input_dirs)}

    return {"run_id":run_id, 
            "code":code, 
            "git_commit":git_commit, 
            "git_clean":git_clean, 
            "script":script, 
            "date":date, 
            "HA2EV":HA2EV, 
            "source_files": source_files, 
            "source_files_relative_to":"Root of graphene-raman",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "h5py": h5py.__version__,
            "hdf5":h5py.version.hdf5_version}

def write_provenance(f, prov, role):
    """
    f: opened HDF file
    prov: dict returned by collect_provenance()
    role: "inputs" or "chain"
    """
    if role not in ["inputs", "chain"]:
        raise ValueError(f"role must be 'inputs' or 'chain', got {role!r}")
    
    g = f.create_group("provenance")
    for key, value in prov.items():
        if key == "source_files":
            continue
        g.attrs[key] = value

    g.attrs["role"] = role

    sg = g.create_group("source_files")
    for key, value in prov["source_files"].items():
        sg.attrs[key] = value

def write_unit_cell(f, uc_save):
    g = f.create_group("input/unit_cell")

    C_nk, nG = qe_io.get_C_nk(uc_save)
    G_red = qe_io.get_G_red(uc_save)
    k_red = qe_io.get_k_red(uc_save)
    eps   = qe_io.get_eigenvalues(uc_save)
    A_cols, Omega = qe_io.get_A_volume(uc_save)
    x_red = qe_io.get_x_red(uc_save)
    ecut  = qe_io.get_ecut(uc_save)

    ds = _put(g, "C_nk", C_nk, "dimensionless", "band, k, G")
    ds.attrs["note"] = "zero-padded beyond nG[k]"
    _put(g, "nG", nG, "count", "k")
    ds = _put(g, "G_red", G_red, "reduced", "k, G, component")
    ds.attrs["not"] = "zero-padded beyond ng[k]"
    _put(g, "k_red", k_red, "reduced", "k, component")
    _put(g, "eps", eps, "hartree", "band, k")
    _put(g, "A_cols", A_cols, "bohr", "cartesian, i")
    _put(g, "Omega", Omega, "bohr^3")
    _put(g, "x_red", x_red, "reduced", "atom, component")
    _put(g, "ecut", ecut, "hartree")

    g.attrs["ngfft"] = qe_io.get_ngfft(uc_save)
    g.attrs["k_order"] = "QE"

def write_supercell(f, name, sc_save, pot_path):
    if name not in ["d", "p"]:
        raise ValueError(f"name must be 'd' or 'p', got {name!r}")
    
    g = f.create_group("input/supercell_"+name)

    A_cols, Omega = qe_io.get_A_volume(sc_save)
    _put(g, "A_cols", A_cols, "bohr", "cartesian, i")
    _put(g, "Omega", Omega, "bohr^3")
    
    x_red = qe_io.get_x_red(sc_save)
    _put(g, "x_red", x_red, "reduced", "atom, component")

    fp = qe_io.read_filplot(pot_path)
    _put(g, "V", fp["V"], "hartree", "ix, iy, iz")

    g.attrs["ngfft"] = fp["ngfft"]
    g.attrs["plot_num"] = fp["plot_num"]

def write_pseudo(f, upf_path):

    ekb_li, fr_li, rgrid, lmax, imax, V_L = pseudo_io.read_upf(upf_path)

    g = f.create_group("input/pseudo")

    _put(g, "ekb_li", ekb_li, "hartree", "l, i")
    _put(g, "fr_li", fr_li, "UPF PP_BETA convention (r*beta), not converted", "l, i, r")
    _put(g, "rgrid", rgrid, "bohr", "r")
    _put(g, "V_L", V_L, "hartree", "r")
    g.attrs["lmax"] = lmax
    g.attrs["imax"] = imax

def write_wannier(f, wdir, k_qe):
    U, k_w90  = wannier_io.read_w90_mat(f"{wdir}/wannier_u.mat")
    U_dis, _ = wannier_io.read_w90_mat(f"{wdir}/wannier_u_dis.mat")
    H_R, R, ndegen, r_R, lattice = wannier_io.read_w90_tb(f"{wdir}/wannier_tb.dat")

    g = f.create_group("input/wannier")
    ds = _put(g, "U", U, "dimensionless", "k, wannier, wannier")
    ds.attrs["note"] = "Wannier90 kpoint order"
    _put(g, "U_dis", U_dis, "dimensionless", "k, band, wannier")
    _put(g, "k_w90", k_w90, "reduced", "k, component")

    perm = _match_kpoint_order(k_w90, k_qe)
    ds = _put(g, "perm", perm, "index", "k")
    ds.attrs["note"] = "_match_kpoint_order(k_w90, k_qe)"

    ds = _put(g, "H_R", H_R, "eV", "R, wannier, wannier")
    ds.attrs["note"] = "NOT divided by ndegen"

    _put(g, "R", R, "reduced", "R, component")
    _put(g, "ndegen", ndegen, "count", "R")
    _put(g, "r_R", r_R, "angstrom", "R, cartesian, wannier, wannier")
    _put(g, "lattice", lattice, "angstrom", "cartesian, i")

    with open(f"{wdir}/wannier_manifest.json") as fh:
        man = json.load(fh)
        g.attrs["manifest"] = json.dumps(man)

def main():

    prov = collect_provenance(list(INPUTS.values()))

    cfg = load_production(verbose=False)
    out_dir = twied_dir(cfg)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "ref_5x5_inputs.h5")

    k_qe = qe_io.get_k_red(INPUTS["unit_cell"])

    with h5py.File(path, "w") as f:
        write_provenance(f, prov, "inputs")
        write_unit_cell(f, INPUTS['unit_cell'])
        write_supercell(f, INPUTS["p"], POT_PATH["p"])
        write_supercell(f, INPUTS["d"], POT_PATH["d"])
        write_pseudo(f, PSEUDO_PATH)
        write_wannier(f, INPUTS["wannier"], k_qe=k_qe)

if __name__ == "__main__":
    main()
