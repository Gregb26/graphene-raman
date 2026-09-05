"""Frozen production parameters: config/production.json (versioned). Single loader for every production script."""
import json
import os

HA2EV = 27.211386245988
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PRODUCTION_JSON = os.path.join(ROOT, "config", "production.json")


def load_production(path=PRODUCTION_JSON, verbose=True):
    with open(path) as f:
        cfg = json.load(f)
    for key in ("R_cut", "grid", "eta_eV", "N_min", "reference_size", "nk_int", "e_window_eV", "ne_per_eta"):
        if key not in cfg:
            raise KeyError(f"production config {path} lacks '{key}'")
    if verbose:
        print(f"[config] {os.path.relpath(path, ROOT)}: R_cut={cfg['R_cut']} grid={cfg['grid']} eta={cfg['eta_eV']} eV "
              f"N_min={cfg['N_min']} ref={cfg['reference_size']} nk_int={cfg['nk_int']}", flush=True)
    return cfg


def dense_paths(cfg, size, scratch="/home/gregb26/links/scratch/qe_tmp", root=ROOT):
    """Dense primitive .save, dense M file and dense Wannier dir for a supercell size."""
    D = cfg["dense"][size]["D"]
    return dict(D=D, p=cfg["dense"][size]["p"],
                uc=f"{scratch}/defect_uc_dense_{D}/defect_uc_dense_{D}.save",
                mfile=os.path.join(root, f"results/M/M_dense_{size}.npy"),
                wdir=os.path.join(root, f"wannier/{D}x{D}"),
                manifest=os.path.join(root, f"wannier/{D}x{D}/wannier_manifest.json"))
