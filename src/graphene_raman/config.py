"""Frozen production parameters: config/production.json (versioned). Single loader for every production script."""
import json
import os

HA2EV = 27.211386245988
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PRODUCTION_JSON = os.path.join(ROOT, "config", "production.json")


def load_production(path=PRODUCTION_JSON, verbose=True):
    with open(path) as f:
        cfg = json.load(f)
    for key in ("R_cut", "grid", "eta_eV", "N_min", "reference_size", "nk_int", "e_window_eV", "ne_per_eta",
                "M_normalization", "results_dir", "matrices_dir", "alignment"):
        if key not in cfg:
            raise KeyError(f"production config {path} lacks '{key}'")
    if verbose:
        print(f"[config] {os.path.relpath(path, ROOT)}: R_cut={cfg['R_cut']} grid={cfg['grid']} eta={cfg['eta_eV']} eV "
              f"N_min={cfg['N_min']} ref={cfg['reference_size']} nk_int={cfg['nk_int']} "
              f"M={cfg['M_normalization']} results={cfg['results_dir']} matrices={cfg['matrices_dir']} "
              f"C_N={cfg['alignment']['version']} ({len(cfg['alignment']['C_N_eV'])} sizes)", flush=True)
    return cfg


def results_dir(cfg, root=ROOT):
    """Directory of the production products only (results/M2_plateau since R10; results/M and results/M2 are frozen). M files: matrices_dir."""
    return os.path.join(root, cfg["results_dir"])


def matrices_dir(cfg, root=ROOT):
    """Directory of the M files (results/M2: raw M2 matrices, read only; R10)."""
    return os.path.join(root, cfg["matrices_dir"])


def epw_dir(cfg, root=ROOT):
    """Directory of the EPW post-processing products (chapter 5): results/epw."""
    return os.path.join(root, cfg.get("epw_results_dir", "results/epw"))


def electron_dir(cfg, root=ROOT):
    """Directory of the chapter 2 electronic-structure data (QE convergence .dat, bands, DOS): results/electron."""
    return os.path.join(root, cfg.get("electron_results_dir", "results/electron"))


def phonon_dir(cfg, root=ROOT):
    """Directory of the chapter 2 phonon data (matdyn dispersion and DOS extracted from EPW runs): results/phonon."""
    return os.path.join(root, cfg.get("phonon_results_dir", "results/phonon"))


def wannier_dir(cfg, D=None, root=ROOT):
    """Directory of the tracked Wannier90 outputs: results/wannier, or results/wannier/<D>x<D> for a grid (D int or "DxD")."""
    base = os.path.join(root, cfg.get("wannier_dir", "results/wannier"))
    if D is None:
        return base
    return os.path.join(base, D if isinstance(D, str) else f"{D}x{D}")


def campaigns_dir(cfg, *parts, root=ROOT):
    """Directory of the campaign data read by the figure scripts (json/npz copied from campagnes/, same relative paths):
    results/campagnes[/<campaign>/<sub>/...]."""
    return os.path.join(root, cfg.get("campaigns_results_dir", "results/campagnes"), *parts)


def alignment_C(cfg, size):
    """Far-field potential offset C_N (eV) of a supercell size (block "alignment", plateau (i), R10), subtracted as M_W(R,R) - C_N on the
    N x N box (defect_mwr). Raises KeyError if the size has no C_N (no silent default)."""
    C = cfg["alignment"]["C_N_eV"]
    if size not in C:
        raise KeyError(f"production config: no alignment C_N for size '{size}' (sizes: {', '.join(C)})")
    return float(C[size])


def dense_paths(cfg, size, scratch="/home/gregb26/links/scratch/qe_tmp", root=ROOT):
    """Dense primitive .save, dense M file and dense Wannier dir for a supercell size."""
    D = cfg["dense"][size]["D"]
    return dict(D=D, p=cfg["dense"][size]["p"],
                uc=f"{scratch}/defect_uc_dense_{D}/defect_uc_dense_{D}.save",
                mfile=os.path.join(matrices_dir(cfg, root), f"M_dense_{size}.npy"),
                wdir=wannier_dir(cfg, D, root=root),
                manifest=os.path.join(wannier_dir(cfg, D, root=root), "wannier_manifest.json"))
