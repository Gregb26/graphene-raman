"""R6 étape 3.0 : lecteurs de M -> config.results_dir(cfg) ; porte v2 ; gardes 'non rejouable'. Idempotent (RES déjà défini => rien)."""
import re, os, sys
os.chdir(os.environ["GRAPHENE_RAMAN"])
V2 = "require_normalization=matrix_io.M_NORM_V2"
TAG = 'M_normalization="v2 : L et NL en norme unit_cell, 2026-09-25"'

# ---------------------------------------------------------------- src
p = "src/electron_defect_interaction/io/matrix_io.py"; s = open(p).read()
if "def check_manifest" not in s:
    s = s.replace('''def load_M_checked(npy_path, require_bloch_norm=SUPERCELL, units=None, require_normalization=None):''',
'''def check_manifest(npy_path, require_bloch_norm=UNIT_CELL, require_normalization=None):
    """Validate the sidecar of an M file WITHOUT loading it (for mmap consumers): same refusals as load_M_checked
    (missing sidecar, bloch_norm, units, and M_normalization when require_normalization is given). Returns the manifest."""
    meta = read_manifest(npy_path)
    if meta is None:
        raise ValueError(f"{npy_path}: no manifest sidecar (.json), unknown Bloch normalization. Refusing to run.")
    if meta.get("bloch_norm") != require_bloch_norm:
        raise ValueError(f"{npy_path}: bloch_norm={meta.get('bloch_norm')!r} but '{require_bloch_norm}' required. Refusing to run.")
    if meta.get("units") != HARTREE:
        raise ValueError(f"{npy_path}: manifest units={meta.get('units')!r}, expected '{HARTREE}'. Refusing to run.")
    if require_normalization is not None:
        ver = str(meta.get("M_normalization", M_NORM_V1))
        if not ver.startswith(require_normalization):
            raise ValueError(f"{npy_path}: M_normalization={ver!r} but {require_normalization!r} required. Refusing to run "
                             f"(pre-R6 file: M^L in supercell norm, factor N_cells missing; use results/M2).")
    return meta


def load_M_checked(npy_path, require_bloch_norm=SUPERCELL, units=None, require_normalization=None):''')
    open(p, "w").write(s); print("matrix_io: check_manifest ajouté")

# ---------------------------------------------------------------- scripts python
PROD = ["make_figures.py", "analyze_M.py", "make_figures_memoire.py", "nkint_check_post.py", "compute_convergence.py", "analyze_Ved.py",
        "resonance_criteria.py", "level2_families.py", "ks_reconstruction_all.py", "sampling_table.py", "make_figures_epw.py", "check_onsite_and_NL.py",
        "validate_ML_grid_7x7.py", "test_local_tmatrix_real.py", "resonance_metrics.py", "_old_vs_new_7x7.py", "mwr_locality_coarse_vs_dense.py",
        "m_rcut_convergence.py", "_mcheck.py", "lnl_frobenius_all.py", "epw_ed_vs_ep.py", "check_M_dense_vs_coarse.py", "tag_vacancy_sublattice.py",
        "summarize_level1_maps.py", "rcut_resigma.py", "_normtest.py", "compute_spectral_wannier.py", "check_ML_coarse_kernel.py",
        "check_M_dense_nb20_vs_nb16.py", "_bz_ratio_LNL.py", "compute_M.py", "compute_M_dense_stages.py",
        "compute_spectral.py", "compute_tmatrix.py", "_eta_scan.py", "migrate_M_norm.py"]
OBSOLETE = {"compute_spectral.py", "compute_tmatrix.py", "_eta_scan.py", "migrate_M_norm.py"}
NEEDS_V2 = {"compute_spectral_wannier.py", "mwr_locality_coarse_vs_dense.py", "m_rcut_convergence.py", "resonance_metrics.py", "rcut_resigma.py",
            "resonance_criteria.py", "test_local_tmatrix_real.py", "check_onsite_and_NL.py", "compute_M.py", "compute_M_dense_stages.py"}
IMPORT_RE = re.compile(r"^(import |from )")

def first_import_block_end(lines):
    """index après la dernière ligne du premier bloc d'imports (docstring de tête ignorée)."""
    i = 0; n = len(lines)
    if lines and lines[0].lstrip().startswith(('"""', "'''")):
        q = lines[0].lstrip()[:3]
        if lines[0].count(q) >= 2 and len(lines[0].strip()) > 3:
            i = 1
        else:
            i = 1
            while i < n and q not in lines[i]:
                i += 1
            i += 1
    while i < n and (not lines[i].strip() or lines[i].lstrip().startswith("#")):
        i += 1
    last = None
    while i < n and (IMPORT_RE.match(lines[i]) or (last is not None and lines[i].startswith(" ") and lines[i - 1].rstrip().endswith((",", "(", "\\")))):
        last = i; i += 1
    return (last + 1) if last is not None else 0

changed = []
for name in PROD:
    p = os.path.join("scripts", name); s0 = open(p).read(); s = s0
    uses_res = "results/M/" in s
    # 1. littéraux -> RES
    s = s.replace('f"results/M/', 'f"{RES}/')
    s = s.replace('"results/M/', 'f"{RES}/')
    s = s.replace("f'results/M/", "f'{RES}/").replace("'results/M/", "f'{RES}/")
    s = s.replace("results/M/", "<results_dir>/")          # docstrings, commentaires, messages
    # 2. définition de RES
    if uses_res and "RES = results_dir(" not in s:
        lines = s.split("\n")
        ci = next((i for i, l in enumerate(lines) if l.startswith("from electron_defect_interaction.config import")), None)
        if ci is not None:
            names = lines[ci].split("import", 1)[1]
            add = [x for x in ("load_production", "results_dir") if x not in names]
            if add:
                lines[ci] = lines[ci].rstrip() + ", " + ", ".join(add)
            lines.insert(ci + 1, "RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)")
        else:
            k = first_import_block_end(lines)
            lines.insert(k, "from electron_defect_interaction.config import load_production, results_dir")
            lines.insert(k + 1, "RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)")
        s = "\n".join(lines)
    # 3. porte v2 sur load_M_checked
    if name in NEEDS_V2:
        s = re.sub(r"load_M_checked\(([^()]*?)\)", lambda m: m.group(0) if "require_normalization" in m.group(1) else f"load_M_checked({m.group(1)}, {V2})", s)
    # 4. cas particuliers
    if name in ("analyze_M.py", "lnl_frobenius_all.py"):
        s = s.replace('''    meta = matrix_io.read_manifest(path); assert meta and meta.get("units") == matrix_io.HARTREE, f"{path}: untagged/non-Hartree sidecar"
    return np.load(path, mmap_mode="r")''', '''    matrix_io.check_manifest(path, require_normalization=matrix_io.M_NORM_V2)     # R6 : sidecar v2 exigé (refus sinon)
    return np.load(path, mmap_mode="r")''')
    if name == "analyze_M.py":
        s = s.replace('''            Xb = np.array(mmap_M(pb)) if matrix_io.read_manifest(pb) else np.load(pb)   # legacy coarse-kernel output (Hartree, no sidecar): plumbing test only''',
                      '''            Xb = np.array(mmap_M(pb))''')
        old = '''        A16 = mmap_M(f"{RES}/M_dense_{S}_nb16.npy"); rng = np.random.default_rng(0); ks = rng.choice(len(kd), 30, replace=False); worst = 0.0
        for i in ks:
            for j in ks:
                sa = np.linalg.svd(np.array(A16[:15, i, :15, j]), compute_uv=False); sb = np.linalg.svd(M16[:15, i, :15, j], compute_uv=False); worst = max(worst, np.abs(sa - sb).max() / sa[0])
        print(f"[pad] {S}: nb20[:16] vs nb16, bands 1..15, 30x30 k pairs: {worst:.2e}", flush=True)
        tests.append((f"non-régression nbnd 16 → 20, {S}", "max écart relatif des valeurs singulières, bandes 1–15", f"{worst:.1e}", "1e-5", "OK" if worst < 1e-5 else "ÉCHEC", "analyze_M.py"))'''
        new = '''        p16 = f"{RES}/M_dense_{S}_nb16.npy"
        if _os.path.exists(p16):
            A16 = mmap_M(p16); rng = np.random.default_rng(0); ks = rng.choice(len(kd), 30, replace=False); worst = 0.0
            for i in ks:
                for j in ks:
                    sa = np.linalg.svd(np.array(A16[:15, i, :15, j]), compute_uv=False); sb = np.linalg.svd(M16[:15, i, :15, j], compute_uv=False); worst = max(worst, np.abs(sa - sb).max() / sa[0])
            print(f"[pad] {S}: nb20[:16] vs nb16, bands 1..15, 30x30 k pairs: {worst:.2e}", flush=True)
            tests.append((f"non-régression nbnd 16 → 20, {S}", "max écart relatif des valeurs singulières, bandes 1–15", f"{worst:.1e}", "1e-5", "OK" if worst < 1e-5 else "ÉCHEC", "analyze_M.py"))
        else:
            print(f"[pad] {S}: nb20 vs nb16 non rejouable (fichier absent : {p16}, artefacts nbnd 16 supprimés au ménage de septembre)", flush=True)
            tests.append((f"non-régression nbnd 16 → 20, {S}", "max écart relatif des valeurs singulières, bandes 1–15", "non rejouable", "1e-5", "fichier absent", "analyze_M.py"))'''
        assert old in s, "analyze_M nb16 block introuvable"; s = s.replace(old, new)
    if name == "resonance_metrics.py":
        s = s.replace('ML_diag_mean = float(np.mean([np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r")',
                      'matrix_io.check_manifest(dp["mfile"].replace("M_dense_", "M_L_dense_"), require_normalization=matrix_io.M_NORM_V2)     # R6\nML_diag_mean = float(np.mean([np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r")')
    if name == "_bz_ratio_LNL.py":
        s = s.replace('    ML = np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r"); MN = np.load(dp["mfile"].replace("M_dense_", "M_NL_dense_"), mmap_mode="r")',
                      '    for _f in (dp["mfile"].replace("M_dense_", "M_L_dense_"), dp["mfile"].replace("M_dense_", "M_NL_dense_")): matrix_io.check_manifest(_f, require_normalization=matrix_io.M_NORM_V2)\n    ML = np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r"); MN = np.load(dp["mfile"].replace("M_dense_", "M_NL_dense_"), mmap_mode="r")')
    if name == "check_onsite_and_NL.py":
        s = re.sub(r'^(    )(ML = np\.load\(f"\{RES\}/M_L_dense_\{S\}\.npy", mmap_mode="r"\))', r'\1for _f in ("M_L_dense", "M_NL_dense", "M_dense"): matrix_io.check_manifest(f"{RES}/{_f}_{S}.npy", require_normalization=matrix_io.M_NORM_V2)\n\1\2', s, flags=re.M)
    if name == "check_M_dense_vs_coarse.py":
        s = s.replace('Mc = np.load(f"{RES}/M_L_{size}.npy") if len(sys.argv) < 3 else np.load(sys.argv[2])',
                      'from electron_defect_interaction.io import matrix_io\nfor _f in ([f"{RES}/M_L_{size}.npy"] if len(sys.argv) < 3 else [sys.argv[2]]) + ([f"{RES}/M_L_dense_{size}.npy"] if len(sys.argv) < 4 else [sys.argv[3]]): matrix_io.check_manifest(_f, require_normalization=matrix_io.M_NORM_V2)\nMc = np.load(f"{RES}/M_L_{size}.npy") if len(sys.argv) < 3 else np.load(sys.argv[2])')
    if name == "check_ML_coarse_kernel.py":
        s = s.replace('A = np.load(f"{RES}/M_L_{size}.npy"); B = np.load(f"{RES}/M_L_dense_{size}_coarsecheck.npy")',
                      'from electron_defect_interaction.io import matrix_io\nfor _f in (f"{RES}/M_L_{size}.npy", f"{RES}/M_L_dense_{size}_coarsecheck.npy"): matrix_io.check_manifest(_f, require_normalization=matrix_io.M_NORM_V2)\nA = np.load(f"{RES}/M_L_{size}.npy"); B = np.load(f"{RES}/M_L_dense_{size}_coarsecheck.npy")')
    if name == "check_M_dense_nb20_vs_nb16.py":
        s = s.replace('A=np.load(f"{RES}/M_dense_{S}_nb16.npy", mmap_mode="r"); B=np.load(f"{RES}/M_dense_{S}.npy", mmap_mode="r")',
                      'import os\nif not os.path.exists(f"{RES}/M_dense_{S}_nb16.npy"): raise SystemExit(f"[non rejouable] {RES}/M_dense_{S}_nb16.npy absent (artefacts nbnd 16 supprimés au ménage de septembre 2026)")\nfrom electron_defect_interaction.io import matrix_io; matrix_io.check_manifest(f"{RES}/M_dense_{S}.npy", require_normalization=matrix_io.M_NORM_V2)\nA=np.load(f"{RES}/M_dense_{S}_nb16.npy", mmap_mode="r"); B=np.load(f"{RES}/M_dense_{S}.npy", mmap_mode="r")')
    if name == "_old_vs_new_7x7.py":
        s = s.replace('old = np.load(f"{RES}/obsolete_grid_7x7/M_L_dense_7x7.npy", mmap_mode="r")',
                      'import os\nif not os.path.isdir(f"{RES}/obsolete_grid_7x7"): raise SystemExit(f"[non rejouable] {RES}/obsolete_grid_7x7 absent (ancien M 7x7 sur grille 216 supprimé au ménage de septembre 2026)")\nold = np.load(f"{RES}/obsolete_grid_7x7/M_L_dense_7x7.npy", mmap_mode="r")')
    if name == "validate_ML_grid_7x7.py":
        s = s.replace('    ref = np.load(f"{RES}/_ML_7x7_ref_b4.npy");',
                      '    import os\n    if not os.path.exists(f"{RES}/_ML_7x7_ref_b4.npy"): raise SystemExit(f"[non rejouable] {RES}/_ML_7x7_ref_b4.npy absent : relancer d\'abord `validate_ML_grid_7x7.py ref` (noyau série corrigé, 09280bb)")\n    ref = np.load(f"{RES}/_ML_7x7_ref_b4.npy");')
    if name in ("_mcheck.py", "_normtest.py"):
        s = s.replace('M=np.load(f"{RES}/M_ed_{N}.npy")', 'matrix_io.check_manifest(f"{RES}/M_ed_{N}.npy", require_normalization=matrix_io.M_NORM_V2); M=np.load(f"{RES}/M_ed_{N}.npy")')
        s = s.replace('M = np.load(f"{RES}/M_ed_{N}.npy")', 'matrix_io.check_manifest(f"{RES}/M_ed_{N}.npy", require_normalization=matrix_io.M_NORM_V2); M = np.load(f"{RES}/M_ed_{N}.npy")')
        if "from electron_defect_interaction.io import matrix_io" not in s:
            s = s.replace("from electron_defect_interaction.io import qe_io", "from electron_defect_interaction.io import qe_io, matrix_io", 1)
    if name == "compute_convergence.py":
        s = s.replace('''    for f in sorted(glob.glob(f"{RES}/{prefix}_*x*.npz")):
        m = SIZE_RE.search(f)
        if m:
            out[int(m.group(1))] = np.load(f)
    return out''', '''    for f in sorted(glob.glob(f"{RES}/{prefix}_*x*.npz")):
        m = SIZE_RE.search(f)
        if m:
            out[int(m.group(1))] = np.load(f)
    if not out:
        raise SystemExit(f"[non rejouable] aucun {RES}/{prefix}_*x*.npz (chaîne compute_spectral/compute_tmatrix obsolète depuis 2026-09-05)")
    return out''')
    if name == "compute_M_dense_stages.py":
        s = s.replace('kernel=a.kernel)', f'kernel=a.kernel, {TAG})')
        s = s.replace('part="M_NL_dense", p=P["p"], D=P["D"])', f'part="M_NL_dense", p=P["p"], D=P["D"], {TAG})')
        s = s.replace('matrix_io.save_M(a.out, M, matrix_io.UNIT_CELL, p=P["p"], D=P["D"], N_kd=int(N_kd))', f'matrix_io.save_M(a.out, M, matrix_io.UNIT_CELL, p=P["p"], D=P["D"], N_kd=int(N_kd), N_cells=int(int(a.size.split("x")[0]) ** 2), {TAG})')
    if name in OBSOLETE and "[obsolète, R6" not in s:
        lines = s.split("\n"); k = first_import_block_end(lines)
        lines.insert(k, 'raise SystemExit("[obsolète, R6 2026-09-26] ce script lit M_ed_*_norm.npy / gamma_*.npz (supprimés) et appartient à la chaîne d\'avant le 2026-09-05 ; chaîne courante : compute_spectral_wannier.py, resonance_metrics.py, resonance_criteria.py")')
        s = "\n".join(lines)
    if s != s0:
        open(p, "w").write(s); changed.append(name)
print("python modifiés :", len(changed)); print(" ".join(changed))

# ---------------------------------------------------------------- submit_*.sh
SH = ["submit_M.sh", "submit_M_dense.sh", "submit_spectral_wannier.sh", "submit_spectral_wannier_dense.sh", "submit_rcut_resigma.sh", "submit_nkint_check.sh",
      "submit_golden_dense.sh", "submit_lnl_frobenius.sh", "submit_tmatrix.sh", "submit_spectral.sh", "_test_mnl_mpi.sh"]
RESLINE = '''export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
RES=$("$PROJ/.venv/bin/python" -c 'from electron_defect_interaction.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')   # results/M2 (results/M gelé, R6)
mkdir -p "$RES/logs"'''
chsh = []
for name in SH:
    p = os.path.join("scripts", name); s0 = open(p).read(); s = s0
    s = s.replace("results/M/logs", "results/M2/logs")            # #SBATCH : littéral obligatoire (pas d'expansion), = results_dir de config/production.json
    s = s.replace('"results/M/', '"$RES/').replace("results/M/", "$RES/")
    if 'RES=$(' not in s:
        lines = s.split("\n")
        mi = next((i for i, l in enumerate(lines) if l.startswith("module restore") or l.startswith("module load")), None)
        # après la dernière ligne module ... consécutive
        j = mi
        while j is not None and j + 1 < len(lines) and lines[j + 1].startswith("module "):
            j += 1
        if j is None:
            raise SystemExit(f"{name}: pas de ligne module")
        lines.insert(j + 1, RESLINE); s = "\n".join(lines)
    if name == "submit_M_dense.sh":
        s = s.replace(' --out-norm "$RES/M_dense_${SIZE}_norm.npy"', '')
    if s != s0:
        open(p, "w").write(s); chsh.append(name)
print("shell modifiés :", len(chsh)); print(" ".join(chsh))
