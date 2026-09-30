"""Chain T validation scripts (scripts/validation/) run under pytest: golden tests of the local t-matrix
against the dense compute_T reference, and local_green_batch against local_green. Each script prints its
own PASS/FAIL and returns 0/1; the assertion is on the return code."""
import pytest
from conftest import local_data
from graphene_raman.config import load_production, matrices_dir


def test_local_tmatrix_golden(repo_cwd):
    """Synthetic golden test: local t (all R) = dense compute_T x N_k, and positivity."""
    import test_local_tmatrix
    assert test_local_tmatrix.main() == 0


def test_local_rcut(repo_cwd):
    """R_cut truncation, extract_V_loc and mwr_locality on the synthetic model."""
    import test_local_rcut
    assert test_local_rcut.main() == 0


@pytest.mark.slow
def test_local_green_batch(repo_cwd):
    """local_green_batch = local_green to 1e-12 (wannier/27x27); scattering_rate_fast vs exact (C8)."""
    import test_local_green_batch
    assert test_local_green_batch.main() == 0


@pytest.mark.cluster
@pytest.mark.slow
@pytest.mark.skipif(not local_data(f"{matrices_dir(load_production(verbose=False), root='')}/M_ed_5x5.npy".lstrip('/'),
                                   "data/graphene/unit_cell/qe/defect_5x5.save"),
                    reason="real golden test: needs the 5x5 M2 matrix and the cluster naming of the .save")
def test_local_tmatrix_real(repo_cwd, monkeypatch):
    """Real golden test (C6) on the coarse 5x5 M2: local t = dense compute_T x N_cells."""
    import runpy
    monkeypatch.setattr("sys.argv", ["test_local_tmatrix_real.py", "5x5"])
    with pytest.raises(SystemExit) as e:
        runpy.run_path(str(repo_cwd / "scripts/validation/test_local_tmatrix_real.py"), run_name="__main__")
    assert e.value.code == 0
