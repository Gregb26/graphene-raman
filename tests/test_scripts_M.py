"""Chain M validation scripts (scripts/validation/) run under pytest: zero-padding of M^L, Kohn-Sham
reconstruction on the unit cell, null-defect check, Wannier interpolation pipeline. The tests that read
the local DFT data (data/graphene, see scripts/_paths.py) are marked needs_data and skipped when absent."""
import pytest
from conftest import local_data, _paths


def test_zero_pad_exact():
    """zero_pad_potential is exact on an analytic potential (no data)."""
    import test_zero_pad_dense
    assert bool(test_zero_pad_dense.test_zero_pad(p=3))


@pytest.mark.needs_data
@pytest.mark.slow
@pytest.mark.skipif(not local_data(_paths.uc("5x5"), _paths.sc_p("5x5"), _paths.sc_d("5x5")), reason="5x5 unit cell and supercell pair absent")
def test_zero_pad_non_regression(repo_cwd):
    """Dense M^L by zero padding vs compute_ML_G on the 5x5 pair (several minutes)."""
    import test_zero_pad_dense
    assert bool(test_zero_pad_dense.test_non_regression(p=3))


@pytest.mark.needs_data
@pytest.mark.slow
@pytest.mark.skipif(not local_data(_paths.uc("5x5"), _paths.sc_p("5x5"), _paths.uc_pot("5x5")), reason="5x5 unit cell, supercell and pp.x potential absent")
def test_ks_reconstruction(repo_cwd):
    """Test A (H = T + V_loc + V^NL = diag eps), B (null defect => M = 0), C (restricted V_p) on the 5x5 data."""
    import test_ks_reconstruction
    assert test_ks_reconstruction.main() == 0


@pytest.mark.needs_data
@pytest.mark.skipif(not local_data(_paths.uc("11x11") + "/wannier_u.mat", _paths.uc("11x11") + "/wannier_tb.dat"), reason="11x11 unit cell with Wannier files absent")
def test_wannier_pipeline(repo_cwd):
    """Parsers, Wannier-gauge round trip, gauge-invariant spectrum, fine-grid smoke test (11x11 local data)."""
    import test_wannier
    assert test_wannier.main() == 0


@pytest.mark.cluster
@pytest.mark.slow
@pytest.mark.skipif(not local_data(_paths.sc_p("6x6"), _paths.sc_p("12x12"), _paths.uc("12x12")), reason="6x6 and 12x12 supercell pairs absent (cluster only)")
def test_pad_vs_full_supercell(repo_cwd, monkeypatch):
    """M^L zero-padded x p on the 6x6 vs the 12x12 supercell (cluster data)."""
    import test_pad_vs_full_supercell
    monkeypatch.setattr("sys.argv", ["test_pad_vs_full_supercell.py",
        "--uc-dense", _paths.uc("12x12"), "--sc-small-p", _paths.sc_p("6x6"), "--pot-small-p", _paths.pot_p("6x6"), "--pot-small-d", _paths.pot_d("6x6"),
        "--sc-large-p", _paths.sc_p("12x12"), "--pot-large-p", _paths.pot_p("12x12"), "--pot-large-d", _paths.pot_d("12x12")])
    assert test_pad_vs_full_supercell.main() == 0
