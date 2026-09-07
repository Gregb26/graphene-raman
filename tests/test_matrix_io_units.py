"""Units gate: M files must be tagged units='hartree'; conversion to eV happens only in load_M_checked."""
import json, os
import numpy as np, pytest
from electron_defect_interaction.io import matrix_io

def _write(tmp_path, meta):
    f = str(tmp_path / "M.npy"); np.save(f, np.ones((2, 3, 2, 3)) * 0.5)
    json.dump(meta, open(str(tmp_path / "M.json"), "w")); return f

def test_tagged_hartree_loads_and_converts(tmp_path):
    f = str(tmp_path / "M.npy"); matrix_io.save_M(f, np.ones((2, 3, 2, 3)) * 0.5, matrix_io.UNIT_CELL)
    assert json.load(open(str(tmp_path / "M.json")))["units"] == "hartree"
    Ha = matrix_io.load_M_checked(f, matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    eV = matrix_io.load_M_checked(f, matrix_io.UNIT_CELL, units=matrix_io.EV)
    assert np.allclose(Ha, 0.5) and np.allclose(eV, 0.5 * matrix_io.HA2EV)

def test_untagged_M_is_refused(tmp_path):
    f = _write(tmp_path, {"bloch_norm": "unit_cell", "shape": [2, 3, 2, 3]})      # legacy sidecar: no units key
    with pytest.raises(ValueError, match="units"):
        matrix_io.load_M_checked(f, matrix_io.UNIT_CELL, units=matrix_io.EV)
    with pytest.raises(ValueError, match="units"):
        matrix_io.load_M_checked(f, matrix_io.UNIT_CELL, units=matrix_io.HARTREE)

def test_wrong_unit_tag_is_refused(tmp_path):
    f = _write(tmp_path, {"bloch_norm": "unit_cell", "units": "eV", "shape": [2, 3, 2, 3]})
    with pytest.raises(ValueError, match="units"):
        matrix_io.load_M_checked(f, matrix_io.UNIT_CELL, units=matrix_io.EV)

def test_caller_must_state_units(tmp_path):
    f = str(tmp_path / "M.npy"); matrix_io.save_M(f, np.zeros((1, 1, 1, 1)), matrix_io.UNIT_CELL)
    with pytest.raises(ValueError, match="explicit"):
        matrix_io.load_M_checked(f, matrix_io.UNIT_CELL)

def test_save_refuses_non_hartree(tmp_path):
    with pytest.raises(ValueError):
        matrix_io.save_M(str(tmp_path / "M.npy"), np.zeros((1, 1, 1, 1)), matrix_io.UNIT_CELL, units="eV")
