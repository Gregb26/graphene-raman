"""
_paths.py
    Single source of truth for the local data paths used by the standalone scripts. Edit the layout
    here once instead of in every script's main(). Override the data root with the EDI_DATA env var.

    Drivers/tests that already take paths on the command line (compute_M_cluster.py,
    test_pad_vs_full_supercell.py, compare_bands_w90_qe.py) do not need this module.

    Available locally (2026-06): unit cells 5x5 / 11x11 / 12x12; supercell pair 5x5 (p, d) with their
    Vks; Wannier + bands.dat only for the 11x11 unit cell.
"""
import os

ROOT = os.environ.get("EDI_DATA", "data/graphene")


def uc(grid):
    """Unit-cell prefix.save dir for a given k-grid size, e.g. uc('5x5')."""
    return f"{ROOT}/unit_cell/qe/defect_unit_cell_{grid}.save"


def sc_p(grid):
    """Pristine supercell prefix.save dir."""
    return f"{ROOT}/supercell/qe/defect_{grid}_p.save"


def sc_d(grid):
    """Defective supercell prefix.save dir."""
    return f"{ROOT}/supercell/qe/defect_{grid}_d.save"


def pot_p(grid):
    """Pristine supercell local potential (pp.x plot_num=1 filplot), inside the .save dir."""
    return f"{sc_p(grid)}/Vks_{grid}_p"


def pot_d(grid):
    """Defective supercell local potential."""
    return f"{sc_d(grid)}/Vks_{grid}_d"


def uc_pot(grid):
    """Unit-cell local potential (pp.x), inside the unit-cell .save dir."""
    return f"{uc(grid)}/Vks_uc_{grid}"


def upf(grid="5x5"):
    """UPF pseudopotential shipped in the unit-cell .save dir."""
    return f"{uc(grid)}/C.upf"
