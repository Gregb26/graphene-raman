"""
pw_utils.py
    Python module containing helper functions to perform operations of planewaves quantities.
"""

import numpy as np

def mask_invalid_G(nG):
    """
    Returns a boolean mask that selects only the active recripocal lattice vector G for each k-point.

    Inputs:
        nG: (nkpt, ) array of int
            Number of active G for each kpoint
    """

    nG_max = np.max(nG)
    id_G_valid = np.arange(nG_max)[None, :] # (1, nG_max)
    keep = id_G_valid < nG[:, None] # (nkpt, nG_max)

    return keep