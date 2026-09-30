"""
Tight-binding model in the Wannier basis (`WannierTB`) and the analytic graphene model of M0.

M1: `make_wannier_tb` builds it from Wannier90's `_tb.dat` (read in io/wannier_io.py). Transformations
of a model: `centres_only` (M1), `pz_block` (M3); `extract_block` cuts a sub-block of the orbital axes.
Conventions: package docstring (electron_photon/__init__.py).
"""

import numpy as np
import dataclasses
from graphene_raman.io.wannier_io import read_w90_tb

@dataclasses.dataclass
class WannierTB:
    """
    Tight-binding model in the Wannier basis (layout of Wannier90's `_tb.dat`).

    Fields:
        lattice : (3, 3) float, Angstrom, lattice[:, i] = a_i
        R_int   : (nR, 3) int, R vectors in reduced coordinates, closed under R -> -R
        R_cart  : (nR, 3) float, Angstrom, R_int @ lattice.T
        ndegen  : (nR,) int, Wigner-Seitz degeneracies (1 for the analytic model)
        H_R     : (nR, nW, nW) complex, eV, H_mn(R) = <0m| H |Rn>
        r_R     : (nR, 3, nW, nW) complex, Angstrom, <0m| r_alpha |Rn>; diagonal of r(0) = centres
        index   : dict, tuple(R) -> row of R in R_int
        minus   : (nR,) int, row of -R, so that X_R[minus] = X(-R)
        t, a_cc : float, eV and Angstrom, hopping and bond length of the toy model (defaults for
                  Wannier90 data; only `ring` uses them, for q0)
    """
    lattice: np.ndarray
    R_int: np.ndarray
    R_cart: np.ndarray
    ndegen: np.ndarray
    H_R: np.ndarray
    r_R: np.ndarray
    index: dict
    minus: np.ndarray
    t: float
    a_cc: float

def make_graphene_tb(t=2.7, a_cc=1.42, c=15.0, shift_B=(0,0,0)):
    """
    Nearest-neighbour tight-binding model of graphene, one p_z orbital per atom (A = 0, B = 1), zero
    on-site energies (Dirac point at 0 eV) and H_AB(R) = -t on the three bonds. shift_B = L moves the B
    centre of cell 0 by a lattice vector, the gauge choice Wannier90 makes: the bonds are unchanged, the
    hop labels move to R - L, and only the Berry term keeps hbar v_cv invariant (EM.md, M0 notes).

    Inputs:
        t       : float, eV, hopping (matrix element -t)
        a_cc    : float, Angstrom, C-C bond length (lattice constant sqrt(3) a_cc)
        c       : float, Angstrom, length of a3 (keeps the vectors 3D)
        shift_B : 3 ints, reduced lattice vector L; (0,0,0) gives nR = 5, (1,0,0) gives nR = 7
    Returns:
        WannierTB with ndegen = 1 and r(0) = diag(tau_A, tau_B + L): centres only
    """

    # lattice vectors (Angstrom): a1 = a (1, 0, 0), a2 = a (1/2, sqrt(3)/2, 0), a3 = (0, 0, c)
    a = np.sqrt(3)*a_cc # graphene lattice constant
    a1 = a*np.array([1,0,0]); a2 = a*np.array([1,np.sqrt(3),0])/2; a3=np.array([0,0,c])
    lattice = np.column_stack((a1, a2, a3)) # lattice[:,i] = a_i

    # bond vectors from A to its three B neighbours, |delta_i| = a_cc, at 30, 150 and -90 degrees
    delta1 = a_cc*np.array([np.sqrt(3)/2, 0.5, 0]); delta2 = a_cc*np.array([-np.sqrt(3)/2, 0.5, 0]); delta3 = a_cc*np.array([0, -1.0, 0])
    delta = np.column_stack((delta1, delta2, delta3)) # delta[:,i] = delta_i

    # Wannier centres: B of cell 0 at the end of the first bond, moved by the lattice vector shift_B
    tau_A = np.zeros(3)
    tau_B = tau_A + delta[:,0] + lattice @ shift_B

    # one hop per bond: R + tau_B = tau_A + delta, solved in reduced coordinates (must be integers)
    hops = []
    for d in range(3):
        R_red = np.linalg.solve(lattice, tau_A + delta[:,d] - tau_B)
        R = np.rint(R_red).astype(int)
        assert np.allclose(R_red, R, atol=1e-8), 'bond does not end on a B site'
        hops.append(R)

    # R list: A -> B hops, their Hermitian partners at -R, and R = 0 for the centres
    # (np.unique removes the duplicated rows and sorts them; the order is irrelevant)
    R_int = np.unique(np.array(hops + [-R for R in hops]+[np.zeros(3, int)]), axis=0)
    R_cart = R_int @ lattice.T
    index = {tuple(R): iR for iR, R in enumerate(R_int)} # integer triplet R -> row of R in R_int
    minus = np.array([index[tuple(-R)] for R in R_int]) # KeyError if some -R is missing

    nR = len(R_int)
    nW = 2 # two p_z orbitals in a unit cell

    H_R = np.zeros((nR, nW, nW), complex)
    r_R= np.zeros((nR, 3, nW, nW), complex)
    ndegen = np.ones(nR, dtype=int)

    A,B = 0,1
    for R in hops:
        H_R[index[tuple(R)], A, B] = -t # <0A|H|RB> : A of cell 0 hops to B of cell R
        H_R[index[tuple(-R)], B, A] = - t # <0B|H|-RA> : Hermitian partner, H_BA(-R) = H_AB(R)^*

    # position operator: Wannier centres on the diagonal of r(0) (the B centre includes the shift)
    i0 = index[(0,0,0)]
    r_R[i0, :, A, A] = tau_A
    r_R[i0,:, B, B] = tau_B

    return WannierTB(lattice=lattice, R_int=R_int, R_cart=R_cart, ndegen=ndegen, H_R=H_R, r_R=r_R, t=t, a_cc=a_cc, index=index, minus=minus)

def make_wannier_tb(path):
    """
    WannierTB from a Wannier90 `_tb.dat` (read by io/wannier_io.read_w90_tb). H_R and r_R stay raw
    (no ndegen division, r not hermitized). t and a_cc take the M0 values: only `ring` uses them (q0).

    Inputs:
        path : str or Path, `seedname_tb.dat`
    Returns:
        WannierTB (Angstrom, eV)
    """
    H_R, R_int, ndegen, r_R, lattice = read_w90_tb(path)
    R_cart = R_int @ lattice.T
    index = {tuple(R): iR for iR, R in enumerate(R_int)} # integer triplet R -> row of R in R_int
    minus = np.array([index[tuple(-R)] for R in R_int]) # KeyError if some -R is missing
    
    return WannierTB(lattice=lattice, R_int=R_int, R_cart=R_cart, ndegen=ndegen, H_R=H_R, r_R=r_R, index=index, minus=minus, t=2.7, a_cc=1.42)

def centres_only(tb):
    """
    Copy of `tb` whose r keeps only the Wannier centres (diagonal of r(0)): A(k) = diag(tau_n), the
    tight-binding (Peierls) approximation. Use as compute_velocity(centres_only(tb), k, 'berry').
    `tb` is not modified; the other fields are shared with it (dataclasses.replace).
    """
    r_R_centres = np.zeros_like(tb.r_R) # (nR, 3, nW, nW)
    iR0 = tb.index[(0,0,0)]
    r_0 = tb.r_R[iR0] # (3, nW, nW)
    r_0_diag = r_0 * np.eye(r_0.shape[-1]) # (3, nW, nW), diagonal
    r_R_centres[iR0] = r_0_diag
    return dataclasses.replace(tb, r_R=r_R_centres)

def extract_block(X, rows, columns):
    """Sub-block rows x columns of the last two axes of X (X[..., rows, columns] would pair them)."""
    return X[..., rows, :][..., columns]

def pz_block(tb, pz):
    """
    Copy of `tb` reduced to the Wannier functions `pz`: H_R and r_R cut to their pz x pz block (their
    centres included); R list and ndegen unchanged. Exact for in-plane light when the mirror decouples
    sigma and p_z (symmetry_report): only r^z between them is dropped.

    Inputs:
        tb : WannierTB
        pz : list of int, Wannier indices to keep (a list, so the blocks are copies)
    Returns:
        WannierTB with nW = len(pz); `tb` is not modified
    """
    H_R_pz = extract_block(tb.H_R, pz, pz)
    r_R_pz = extract_block(tb.r_R, pz, pz)

    return dataclasses.replace(tb, H_R=H_R_pz, r_R=r_R_pz)
