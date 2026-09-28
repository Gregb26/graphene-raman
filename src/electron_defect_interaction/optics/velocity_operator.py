"""
Wannier-interpolated velocity operator (M0: tight-binding toy model and Fourier machinery).

Goal of the module: build, from a Wannier tight-binding description (H(R), r(R)) as written by
Wannier90 in `_tb.dat`, the interband velocity matrix

    hbar v(k) = V(k)^dagger [ dH(k)/dk + i [H(k), A(k)] ] V(k),

where V(k) diagonalizes H(k) and A(k) is the Berry connection in the Wannier gauge. In the M0
phase the tight-binding data come from an analytic graphene model (`make_graphene_tb`) that has
exactly the same layout as what the M1 reader of the real `_tb.dat` will return, so every
downstream function is first validated on a model with known answers.

Conventions (fixed once for all, shared with the M1 reader):
    - Units: lengths in Angstrom, energies in eV, wavevectors in 1/Angstrom, hbar*v in eV*Angstrom.
    - Vectors are always 3D (z component = 0 for the 2D sheet), as in `_tb.dat`.
    - Matrix elements in the Wannier basis: X_mn(R) = <0m| X |Rn>, i.e. the row m is the Wannier
      function of the home cell 0 and the column n the one of cell R.
    - Fourier transform in the lattice gauge (no intracell positions tau in the phases):
          X(k) = sum_R exp(+i k.R) X(R) / ndegen(R)
      A single routine (`fourier`) does it for H, for i R_mu H (the k-derivative) and for r (the
      Berry connection A), so that sign and ndegen conventions cannot diverge between them.
    - Array layout: `lattice` and the reciprocal matrix `B` store their basis vectors in COLUMNS
      (lattice[:, i] = a_i, B[:, j] = b_j). Lists of vectors (R, k, bond vectors) are stored in
      ROWS, shape (N, 3). Hence a single vector converts as x_cart = lattice @ n_red, while a list
      of vectors converts as X_cart = N_red @ lattice.T.
"""

import numpy as np
from dataclasses import dataclass

@dataclass
class WannierTB:
    """
    Tight-binding model in the Wannier basis, in the layout of a Wannier90 `_tb.dat` file.

    Fields:
        lattice : (3, 3) float, Angstrom
            Primitive lattice vectors stored in columns, lattice[:, i] = a_i.
        R_int : (nR, 3) int
            Lattice vectors R in reduced (integer) coordinates, one per row; R_cart = R_int @ lattice.T.
            The list must be closed under R -> -R, otherwise H(k) cannot be Hermitian.
        ndegen : (nR,) int
            Wigner-Seitz degeneracy of each R. Every Fourier sum divides by it:
            X(k) = sum_R e^{ik.R} X(R) / ndegen(R). Equal to 1 for the analytic model.
        H_R : (nR, nW, nW) complex, eV
            Hamiltonian matrix elements H_mn(R) = <0m| H |Rn>.
        r_R : (nR, 3, nW, nW) complex, Angstrom
            Position matrix elements r_mn,alpha(R) = <0m| r_alpha |Rn>; axis 1 is the Cartesian
            component alpha. Its diagonal at R = 0 holds the Wannier centres.
        t : float, eV
            Nearest-neighbour hopping of the toy model. Not a Wannier90 quantity: the M1 reader
            will need a default value for it.
        a_cc : float, Angstrom
            Carbon-carbon bond length of the toy model (same remark as for t).
    """
    lattice: np.ndarray
    R_int: np.ndarray
    ndegen: np.ndarray
    H_R: np.ndarray
    r_R: np.ndarray
    t: float
    a_cc: float

def make_graphene_tb(t=2.7, a_cc=1.42, c=15.0, shift_B=(0,0,0)):
    """
    Nearest-neighbour tight-binding model of graphene, one p_z orbital per atom (nW = 2).

    Orbital 0 sits on sublattice A, orbital 1 on sublattice B. The model is built geometrically:
    each A atom has three B neighbours at the bond vectors delta_j (|delta_j| = a_cc), and for each
    bond we ask in which cell R the B orbital at the end of the arrow is labelled. With Wannier
    centres tau_A and tau_B, the B orbital of cell R sits at R + tau_B, so the hop A(0) -> B(R)
    along delta requires

        R + tau_B = tau_A + delta   =>   R = tau_A + delta - tau_B,

    which must be a lattice vector (checked). That hop gets H_AB(R) = -t and its Hermitian partner
    H_BA(-R) = -t. The on-site energies are zero, so the Dirac point sits at 0 eV.

    shift_B is a pure relabelling (a gauge choice), not a change of physics. It moves the centre of
    "the B orbital of cell 0" by the lattice vector L = shift_B, i.e. it declares that orbital to be
    the B atom physically located at tau_B + L. Atoms and bonds are unchanged; only the integer
    labels R of the hops move: with S = {(0,0,0), (-1,0,0), (0,-1,0)} the cells of the B neighbours
    of A when L = 0, H_AB lives on S - L = {s - L} and H_BA on L - S = {L - s}. Wannier90 makes this
    kind of choice by itself when it places a Wannier centre in a neighbouring cell. In the lattice
    gauge, H(k) then picks up a k-dependent phase on the B row and column, the eigenvalues do not
    change but dH/dk alone does; the Berry term i[H, A] compensates it exactly, which is what the
    gauge-shift tests check.

    Inputs:
        t : float, eV
            Nearest-neighbour hopping (positive; the matrix element is -t).
        a_cc : float, Angstrom
            Carbon-carbon bond length. The lattice constant is a = sqrt(3) a_cc.
        c : float, Angstrom
            Length of the out-of-plane lattice vector a3. It plays no role in the physics of the
            sheet but keeps everything 3D, as in `_tb.dat`.
        shift_B : 3-tuple of ints
            Lattice vector, in reduced coordinates, by which the centre of the B orbital is
            relabelled. (0,0,0) gives nR = 5; (1,0,0) gives nR = 7 with identical bonds.
    Returns:
        WannierTB
            ndegen = 1 for every R: the model is defined exactly on its R list, there is no
            Wigner-Seitz truncation of a larger box. r(R) is non-zero only at R = 0, where
            r(0) = diag(tau_A, tau_B + L): the model has no off-diagonal position elements, so in
            M0 the full Berry connection and the "centres only" one coincide.
    """

    # build lattice vectors in cartesian coords (Angstrom):
    # a1 = a (1, 0, 0), a2 = a (1/2, sqrt(3)/2, 0), a3 = (0, 0, c)
    a = np.sqrt(3)*a_cc # graphene lattice constant
    a1 = a*np.array([1,0,0]); a2 = a*np.array([1,np.sqrt(3),0])/2; a3=np.array([0,0,c])
    lattice = np.column_stack((a1, a2, a3)) # lattice[:,i] = a_i

    # cartesian vectors that connect A to its three B neighbours, |delta_i| = a_cc.
    # They point at 30, 150 and -90 degrees from the x axis.
    delta1 = a_cc*np.array([np.sqrt(3)/2, 0.5, 0]); delta2 = a_cc*np.array([-np.sqrt(3)/2, 0.5, 0]); delta3 = a_cc*np.array([0, -1.0, 0])
    delta = np.column_stack((delta1, delta2, delta3)) # delta[:,i] = delta_i

    # Wannier centres (atom positions) in cartesian coords. The B orbital of cell 0 is placed at the
    # end of the first bond, then moved by the lattice vector shift_B (reduced -> cartesian with
    # lattice @ shift_B, since the a_i are the columns of lattice).
    tau_A = np.zeros(3)
    tau_B = tau_A + delta[:,0] + lattice @ shift_B

    # Hopping : one hop per bond. Solve R + tau_B = tau_A + delta for R.
    # np.linalg.solve(lattice, x) returns the coefficients n such that lattice @ n = x, i.e. the
    # reduced coordinates of x. They must be integers, otherwise the bond does not end on a B site
    # (wrong tau or delta).
    hops = []
    for d in range(3):
        R_red = np.linalg.solve(lattice, tau_A + delta[:,d] - tau_B)
        R = np.rint(R_red).astype(int)
        assert np.allclose(R_red, R, atol=1e-8), 'bond does not end on a B site'
        hops.append(R)

    # R list : A-B hops, their Hermitian partners at B-A at -R and R=0 for centers.
    # np.unique(axis=0) removes duplicated rows (R = 0 appears several times when shift_B = 0) and
    # sorts them; the order of the R list is irrelevant for the Fourier sums.
    R_int = np.unique(np.array(hops + [-R for R in hops]+[np.zeros(3, int)]), axis=0)
    index = {tuple(R): iR for iR, R in enumerate(R_int)} # integer triplet R -> row of R in R_int
    nR = len(R_int)
    nW = 2 # two p_z orbitals in a unit cell

    H_R = np.zeros((nR, nW, nW), complex)
    r_R= np.zeros((nR, 3, nW, nW), complex)
    ndegen = np.ones(nR, dtype=int)

    A,B = 0,1
    for R in hops:
        H_R[index[tuple(R)], A, B] = -t # <0A|H|RB> : A of cell 0 hops to B of cell R
        H_R[index[tuple(-R)], B, A] = - t # <0B|H|-RA> : Hermitian partner, H_BA(-R) = H_AB(R)^*

    # Position operator: only the diagonal of R = 0 is non-zero and holds the Wannier centres
    # (the B centre includes the shift, which is what makes the Berry term compensate the gauge).
    i0 = index[(0,0,0)]
    r_R[i0, :, A, A] = tau_A
    r_R[i0,:, B, B] = tau_B

    return WannierTB(lattice=lattice, R_int=R_int, ndegen=ndegen, H_R=H_R, r_R=r_R, t=t, a_cc=a_cc)

def reciprocal(lattice):
    """
    Primitive reciprocal lattice vectors, defined by a_i . b_j = 2 pi delta_ij.

    With the a_i in the columns of `lattice` (call it A), the matrix of dot products a_i . b_j is
    A^T B, so A^T B = 2 pi I gives B = 2 pi (A^T)^{-1}, with the b_j in the columns of B.
    np.linalg.inv (not pinv) is used on purpose: a singular lattice must raise an error instead of
    returning a meaningless pseudo-inverse.

    Inputs:
        lattice : (3, 3) float, Angstrom, lattice[:, i] = a_i
    Returns:
        B : (3, 3) float, 1/Angstrom, B[:, j] = b_j
    """

    B = 2*np.pi*np.linalg.inv(lattice.T) # B[:,i] = b_i

    # consistency check written in the completeness form sum_i a_i b_i^T = A B^T = 2 pi I, which is
    # equivalent to A^T B = 2 pi I for square matrices (a left inverse is also a right inverse)
    assert np.allclose(lattice @ B.T, 2*np.pi*np.eye(3)), 'not lattice vectors'

    return B

def k_grid(B, N, shift=0.5):
    """
    Uniform N x N k-point grid covering one reciprocal unit cell of the 2D lattice.

        k = ((i + s)/N) b1 + ((j + s)/N) b2,    i, j = 0 ... N-1,    s = shift

    The reduced coordinates therefore lie in (0, 1). Any complete reciprocal cell is equivalent for
    Brillouin-zone sums: in the lattice gauge H(k + G) = H(k) exactly, because e^{iG.R} = 1 for
    every lattice vector R, and the same holds for dH and A.

    The default half-step shift s = 0.5 keeps K and K' off the grid when N is a multiple of 3
    (e.g. N = 1800). At K the two bands are degenerate at the chemical potential mu = 0, so the
    occupation eps < mu would be ambiguous there and eigh would return an arbitrary eigenbasis.

    Inputs:
        B : (3, 3) float, 1/Angstrom
            Reciprocal vectors in columns, B[:, j] = b_j (from `reciprocal`).
        N : int
            Number of points along b1 and along b2 (N^2 points in total).
        shift : float
            Offset of the grid, in units of the grid step.
    Returns:
        k_red : (N^2, 3) float
            Reduced coordinates, third component 0. Row order is 'ij' (C order): j (along b2)
            varies fastest, then i (along b1).
        k_cart : (N^2, 3) float, 1/Angstrom
            Cartesian k points in rows, k_cart = k_red @ B.T (row version of k = B @ kappa).
        K : (3,) float, 1/Angstrom
            Dirac point K = (2 b1 + b2)/3, where H_AB(K) = 0.
    """

    u = (np.arange(N) + shift)/N # reduced coordinates along one axis, in (0, 1)

    # all (u1, u2) pairs: two (N, N) arrays, flattened to N^2 aligned entries
    u1, u2 = np.meshgrid(u,u, indexing='ij')
    u3 = np.zeros(N**2)

    k_red = np.column_stack((u1.ravel(), u2.ravel(), u3))
    k_cart = k_red @ B.T # rows of reduced coords -> rows of cartesian vectors

    K = (2*B[:,0] + B[:,1])/3 # Dirac point

    return k_red, k_cart, K

def dagger(X):
    """
    Hermitian conjugate over the last two axes: (X^dagger)_mn = conj(X_nm).

    Only the two orbital (or band) axes are swapped; leading axes (k points, Cartesian component)
    are left untouched, so it works on (Nk, nW, nW) as well as on (Nk, 3, nW, nW) arrays.
    """

    return np.conj(np.swapaxes(X, -2, -1))

def hermitize(X):
    """
    Hermitian part of X over the last two axes: (X + X^dagger)/2.

    Meant for the Berry connection A(k): the position matrix elements r(R) written by Wannier90
    satisfy r_mn(R) = r_nm(-R)^* only approximately, so A(k) is not exactly Hermitian. It is applied
    by the caller after `fourier`, not inside it, so that an unexpected non-Hermiticity shows up in
    the tests instead of being silently removed. It leaves an already Hermitian input (the M0
    model) unchanged and is idempotent.
    """

    return (X + dagger(X))/2


def fourier(X_R, R_cart, ndegen, k, deriv=None):
    """
    Lattice-gauge Fourier transform from the Wannier R basis to arbitrary k points.

        X(k)        = sum_R exp(+i k.R) X(R) / ndegen(R)              (deriv falsy)
        dX(k)/dk_mu = sum_R i R_mu exp(+i k.R) X(R) / ndegen(R)       (deriv truthy)

    This is the only Fourier routine of the module; it is used for three objects:
        H(k)  : X_R = H_R (nR, nW, nW)     -> (Nk, nW, nW),     eV
        dH(k) : X_R = H_R, deriv=True       -> (Nk, 3, nW, nW),  eV*Angstrom (units of hbar*v)
        A(k)  : X_R = r_R (nR, 3, nW, nW)  -> (Nk, 3, nW, nW),  Angstrom
    The phase uses Cartesian k and R. It is the same phase as Hwr_to_Hwk in reduced coordinates,
    since k_cart . R_cart = 2 pi k_red . R_int.

    Implementation: the trailing axes of X_R are flattened to (nR, M) so that the sum over R is a
    single matrix product (Nk, nR) @ (nR, M), then reshaped back to (Nk, ...). The derivative
    multiplies X(R) by i R_mu BEFORE the sum, which inserts the Cartesian axis mu right after the
    k axis.

    Memory: the phase matrix is (Nk, nR) complex128, i.e. 16 Nk nR bytes (about 38 GB for 1800^2
    k points and 741 R vectors). Large grids must be processed by the caller in blocks of 1e4 to
    1e5 k points.

    Inputs:
        X_R : (nR, ...) complex
            Real-space matrix elements. Only a (nR, nW, nW) array can be differentiated.
        R_cart : (nR, 3) float, Angstrom
            Cartesian R vectors in rows, R_cart = R_int @ lattice.T.
        ndegen : (nR,) int
            Wigner-Seitz degeneracies.
        k : (Nk, 3) float, 1/Angstrom
            Cartesian k points in rows. A single k point must be passed as a (1, 3) array,
            e.g. K[None, :].
        deriv : bool
            If truthy, return the Cartesian gradient dX/dk instead of X(k).
    Returns:
        X_k : (Nk, ...) complex, or (Nk, 3, nW, nW) for the derivative
    """

    # shapes
    nR = X_R.shape[0]; Nk = k.shape[0]

    # compute phase scaled by ndegen: phase[k, R] = e^{ik.R} / ndegen(R)
    phase = np.exp(1j * k @ R_cart.T) / ndegen[None, :] # (nK, nR)

    # Fourier transform
    if deriv:
        # a 4D r_R would broadcast silently into a wrong shape, hence the check
        assert X_R.ndim == 3 # differentiate hamiltonian only and not position operator
        # (nR, 1, nW, nW) * (nR, 3, 1, 1) -> (nR, 3, nW, nW): element [R, mu] = i R_mu X(R)
        X_R = X_R[:, None, :, :] * 1j*R_cart[:,:, None, None]

    # trailing shape, read AFTER the derivative (which adds the Cartesian axis)
    tail = X_R.shape[1:]

    # flatten the trailing axes, sum over R as one matrix product, restore the trailing axes.
    # Flattening and restoring both use C order, so no axes are mixed.
    X_R = X_R.reshape(nR, -1) # (nR, M)
    X_k = phase @ X_R # (nK, M)
    X_k = X_k.reshape(Nk, *tail)

    return X_k
