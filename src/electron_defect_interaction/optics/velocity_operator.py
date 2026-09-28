"""
Wannier-interpolated velocity operator for the electron-photon coupling (EM series, M0).

From a Wannier tight-binding model (H(R), r(R), layout of Wannier90's `_tb.dat`), computes

    hbar v(k) = V(k)^dagger [ dH(k)/dk + i [H(k), A(k)] ] V(k),

with V(k) diagonalizing H(k) and A(k) the Berry connection in the Wannier gauge. Chain:
fourier (H, dH, A) -> hermitize (A) -> velocity, block by block on the k points; `ring` gives the
resonant k points around K. In M0 the data come from the analytic graphene model
`make_graphene_tb`, with the same layout as the future M1 reader. Plan, derivations and reference
values: memoire/EM/EM.md.

Conventions:
    - Units: Angstrom, eV, 1/Angstrom; hbar v in eV*Angstrom.
    - Vectors are always 3D (z = 0), as in `_tb.dat`.
    - X_mn(R) = <0m| X |Rn> and X(k) = sum_R e^{+ik.R} X(R) / ndegen(R) (lattice gauge, no tau in
      the phases), computed by the single routine `fourier` for H, dH/dk and A.
    - `lattice` and `B` hold their basis vectors in COLUMNS; lists of vectors (R, k) are ROWS,
      shape (N, 3), so that X_cart = X_red @ lattice.T.
    - Occupations are decided by energy (eps < mu), never by band index.
"""

import numpy as np
from dataclasses import dataclass

@dataclass
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
        t, a_cc : float, eV and Angstrom, hopping and bond length of the toy model (not Wannier90
                  quantities: the M1 reader will need defaults)
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

@dataclass
class GridTB:
    """
    k-grid quantities for the callers (tests, driver), built by `make_grid_tb`. The computing
    functions take plain k arrays instead, since they also run on blocks, rings and paths.

    Fields:
        B      : (3, 3) float, 1/Angstrom, B[:, j] = b_j
        k_red  : (N^2, 3) float, reduced k points
        k_cart : (N^2, 3) float, 1/Angstrom
        K      : (3,) float, 1/Angstrom, Dirac point (2 b1 + b2)/3
    """
    B: np.ndarray
    k_red: np.ndarray
    k_cart: np.ndarray
    K: np.ndarray


def make_graphene_tb(t=2.7, a_cc=1.42, c=15.0, shift_B=(0,0,0)):
    """
    Nearest-neighbour tight-binding model of graphene, one p_z orbital per atom (A = 0, B = 1).

    Built geometrically: the hop from A(0) along the bond delta lands on the B orbital of cell
    R = tau_A + delta - tau_B (checked to be a lattice vector), with H_AB(R) = H_BA(-R) = -t and
    zero on-site energies (Dirac point at 0 eV). shift_B = L relabels which B atom belongs to cell 0,
    a gauge choice like the one Wannier90 makes when it puts a centre in a neighbouring cell: the
    bonds are unchanged, the hop labels move to R - L, and only the Berry term keeps hbar v_cv
    invariant (EM.md, implementation notes).

    Inputs:
        t       : float, eV, hopping (matrix element -t)
        a_cc    : float, Angstrom, C-C bond length (lattice constant sqrt(3) a_cc)
        c       : float, Angstrom, length of a3 (keeps the vectors 3D)
        shift_B : 3 ints, reduced lattice vector L; (0,0,0) gives nR = 5, (1,0,0) gives nR = 7
    Returns:
        WannierTB with ndegen = 1 and r(0) = diag(tau_A, tau_B + L): centres only, so the full and
        "centres only" Berry connections coincide in M0.
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

def reciprocal(lattice):
    """
    Reciprocal vectors from a_i . b_j = 2 pi delta_ij, i.e. B = 2 pi (lattice^T)^{-1}.

    Inputs:
        lattice : (3, 3) float, Angstrom, lattice[:, i] = a_i
    Returns:
        B : (3, 3) float, 1/Angstrom, B[:, j] = b_j
    """

    B = 2*np.pi*np.linalg.inv(lattice.T) # B[:,i] = b_i; inv, not pinv: a singular lattice must fail

    return B

def k_grid(B, N, shift=0.5):
    """
    Shifted N x N grid over one reciprocal cell, k = ((i + s)/N) b1 + ((j + s)/N) b2.

    Any complete cell is equivalent for Brillouin-zone sums (H(k + G) = H(k) in the lattice gauge).
    The half-step shift keeps K and K' off the grid when N is a multiple of 3: the bands are
    degenerate there and the occupation would be ambiguous.

    Inputs:
        B     : (3, 3) float, 1/Angstrom, B[:, j] = b_j
        N     : int, number of points along b1 and along b2
        shift : float, offset in units of the grid step
    Returns:
        k_red  : (N^2, 3) float, reduced coordinates, j (along b2) varying fastest
        k_cart : (N^2, 3) float, 1/Angstrom, k_red @ B.T
        K      : (3,) float, 1/Angstrom, Dirac point (2 b1 + b2)/3
    """

    u = (np.arange(N) + shift)/N # reduced coordinates along one axis, in (0, 1)

    # all (u1, u2) pairs: two (N, N) arrays, flattened to N^2 aligned entries
    u1, u2 = np.meshgrid(u,u, indexing='ij')
    u3 = np.zeros(N**2)

    k_red = np.column_stack((u1.ravel(), u2.ravel(), u3))
    k_cart = k_red @ B.T # rows of reduced coords -> rows of cartesian vectors

    K = (2*B[:,0] + B[:,1])/3 # Dirac point

    return k_red, k_cart, K

def make_grid_tb(tb, N):
    """
    GridTB of `tb` (only its lattice is used): reciprocal vectors, shifted N x N grid (`k_grid`)
    and Dirac point.
    """
    B = reciprocal(tb.lattice)
    k_red, k_cart, K = k_grid(B,N)

    return GridTB(B=B, k_red=k_red, k_cart=k_cart, K=K)

def dagger(X):
    """
    Hermitian conjugate over the last two axes; the leading axes (k, mu) are untouched.
    """

    return np.conj(np.swapaxes(X, -2, -1))

def hermitize(X):
    """
    Hermitian part (X + X^dagger)/2 over the last two axes.

    For A(k): the r(R) of Wannier90 are Hermitian only approximately (finite differences). Applied
    by the caller after `fourier`, so that the raw non-Hermiticity stays visible to the tests.
    """

    return (X + dagger(X))/2


def fourier(X_R, R_cart, ndegen, k, deriv=None):
    """
    Lattice-gauge Fourier transform to arbitrary k points, the single one of the module:

        X(k)        = sum_R e^{+ik.R} X(R) / ndegen(R)
        dX(k)/dk_mu = sum_R i R_mu e^{+ik.R} X(R) / ndegen(R)      (deriv=True)

    Used for H(k) (Nk, nW, nW) in eV, dH/dk (Nk, 3, nW, nW) in eV*Angstrom and A(k) (Nk, 3, nW, nW)
    in Angstrom. Same phase as Hwr_to_Hwk (k_cart . R_cart = 2 pi k_red . R_int). The phase matrix
    takes 16 Nk nR bytes (38 GB for 1800^2 k and 741 R): large grids go through in blocks of
    1e4-1e5 k points.

    Inputs:
        X_R    : (nR, ...) complex; only (nR, nW, nW) can be differentiated
        R_cart : (nR, 3) float, Angstrom
        ndegen : (nR,) int
        k      : (Nk, 3) float, 1/Angstrom (a single k point as k[None, :])
        deriv  : bool, return dX/dk instead of X(k)
    Returns:
        X_k : (Nk, ...) complex, or (Nk, 3, nW, nW) with deriv
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

    # flatten the trailing axes, sum over R as one matrix product, restore the trailing axes (C order)
    X_R = X_R.reshape(nR, -1) # (nR, M)
    X_k = phase @ X_R # (nK, M)
    X_k = X_k.reshape(Nk, *tail)

    return X_k

def velocity(H_k, dH_k, A_k=None):
    """
    Band energies, eigenvectors and hbar v_mn(k) in the band basis,

        hbar v^mu = V^dagger [ dH/dk_mu + i [H, A_mu] ] V      (Wang et al., PRB 74, 195118 (2006)).

    The Berry term i[H, A] accounts for where the Wannier functions actually are. It vanishes on the
    diagonal (hbar v_nn = d eps_n/dk either way) and makes the interband elements independent of the
    home cell of each Wannier function. No consistency check inside (they would cost a
    diagonalization per block; they live in the tests). Interband elements between degenerate bands
    (at K) are ill-defined: grids and rings avoid K.

    Inputs:
        H_k  : (Nk, nW, nW) complex, eV
        dH_k : (Nk, 3, nW, nW) complex, eV*Angstrom (not modified)
        A_k  : (Nk, 3, nW, nW) complex, Angstrom, Hermitian; None drops the Berry term
               (eq. (2.5.7) of the thesis as written)
    Returns:
        eps : (Nk, nW) float, eV, ascending at every k
        V   : (Nk, nW, nW) complex, eigenvectors in columns (Wannier -> band rotation); the phase of
              each column is arbitrary at every k, so only |hbar v_mn|^2 is comparable between calls
        hv  : (Nk, 3, nW, nW) complex, eV*Angstrom, Hermitian, mu on axis 1
              (hbar v_F = 5.751 eV*Angstrom <-> v_F = 8.7e5 m/s)
    """
    nk, nW, _ = H_k.shape

    # diagonalize H(k) at every k at once: eps ascending, eigenvectors in the columns of V
    eps, V = np.linalg.eigh(H_k) # (nk, nW) floats, (nk, nW, nW) complex

    # velocity operator (times hbar) in the Wannier basis, W = dH + i [H, A];
    # copy(): W = dH_k would alias the caller's array and the += below would overwrite it
    W = dH_k.copy() # (nk, 3, nW, nW) complex
    if A_k is not None:
        assert A_k.shape == (nk, 3, nW, nW), 'wrong shape for A_k'
        # H_k[:, None]: same H(k) for the three components A_mu; i makes the commutator Hermitian
        W += 1j*(H_k[:, None, ...] @ A_k - A_k @ H_k[:, None, ...]) # add berry connection term

    assert W.shape == (nk, 3, nW, nW), 'wrong shape for W'

    # rotation to the band basis, with the same V(k) for the three components mu
    hv = dagger(V)[:, None, ...] @ W @ V[:, None, ...] # (nk, 3, nW, nW) complex
    assert hv.shape == (nk, 3, nW, nW), 'wrong shape for hv'

    return eps, V, hv

def _gap(tb, K, q, e, mu):
    """
    Gap eps_c - eps_v at k = K + q e along each ray, for trial distances q: the function whose
    root `ring` looks for. c and v are the bands just above and below mu, found by energy (fixed
    indices would fail with the 5-band data). Never called at q = 0, where the bands are degenerate
    at mu.

    Inputs:
        tb : WannierTB
        K  : (3,) float, 1/Angstrom
        q  : (ntheta,) float, 1/Angstrom, trial distance per ray
        e  : (ntheta, 3) float, unit vectors of the rays
        mu : float, eV
    Returns:
        gap : (ntheta,) float, eV
    """

    # one k point per ray
    k = K[None, :] + q[:, None] * e # (ntheta, 3)
    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, k) # (ntheta, nW, nW)
    eps = np.linalg.eigvalsh(H_k) # (ntheta, nW), ascending
    n_occ = np.sum(eps < mu, axis=1) # (ntheta,), number of occupied states (with energy less than chemical potential) per k
    v = n_occ - 1 # index of valence band
    c = n_occ  # index of conduction band

    # eps[i, v[i]] row by row; [:, 0] turns (ntheta, 1) into (ntheta,)
    eps_v = np.take_along_axis(eps, v[:, None], axis=1)[:,0]
    eps_c = np.take_along_axis(eps, c[:, None], axis=1)[:,0]

    gap = eps_c - eps_v # (ntheta, )

    return gap

def ring(tb, K, hw, mu=0.0, ntheta=720, tol=1e-12, maxsteps=100):
    """
    Resonant ring around K: one point per direction theta where eps_c - eps_v = hw,

        k(theta) = K + q(theta) (cos theta, sin theta, 0).

    These states absorb light of energy hw exactly (no broadening); they serve for the angular
    checks of the Berry term, the ring averages of M3, the figure of section 2.5 and the k list of
    EM2. A circle of radius q0 = hw / (2 hbar v_F) for a perfect cone, a rounded triangle with
    trigonal warping. Vectorized bisection on [0, 2 q0], valid while hw stays well below the gap 2t
    at M (checked by the bracket assert): ~39 steps, gap exact to ~1e-11 eV.

    Inputs:
        tb       : WannierTB (t and a_cc only set q0)
        K        : (3,) float, 1/Angstrom, centre of the ring (Dirac point)
        hw       : float, eV, photon energy
        mu       : float, eV, chemical potential (0 in M0, E_D in M1)
        ntheta   : int, multiple of 6 (theta + pi and theta + 120 degrees stay on the grid)
        tol      : float, 1/Angstrom, final bisection width
        maxsteps : int, cap on the number of bisection steps
    Returns:
        k     : (ntheta, 3) float, 1/Angstrom
        theta : (ntheta,) float, rad, equally spaced, endpoint excluded (plain means are angular averages)
        q     : (ntheta,) float, 1/Angstrom, |k - K|
        q0    : float, 1/Angstrom, Dirac estimate hw / (2 hbar v_F)
    """
    assert ntheta % 3 == 0 and ntheta % 2 == 0, 'number of angles must be a multiple of 2 and 3'

    # directions: unit vectors (cos theta, sin theta, 0), one per row
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False) # (ntheta,)
    e = np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))) # (ntheta, 3)

    # bracket: q = 0 (gap 0 < hw, known, not evaluated) and twice the Dirac estimate
    hv_F = 3*tb.t * tb.a_cc / 2
    q0 = hw / (2*hv_F) # float
    lo = np.zeros(ntheta); hi = 2*q0*np.ones(ntheta) # (ntheta, )
    assert np.all(_gap(tb, K, hi, e, mu) > hw ), 'hw too large: gap below hw at 2 q0 on some ray (ring not closed?)'

    # bisection on all rays at once: gap(mid) < hw -> crossing further out (lo <- mid), else hi <- mid
    nsteps = 0
    while (hi - lo).max() > tol and nsteps < maxsteps:
        mid = (hi + lo) / 2
        gap_mid = _gap(tb, K, mid, e, mu)
        cond = gap_mid < hw # (ntheta,) bool

        lo = np.where(cond, mid, lo); hi = np.where(cond, hi, mid)
        nsteps += 1
    assert nsteps < maxsteps, 'bisection did not converge within maxsteps'

    q = (hi + lo) / 2
    k = K[None, :] + q[:, None]*e # (ntheta, 3)

    assert k.shape == (ntheta, 3), 'wrong shape for k'
    assert theta.shape == (ntheta,), 'wrong shape for theta'
    assert q.shape == (ntheta,), 'wrong shape for q'

    return k, theta, q, q0
