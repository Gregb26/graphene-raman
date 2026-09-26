"""
deltav_pw.py
    ROLE: direct application of the defect potential dV = V_d - V_p (local part on the supercell FFT grid, non-local part
    through the Kleinman-Bylander projectors of the removed atom) to Gamma-point supercell states given on the supercell
    plane-wave / FFT grid. This is the reference against which every production M is gated (R6 gate A.2, scripts/
    gate_M_normalization.py): for pure folded Bloch states <nk|dV|n'k'>_sc must equal M[nk,n'k'] / N_cells (L and NL
    separately, <= 1e-6 eV). Promoted from the R5 campaign module r5_deltav_pw.py (graphene/qe/defects/R5_base_vs_M/,
    2026-09-25) into src/ by R6, unchanged except this header and the import of sc_projection.
    Same kernels and conventions as the production M chain:
      * local part    : compute_ML_R's real-space product  <bra|dV^L|ket> = sum_r conj(Psi_bra) dV Psi_ket dvol, dV = V_d - V_p
                        from pp.x plot_num 1 (qe_io.get_pot, subtract_mean=False as in production), dvol = Omega_sc / N_grid;
      * non-local part: compute_M_NL's Kleinman-Bylander sum  B_lim = 4 pi/sqrt(Omega) sum_g conj(D(g)) F_li(|g|) Y_lm(g_hat) e^{-i g.tau},
                        <bra|V^NL_atom|ket> = sum_li D_li B_bra conj(B_ket), evaluated with the SUPERCELL plane waves (k = 0, full
                        sphere), the supercell volume in the prefactor and ONLY the removed atom; M^NL = M_d - M_p => dV^NL = - V^NL_removed.
    With |nk>_sc = psi_nk / sqrt(N_cells) these expectation values equal M[nk, n'k'] / N_cells (unit-cell Bloch norm): check_pure_bloch.
"""
import numpy as np

from electron_defect_interaction.wavefunctions.sc_projection import grid_to_real, bloch_state_grid
from electron_defect_interaction.defects.non_local import build_K_vectors, compute_phase, compute_angular_part
from electron_defect_interaction.io.pseudo_io import read_upf, fq_from_fr


def expect_local(psi_bra, psi_ket, dV, Omega_sc):
    """<bra|dV|ket> = sum_r conj(psi_bra) dV psi_ket * dvol; psi_* real-space arrays from grid_to_real (unit norm), dV same shape."""
    dvol = float(Omega_sc) / dV.size
    return complex(np.vdot(psi_bra, dV * psi_ket) * dvol)


class RemovedAtomProjector:
    """
    KB projectors of one atom at tau_cart (Bohr) on the full supercell sphere mill_full (nfull, 3), B_sc (3, 3) columns b_i (1/Bohr).
    W[g, l, i, m] = F_li(|g|) Y_lm(g_hat) e^{-i g.tau}; projections(D) = pref * sum_g conj(D) W with pref = 4 pi / sqrt(Omega_sc);
    expect(B_bra, B_ket) = - sum_lim D_li B_bra conj(B_ket) (minus: the atom is removed in the defective cell).
    Uses build_K_vectors / compute_phase / compute_angular_part / read_upf / fq_from_fr of the production chain (nkpt = 1, natom = 1).
    """

    def __init__(self, mill_full, B_sc, tau_cart, upf_path, Omega_sc, ecut_Ha, nq=2000, energy_scale=1.0):
        """energy_scale multiplies expect(): 1.0 -> Hartree (read_upf's D_li), HA2EV -> eV."""
        from scipy.interpolate import CubicSpline
        ekb_li, fr_li, rgrid, lmax, imax, _ = read_upf(upf_path)
        nfull = len(mill_full)
        k_red = np.zeros((1, 3)); G_red = np.asarray(mill_full, float)[None]; keep = np.ones((1, nfull), bool)
        K, K_norm, K_hat = build_K_vectors(k_red, G_red, keep, np.asarray(B_sc, float))
        qmax = 2 * np.sqrt(2 * float(ecut_Ha)); q = np.linspace(0, qmax, nq)
        fq = fq_from_fr(rgrid, fr_li, q); Fq = CubicSpline(q, fq, axis=-1, extrapolate=False)
        if K_norm.max() >= qmax:
            raise ValueError("removed-atom projector: |g| beyond the q grid")
        F = Fq(K_norm)[:, :, 0, :]                                             # (l, i, nfull)
        ph = compute_phase(K, np.asarray(tau_cart, float)[None])[0, 0]         # (nfull,)
        Y = compute_angular_part(K_hat, lmax)[0]                               # (nfull, l, m)
        self.W = np.einsum("lig,glm,g->glim", F, Y, ph, optimize=True)         # (nfull, l, i, m)
        self.pref = 4 * np.pi / np.sqrt(float(Omega_sc)); self.ekb = ekb_li
        self.lmax = lmax; self.imax = imax; self.K_norm_max = float(K_norm.max()); self.qmax = qmax; self.escale = float(energy_scale)

    def projections(self, D):
        """D (ns, nfull) full-sphere coefficients -> B (ns, l, i, m)."""
        return self.pref * np.einsum("ng,glim->nlim", np.conj(np.atleast_2d(D)), self.W, optimize=True)

    def expect(self, B_bra, B_ket):
        """<bra|dV^NL|ket> = - energy_scale * sum_li D_li sum_m B_bra[l,i,m] conj(B_ket[l,i,m])."""
        return -self.escale * complex(np.einsum("li,lim,lim->", self.ekb, B_bra, np.conj(B_ket)))


def check_pure_bloch(pairs, C_nkg, nG, flat_idx, ngfft, Omega_sc, dV, M_L, M_NL, N_cells, proj, flat_full, workers=8):
    """
    For (n, k, n', k') in pairs: <nk|dV^L|n'k'>_sc and <nk|dV^NL|n'k'>_sc from the grid kernels against M_L[n,k,n',k']/N_cells and
    M_NL[...]/N_cells (M in the unit of dV). Returns a list of dicts with both values and the differences.
    """
    out = []
    for (n, k, n2, k2) in pairs:
        A1 = bloch_state_grid(C_nkg, nG, flat_idx, n, k, ngfft)
        same = (n2, k2) == (n, k)
        A2 = A1 if same else bloch_state_grid(C_nkg, nG, flat_idx, n2, k2, ngfft)
        p1 = grid_to_real(A1, ngfft, Omega_sc, workers); p2 = p1 if same else grid_to_real(A2, ngfft, Omega_sc, workers)
        vl = expect_local(p1, p2, dV, Omega_sc)
        B1 = proj.projections(A1[flat_full][None])[0]; B2 = B1 if same else proj.projections(A2[flat_full][None])[0]
        vnl = proj.expect(B1, B2)
        ml = complex(M_L[n, k, n2, k2]) / N_cells; mnl = complex(M_NL[n, k, n2, k2]) / N_cells
        out.append(dict(pair=(int(n), int(k), int(n2), int(k2)), VL_direct=vl, ML_over_N=ml, dL=abs(vl - ml),
                        VNL_direct=vnl, MNL_over_N=mnl, dNL=abs(vnl - mnl), norm1=float(np.vdot(A1, A1).real)))
    return out
