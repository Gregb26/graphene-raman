import numpy as np
from electron_defect_interaction.io import qe_io
from electron_defect_interaction.defects.many_body.single_defect import compute_T
HA2EV = 27.211386245988


def medgamma(N, norm):
    uc = f"data/graphene/unit_cell/qe/defect_{N}.save"
    M = np.load(f"results/M/M_ed_{N}.npy")
    eigs = qe_io.get_eigenvalues(uc, shift_Fermi=True) * HA2EV
    nb, nk = eigs.shape
    NN = nb * nk
    Mn = M / nk if norm else M
    ef = eigs.reshape(NN)
    g = np.zeros(NN)
    for s in range(0, NN, 32):
        e = ef[s:s + 32]
        T = compute_T(e, eigs, Mn, 0.02).reshape(len(e), NN, NN)
        for j in range(len(e)):
            g[s + j] = -2 * T[j, s + j, s + j].imag
    return np.median(np.abs(g)) * 1e3


print(f"{'N':>5} {'Gamma raw (meV)':>18} {'Gamma M/nk (meV)':>18}")
for N in ["5x5", "6x6", "7x7", "8x8"]:
    print(f"{N:>5} {medgamma(N, False):>18.3f} {medgamma(N, True):>18.4f}", flush=True)
