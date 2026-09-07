"""Write `vacancy_sublattice` (A/B), `vacancy_s_red` and the A/B-equivalence note into every M sidecar of each size."""
import glob, json, numpy as np, sys
from electron_defect_interaction.io import qe_io
NOTE = "vacancy on sublattice {sub}; A and B are equivalent by the honeycomb inversion symmetry (unrelaxed isolated vacancy): same Gamma, M identical up to pz(A)<->pz(B) and a BZ rotation"
for S in (sys.argv[1:] or ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]):
    N = int(S.split("x")[0]); scp = f"data/graphene/supercell/qe/defect_{S}_p.save"; scd = f"data/graphene/supercell/qe/defect_{S}_d.save"
    xp = np.mod(qe_io.get_x_red(scp), 1.0); xd = np.mod(qe_io.get_x_red(scd), 1.0)
    dmin = np.array([np.min(np.linalg.norm(np.mod(xd - p + 0.5, 1) - 0.5, axis=1)) for p in xp]); s_vac = xp[int(np.argmax(dmin))]; s_uc = np.mod(N * s_vac, 1.0)
    sub = "A" if np.allclose(s_uc[:2], 1 / 3, atol=1e-3) else ("B" if np.allclose(s_uc[:2], 2 / 3, atol=1e-3) else "?")
    files = [f for f in glob.glob(f"results/M/M_*_{S}*.json") + glob.glob(f"results/M/M_{S}*.json") if "obsolete" not in f]
    for f in files:
        d = json.load(open(f)); d["vacancy_sublattice"] = sub; d["vacancy_s_red"] = [float(x) for x in s_vac]; d["vacancy_note"] = NOTE.format(sub=sub); json.dump(d, open(f, "w"), indent=2)
    print(f"[{S}] sublattice {sub}, s_red {np.round(s_vac, 5).tolist()} -> {len(files)} sidecars tagged")
