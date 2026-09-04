"""Level-1 summary from specwd logs: R_cut convergence, grid/eta plateau, resonance. Usage: python scripts/summarize_level1_maps.py 5x5:JOB 7x7:JOB ..."""
import re, sys, numpy as np
maps = {}
for arg in sys.argv[1:]:
    S, J = arg.split(":"); rows = {}
    for line in open(f"results/M/logs/specwd_{J}.out"):
        m = re.match(r"^\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s*$", line)
        if m: rows[(int(m[1]), int(m[2]), float(m[3]))] = (float(m[4]), float(m[5]))
    maps[S] = rows
sizes = list(maps)
print("median Gamma*N_cells (meV) vs R_cut at (grid 240, eta 0.02) | (grid 120, eta 0.01)")
print(f"{'size':>5} " + " ".join(f"{'Rc'+str(r):>8}" for r in range(4)) + "   rel.change 0->1, 1->2, 2->3   |  " + " ".join(f"{'Rc'+str(r):>8}" for r in range(4)))
for S in sizes:
    a = [maps[S].get((r, 240, 0.02), (np.nan,))[0] for r in range(4)]; b = [maps[S].get((r, 120, 0.01), (np.nan,))[0] for r in range(4)]
    ch = [f"{(a[i+1]-a[i])/a[i]*100:+.1f}%" for i in range(3)]
    print(f"{S:>5} " + " ".join(f"{x:8.2f}" for x in a) + "   " + ", ".join(ch) + "   |  " + " ".join(f"{x:8.2f}" for x in b))
print("\ngrid/eta plateau at R_cut 3: max rel change when grid doubles (120->240) and when eta halves (0.02->0.01), per size")
for S in sizes:
    g = max(abs(maps[S][(3, 240, e)][0] - maps[S][(3, 120, e)][0]) / maps[S][(3, 240, e)][0] for e in (0.05, 0.02, 0.01))
    h = max(abs(maps[S][(3, N, 0.01)][0] - maps[S][(3, N, 0.02)][0]) / maps[S][(3, N, 0.02)][0] for N in (120, 240))
    print(f"{S:>5}  grid 120->240: {g*100:.2f}%   eta 0.02->0.01: {h*100:.2f}%")
print("\nresonance E-E_D (eV) reported at R_cut 3 (max on-shell rate within +-1.5 eV of E_D):")
for S in sizes:
    print(f"{S:>5} " + " ".join(f"{maps[S][(3, N, e)][1]:6.3f}" for N in (60, 120, 240) for e in (0.05, 0.02, 0.01)))
print("\nLevel 2 (fixed R_cut 3, grid 240, eta 0.02): Gamma*N_cells vs N")
for S in sizes: print(f"{S:>5} {maps[S][(3, 240, 0.02)][0]:8.2f} meV")
