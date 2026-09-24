#!/bin/bash
#SBATCH --job-name=test_mnl_mpi
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=00:45:00
#SBATCH --output=results/M/logs/%x_%j.out
#SBATCH --error=results/M/logs/%x_%j.err

# #3 root-cause test: does the chunked-Allreduce fix make compute_M_NL_mpi work at 9x9
# (nk=81, the heap-corruption regime)? Run the MPI M^NL (faulthandler on), the serial
# M^NL, and compare. If the MPI step still aborts, faulthandler localizes the line.

PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}
cd "$PROJ" || exit 1
module restore qe
module load mpi4py/4.0.3 scipy-stack
PY="$PROJ/.venv/bin/python"

N=9x9
UC=data/graphene/unit_cell/qe/defect_${N}.save
SCP=data/graphene/supercell/qe/defect_${N}_p.save
SCD=data/graphene/supercell/qe/defect_${N}_d.save
OUT=results/M/_test_mnl
mkdir -p "$OUT"

echo "[$(date)] MPI M^NL (32 ranks, chunked Allreduce + faulthandler, FLEXIBLAS=IMKL)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# Production candidate: Intel MKL backend (fast AND robust under MPI). The default OpenBLAS
# backend corrupts the heap under srun at nk>=81; NETLIB avoids it but is slow.
export FLEXIBLAS=IMKL
srun --cpu-bind=cores "$PY" -u scripts/_diag_mnl_mpi.py \
    "$UC" "$SCP" "$SCD" "$UC/C.upf" "$OUT/M_NL_mpi.npy"
MPI_RC=$?
echo "[$(date)] MPI step exit code: $MPI_RC"

echo "[$(date)] serial M^NL (reference)"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 FLEXIBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
"$PY" -u scripts/compute_M.py --stage nl \
    --uc "$UC" --sc-p "$SCP" --sc-d "$SCD" --upf "$UC/C.upf" \
    --bands all --out "$OUT/M_NL_serial.npy"

echo "[$(date)] compare MPI vs serial"
"$PY" - <<'PYEOF'
import numpy as np, os
mpi_f = "results/M/_test_mnl/M_NL_mpi.npy"
ser_f = "results/M/_test_mnl/M_NL_serial.npy"
if not os.path.exists(mpi_f):
    print("RESULT: FAIL (MPI M^NL never written -> it crashed; see faulthandler trace above)")
    raise SystemExit(1)
a = np.load(mpi_f); b = np.load(ser_f)
print(f"mpi {a.shape}  serial {b.shape}")
amax = np.max(np.abs(a - b)); rel = amax / max(np.max(np.abs(b)), 1e-300)
ok = np.allclose(a, b, rtol=1e-6, atol=1e-9)
print(f"max|mpi-serial|={amax:.3e}  rel={rel:.3e}  allclose={ok}")
print("RESULT:", "PASS" if ok else "FAIL")
PYEOF
echo "[$(date)] done"
