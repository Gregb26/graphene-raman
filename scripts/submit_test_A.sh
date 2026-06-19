#!/bin/bash
#SBATCH --job-name=test_A_ks
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Test A (KS Hamiltonian reconstruction) for all unit cells, saving per-cell
# deviation arrays diag(H)-eps to results/test_A/ for plotting.

PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1

# Module stack that provides h5py (hdf5-mpi), scipy, netCDF4 and mpi4py for the venv.
# mpi4py must be loaded before the venv python is used.
module restore qe
module load mpi4py/4.0.3 scipy-stack

# Single process, threaded BLAS over the allocated cores.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

"$PROJ/.venv/bin/python" -u scripts/run_test_A_batch.py
