#!/bin/bash
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=pp.x
#SBATCH --mail-type=NONE          # Mail events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G               # memory; default unit is megabytes
#SBATCH --time=00-01:00          # time (DD-HH:MM)
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module restore qe

export OMP_NUM_THREADS=1

srun pp.x < pp.in > pp.out 
