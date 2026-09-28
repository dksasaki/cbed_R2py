#!/bin/bash
#SBATCH --job-name=cbed
#SBATCH --partition=sharing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --exclude=d3032,d3232,d3203

N_WORKERS=${1:-128}
THREADS_PER_WORKER=${2:-1}

pixi run python scripts/cbed_wrapper.py "$N_WORKERS" "$THREADS_PER_WORKER"
