#!/bin/bash
#SBATCH --job-name=cbed
#SBATCH --partition=sharing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --exclude=d3032,d3232,d3203
pixi run python scripts/cbed_wrapper.py 192 5
