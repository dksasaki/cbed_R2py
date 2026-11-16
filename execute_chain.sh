#!/bin/bash
# Run all chunks one after another, then stitch them.
# Submits one chunk (`sbatch execute.sh N_WORKERS CHUNK`), waits until it has
# finished, then submits the next: only one job is in the queue at a time.
#
# This script keeps running until the last chunk is done, so start it where
# it can run that long, from root_dir, e.g.
#
#   srun --ntasks=1 --partition=short --time=24:00:00 ./execute_chain.sh
#
# Chunks that already have an output file are skipped, so after a failure
# (or if this script is stopped) just run it again.

N_WORKERS=128
TIME=01:00:00       # per chunk (execute.sh asks for 1 h); set from a measured chunk

# number of chunks = [chunks] nx * ny in config.toml (lines "nx = 50", "ny = 1")
NX=$(grep -E '^nx *=' config.toml | awk '{print $3}')
NY=$(grep -E '^ny *=' config.toml | awk '{print $3}')
N_CHUNKS=$((NX * NY))

mkdir -p logs

for chunk in $(seq 1 $N_CHUNKS); do
    output=$(printf 'data/cache/cbed_mom6_%03d.nc' $chunk)
    if [ -f "$output" ]; then
        echo "chunk $chunk: done already, skipped"
        continue
    fi

    echo "$(date '+%F %T') chunk $chunk of $N_CHUNKS: submitted"
    # --wait: sbatch returns only when the job has finished
    sbatch --wait --time=$TIME \
           --output=logs/chunk_${chunk}_%j.out \
           execute.sh $N_WORKERS $chunk

    if [ -f "$output" ]; then
        echo "$(date '+%F %T') chunk $chunk: finished"
    else
        echo "$(date '+%F %T') chunk $chunk: FAILED (no $output), see logs/chunk_${chunk}_*.out"
    fi
done

# stitch only if every chunk has its output file
missing=0
for chunk in $(seq 1 $N_CHUNKS); do
    [ -f "$(printf 'data/cache/cbed_mom6_%03d.nc' $chunk)" ] || missing=$((missing + 1))
done
if [ $missing -gt 0 ]; then
    echo "$missing chunk(s) missing: not stitching. Run this script again to redo them."
    exit 1
fi

echo "$(date '+%F %T') stitching"
sbatch --wait --cpus-per-task=1 --mem=16G --time=01:00:00 \
       --output=logs/stitch_%j.out \
       --wrap="pixi run python scripts/cbed_stitch.py"
echo "$(date '+%F %T') done: data/cbed_results.nc"
