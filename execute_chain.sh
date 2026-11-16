#!/bin/bash
# Run all chunks one after another, then stitch them.
# Each chunk is `sbatch execute.sh N_WORKERS CHUNK` and starts when the
# previous one has finished. Run on the login node, from root_dir:
#
#   ./execute_chain.sh
#
# Chunks that already have an output file are skipped, so after a failure
# just run it again.

N_WORKERS=128
TIME=08:00:00       # per chunk (execute.sh asks for only 1 h)

# number of chunks = [chunks] nx * ny in config.toml (lines "nx = 5", "ny = 1")
NX=$(grep -E '^nx *=' config.toml | awk '{print $3}')
NY=$(grep -E '^ny *=' config.toml | awk '{print $3}')
N_CHUNKS=$((NX * NY))

mkdir -p logs
previous=""         # job id of the previous chunk
chunk_jobs=""       # all chunk job ids, for the stitch job

for chunk in $(seq 1 $N_CHUNKS); do
    output=$(printf 'data/cache/cbed_mom6_%03d.nc' $chunk)
    if [ -f "$output" ]; then
        echo "chunk $chunk: done already, skipped"
        continue
    fi

    # wait for the previous chunk (the first one starts right away)
    if [ -n "$previous" ]; then
        wait_for="--dependency=afterany:$previous"
    else
        wait_for=""
    fi

    job=$(sbatch --parsable $wait_for --time=$TIME \
                 --output=logs/chunk_${chunk}_%j.out \
                 execute.sh $N_WORKERS $chunk)
    echo "chunk $chunk of $N_CHUNKS: job $job"

    previous=$job
    chunk_jobs="$chunk_jobs:$job"
done

# stitch once every chunk has finished successfully
if [ -n "$chunk_jobs" ]; then
    job=$(sbatch --parsable --dependency=afterok$chunk_jobs \
                 --cpus-per-task=1 --mem=16G --time=01:00:00 \
                 --output=logs/stitch_%j.out \
                 --wrap="pixi run python scripts/cbed_stitch.py")
    echo "stitch: job $job"
fi
