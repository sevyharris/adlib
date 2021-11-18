#!/bin/bash

#SBATCH --time=unlimited
#SBATCH --job-name=QE_bulk_relax
#SBATCH --partition=west
#SBATCH --mem=20Gb


start=$SECONDS
python relax.py | tee bash_log.txt

duration=$(( SECONDS - start))
echo "Completed in $duration seconds" | tee -a bash_log.txt

