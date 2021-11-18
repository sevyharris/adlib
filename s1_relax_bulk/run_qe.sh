#!/bin/bash

#SBATCH --time=unlimited
#SBATCH --job-name=QE_bulk_relax
#SBATCH --partition=west
#SBATCH --mem=20Gb

module list

start=$SECONDS
python relax.py

duration=$(( SECONDS - start))
echo "Completed in $duration seconds"

