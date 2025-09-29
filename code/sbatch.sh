#!/bin/bash

module load R

Rscript integration.R

#sbatch --cpus-per-task=20 --mem=50g --time=24:00:00 sbatch.sh