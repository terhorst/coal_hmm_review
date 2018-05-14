#!/bin/bash -ex
#SBATCH -C haswell
#SBATCH -A m2871
#SBATCH -q premium
#SBATCH -N 8
#SBATCH -t 12:00:00
export PYTHONPATH=.
OUTPUT_DIRECTORY=/project/projectdirs/m2871/terhorst/coal_hmm_review/pipeline.output
# rm -rf $OUTPUT_DIRECTORY
mkdir -p $OUTPUT_DIRECTORY
HPC=1 luigi --module tasks \
     PlotAllCombined \
     --PlotAllCombined-N 100 \
     --PlotAllCombined-n-replicates 20 \
     --GlobalConfig-n-contigs 8 \
     --GlobalConfig-chromosome-length 125_000_000 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 512 $@
