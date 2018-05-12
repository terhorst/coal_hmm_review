#!/bin/bash -ex
#-SBATCH -C haswell
#SBATCH -A m2871
#SBATCH -q premium
#SBATCH -N 8
export PYTHONPATH=.
OUTPUT_DIRECTORY=/project/projectdirs/m2871/terhorst/coal_hmm_review/pipeline.output
# rm -rf $OUTPUT_DIRECTORY
mkdir -p $OUTPUT_DIRECTORY
# luigi --module tasks \
#      EstimateManyReplicates \
#      --EstimateManyReplicates-N 20 \
#      --EstimateManyReplicates-n-replicates 10 \
#      --GlobalConfig-chromosome-length 100000000 \
#      --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#      --workers 24 --local-scheduler $@
luigi --module tasks \
     PlotAllCombined \
     --PlotAllCombined-N 100 \
     --PlotAllCombined-n-replicates 1 \
     --GlobalConfig-n-contigs 10 \
     --GlobalConfig-chromosome-length 100_000_000 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 128 $@
# luigi --module tasks \
#      PlotDical \
#      --PlotDical-N 10 \
#      --PlotDical-seed 1 \
#      --PlotDical-demography migration \
#      --GlobalConfig-n-contigs 1 \
#      --GlobalConfig-chromosome-length 10_000_000 \
#      --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#      --workers 12 --local-scheduler $@

# cp $OUTPUT_DIRECTORY/*/*/*.pdf ~/Dropbox/plots/bench
