#!/bin/bash -ex
export PYTHONPATH=.
OUTPUT_DIRECTORY=/scratch/project/projectdirs/m2871/terhorst/coal_hmm_review/pipeline.output.jun
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
     --PlotAllCombined-N 50 \
     --PlotAllCombined-n-replicates 1 \
     --GlobalConfig-n-contigs 8 \
     --GlobalConfig-chromosome-length 125_000_000 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 12 --local-scheduler $@
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
