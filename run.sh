#!/bin/bash -ex
OUTPUT_DIRECTORY=/scratch/terhorst/benchmark/output3
# rm -rf $OUTPUT_DIRECTORY   ## THIS IS VERY DANGEROUS
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
# luigi --module tasks \
#      EstimateManyReplicates \
#      --EstimateManyReplicates-N 20 \
#      --EstimateManyReplicates-n-replicates 10 \
#      --GlobalConfig-chromosome-length 100000000 \
#      --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#      --workers 24 --local-scheduler $@
luigi --module tasks \
     PlotDical \
     --PlotDical-N 10 \
     --PlotDical-seed 1 \
     --PlotDical-demography migration \
     --GlobalConfig-n-contigs 1 \
     --GlobalConfig-chromosome-length 10_000_000 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 12 --local-scheduler $@
# luigi --module tasks \
#       EstimateSizeHistorySMC \
#       --EstimateSizeHistorySMC-seed 1 \
#       --EstimateSizeHistorySMC-N 10 \
#       --EstimateSizeHistorySMC-demography migration \
#       --GlobalConfig-n-contigs 1 \
#       --GlobalConfig-chromosome-length 10_000_000 \
#       --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#       --workers 24 --local-scheduler $@

# cp $OUTPUT_DIRECTORY/*/*/*.pdf ~/Dropbox/plots/bench
