#!/bin/bash -ex
OUTPUT_DIRECTORY=/scratch/terhorst/benchmark/output
# rm -rf $OUTPUT_DIRECTORY/*
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
# luigi --module tasks \
#      EstimateManyReplicates \
#      --EstimateManyReplicates-N 20 \
#      --EstimateManyReplicates-n-replicates 10 \
#      --GlobalConfig-chromosome-length 100000000 \
#      --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#      --workers 24 --local-scheduler $@
for demo in bottleneck recent_growth constant; do 
    luigi --module tasks \
         PlotCombined \
         --PlotCombined-N 50 \
         --PlotCombined-demography $demo \
         --PlotCombined-n-replicates 10 \
         --GlobalConfig-n-contigs 10 \
         --GlobalConfig-chromosome-length 10_000_000 \
         --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
         --workers 12 --local-scheduler $@
done

cp $OUTPUT_DIRECTORY/*/*/*.pdf ~/Dropbox/plots/bench
