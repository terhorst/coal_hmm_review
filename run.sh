#!/bin/bash -ex
OUTPUT_DIRECTORY=/scratch/terhorst/benchmark/output
rm -rf $OUTPUT_DIRECTORY/*
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
     EstimateSizeHistoryDical \
     --EstimateSizeHistoryDical-N 10 \
     --EstimateSizeHistoryDical-demography constant \
     --EstimateSizeHistoryDical-seed 1 \
     --GlobalConfig-chromosome-length 1000000 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 24 --local-scheduler $@
