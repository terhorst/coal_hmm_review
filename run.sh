#!/bin/bash -ex
OUTPUT_DIRECTORY=/scratch/terhorst/benchmark/output
rm -rf $OUTPUT_DIRECTORY/*
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
luigi --module tasks \
     EstimateSizeHistoryMSMC \
     --EstimateSizeHistoryMSMC-N 20 \
     --EstimateSizeHistoryMSMC-demography constant \
     --EstimateSizeHistoryMSMC-seed 1 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 24 --local-scheduler $@
# luigi --module tasks \
#      EstimateManyReplicates \
#      --EstimateManyReplicates-N 20 \
#      --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#      --workers 24 --local-scheduler $@
