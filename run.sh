#!/bin/bash -ex
OUTPUT_DIRECTORY=/scratch/terhorst/benchmark/output
rm -rf $OUTPUT_DIRECTORY/*
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
luigi --module tasks \
     EstimateManyReplicates \
     --EstimateManyReplicates-N 10 \
     --EstimateManyReplicates-n-replicates 10 \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --workers 24 --local-scheduler $@
