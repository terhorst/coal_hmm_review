#!/bin/bash
INPUT_DIRECTORY=$PWD/input.data
OUTPUT_DIRECTORY=$PWD/pipeline.output
mkdir -p $OUTPUT_DIRECTORY
PYTHONPATH=. luigi --module estimate \
    EstimateSizeHistory \
    --EstimateSizeHistory-population Amhara \
    --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
    --GlobalConfig-input-directory $INPUT_DIRECTORY \
    --workers 1 --local-scheduler $@
