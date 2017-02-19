#!/bin/bash
OUTPUT_DIRECTORY=$PWD/pipeline.output
mkdir -p $OUTPUT_DIRECTORY
PYTHONPATH=. luigi --module estimate \
    EstimateSizeHistory \
    --EstimateSizeHistory-population Amhara \
    --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
    --workers 30 --local-scheduler
