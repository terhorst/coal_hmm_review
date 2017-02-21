#!/bin/bash
INPUT_DIRECTORY=/scratch/terhorst/tishkoff
OUTPUT_DIRECTORY=/scratch/terhorst/tishkoff/pipeline.output
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
luigi --module tasks \
    EstimateAllPairwise \
    --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
    --GlobalConfig-input-directory $INPUT_DIRECTORY \
    --workers 1 --local-scheduler $@
# luigi --module tasks PairwiseMomiAnalysisFromOriginalData \
#     --populations '["Amhara", "Dizi"]' \
#     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#     --GlobalConfig-input-directory $INPUT_DIRECTORY \
#     --workers 1 --local-scheduler $@
