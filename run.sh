#!/bin/bash -ex
INPUT_DIRECTORY=/scratch/terhorst/tishkoff
OUTPUT_DIRECTORY=$HOME/Dropbox.new/Dropbox/Berkeley/Research/african_tishkoff/pipeline.output
export PYTHONPATH=.
mkdir -p $OUTPUT_DIRECTORY
luigi --module tasks \
     ConvertAllSMC \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --GlobalConfig-input-directory $INPUT_DIRECTORY \
     --workers 24 --local-scheduler $@
luigi --module tasks \
     EstimateAllSizeHistories \
     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
     --GlobalConfig-input-directory $INPUT_DIRECTORY \
     --workers 2 --local-scheduler $@
# luigi --module estimate.treemix \
#     TreeMixData \
#     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#     --GlobalConfig-input-directory $INPUT_DIRECTORY \
#     --workers 24 --local-scheduler $@
# luigi --module tasks PairwiseMomiAnalysisFromOriginalData \
#     --populations '["Amhara", "Dizi"]' \
#     --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
#     --GlobalConfig-input-directory $INPUT_DIRECTORY \
#     --workers 1 --local-scheduler $@
