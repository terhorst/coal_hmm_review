#!/bin/bash
OUTPUT_DIRECTORY=/scratch/terhorst/tishkoff/pipeline.output
mkdir -p $OUTPUT_DIRECTORY
PYTHONPATH=. luigi --module convert_data \
    ConvertAll \
    --GlobalConfig-output-directory $OUTPUT_DIRECTORY \
    --workers 30 --local-scheduler
