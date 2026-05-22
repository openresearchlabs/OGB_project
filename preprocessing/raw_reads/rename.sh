#!/bin/bash

# Directory containing the fastq.gz files
FASTQ_DIR="/scratch/project_2010455/USERS/nitin/5FINALRUN/RESULTS/analysis_ready_fastqs"

# Loop over files and rename them
for file in "$FASTQ_DIR"/*.fastq.gz; do
    if [[ $file == *"_1.merged.fastq.gz" ]]; then
        # Rename files with '_1.merged.fastq.gz' to '{sample}_1.merged.fastq.gz'
        mv "$file" "${file/_1.merged.fastq.gz/.merged_1.fastq.gz}"
    fi
done