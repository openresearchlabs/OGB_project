#!/bin/bash

# Output file name
output_file="/scratch/project_2010455/RESULTS_15.6/reads_count_postprocess.tsv"

# Create header for the output file
echo -e "Sample_ID\tReads_Count" > "$output_file"

# Loop through both types of fastq.gz files
for file in /scratch/project_2010455/USERS/nitin/RESULTS1/analysis_ready_fastqs/*1.merged.fastq.gz /scratch/project_2010455/USERS/nitin/RESULTS1/analysis_ready_fastqs/*unmapped_1.fastq.gz; do
    # Check if the file exists
    if [[ -e "$file" ]]; then
        # Extract sample ID from file name
        if [[ "$file" == *"unmapped_1.fastq.gz" ]]; then
            full_sample_id=$(basename "$file" .unmapped_1.fastq.gz)
            sample_id=${full_sample_id%%_*}  # Extract the part before the first underscore
        else
            full_sample_id=$(basename "$file" 1.merged.fastq.gz)
            sample_id=${full_sample_id%%_*}  # Extract the part before the first underscore
        fi
        
        # Count the number of reads in the file (total lines divided by 4)
        reads_count=$(zcat "$file" | wc -l)
        reads_count=$((reads_count / 4))
        
        # Output sample ID and reads count to the output file
        echo -e "$sample_id\t$reads_count" >> "$output_file"
    fi
done

echo "Reads count for each sample has been saved to $output_file."
