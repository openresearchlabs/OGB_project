#!/bin/bash

# Define the directory containing the fastq.gz files
FASTQ_DIR="/scratch/project_2010455/USERS/nitin/5FINALRUN/RESULTS/analysis_ready_fastqs"

# Use globbing to list all matching files and extract the sample names
# Note: This assumes that all fastq.gz files have the format "{sample}.unmapped_1.fastq.gz"
# If the file naming convention is different, you may need to adjust the pattern accordingly
SAMPLES=("$FASTQ_DIR"*_1.fastq.gz)
module load humann/3.8

# Normalize
renorm_cpm() {
    input="RESULTS/humann3/raw/${sample}_1_genefamilies.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_cpm.tsv"
    humann_renorm_table --input "$input" --output "$output" --units cpm --update-snames
}

renorm_relab() {
    input="RESULTS/humann3/raw/${sample}_1_genefamilies.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_relab.tsv"
    humann_renorm_table --input "$input" --output "$output" --units relab --update-snames
}

# Regroup
regroup_cpm() {
    input="RESULTS/humann3/processed/${sample}_1_genefamilies_cpm.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_cpm_regroup.tsv"
    humann_regroup_table --input "$input" --output "$output" --groups uniref90_rxn
}

regroup_relab() {
    input="RESULTS/humann3/processed/${sample}_1_genefamilies_relab.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_relab_regroup.tsv"
    humann_regroup_table --input "$input" --output "$output" --groups uniref90_rxn
}

# Rename
rename_cpm() {
    input="RESULTS/humann3/processed/${sample}_1_genefamilies_cpm_regroup.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_cpm_regroup_rename.tsv"
    humann_rename_table --input "$input" --output "$output" --names metacyc-rxn
}

rename_relab() {
    input="RESULTS/humann3/processed/${sample}_1_genefamilies_relab_regroup.tsv"
    output="RESULTS/humann3/processed/${sample}_1_genefamilies_relab_regroup_rename.tsv"
    humann_rename_table --input "$input" --output "$output" --names metacyc-rxn
}


for filepath in "${SAMPLES[@]}"; do
    # Extract the sample name from the file name
    sample=$(basename "$filepath" _1.fastq.gz)
    renorm_cpm
    renorm_relab
    regroup_cpm
    regroup_relab
    rename_cpm
    rename_relab
    # Print the extracted sample name
    echo "$sample"
done


##MERGE
humann_join_tables --input RESULTS/humann3/processed/ --output RESULTS/humann3/final/genefamilies_relab_regroup.txt --file_name genefamilies_relab_regroup.tsv
humann_join_tables --input RESULTS/humann3/processed/ --output RESULTS/humann3/final/genefamilies_cpm_regroup.txt --file_name genefamilies_cpm_regroup.tsv
humann_join_tables --input RESULTS/humann3/raw/ --output RESULTS/humann3/final/pathabundance.txt --file_name _pathabundance
humann_join_tables --input RESULTS/humann3/raw/ --output RESULTS/humann3/final/pathcoverage.txt --file_name _pathcoverage
humann_join_tables --input RESULTS/humann3/processed/ --output RESULTS/humann3/final/genefamilies_relab.txt --file_name genefamilies_relab.tsv
humann_join_tables --input RESULTS/humann3/processed/ --output RESULTS/humann3/final/genefamilies_cpm.txt --file_name genefamilies_cpm.tsv
humann_join_tables --input RESULTS/humann3/raw/ --output RESULTS/humann3/final/genefamilies.txt --file_name _genefamilies


#split stratification
humann_split_stratified_table --input RESULTS/humann3/final/genefamilies_relab.txt --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
humann_split_stratified_table --input RESULTS/humann3/final/genefamilies_cpm.txt --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
humann_split_stratified_table --input RESULTS/humann3/final/pathcoverage.txt --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
humann_split_stratified_table --input RESULTS/humann3/final/pathabundance.txt --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
humann_split_stratified_table --input RESULTS/humann3/final/genefamilies_relab_regroup.txt  --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
humann_split_stratified_table --input RESULTS/humann3/final/genefamilies_cpm_regroup.txt  --output /scratch/project_2009677/RESULTS/humann3/processed/split_strat
#only extract UNIREF without organism
