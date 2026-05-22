#!/bin/bash

module load snakemake/7.17.1
#module load snakemake
source ~/.bashrc
source ~/.bash_profile
module load humann/3.8 

mkdir -p ./RESULTS/humann3/
snakemake -s ./workflow/Snakefile_humann --cluster " sbatch --account=project_2009677 --time=72:00:00 \
--mem 150G --partition=small -e LOGS/human_err_%A_%a.txt -o LOGS/human_out_%A_%a.txt \
--ntasks=1 --cpus-per-task=32 --parsable" -j 43 \
 --keep-incomplete --keep-going --stats ./comp_jobs_stats_humann.json \
  --rerun-incomplete --use-envmodules \
   --latency-wait 60