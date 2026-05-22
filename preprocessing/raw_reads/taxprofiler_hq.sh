#!/bin/bash
#SBATCH --job-name=taxprofiler
#SBATCH --time=72:00:00
#SBATCH --partition=large        #for multinode, large partition should be use
#SBATCH --nodes=10                #max nodes allowed 26, increase here if needed for faster
#SBATCH --account=project_2010455
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50G       #increase if needed, reduce if too much. check with seff -j <jobID>


export SINGULARITY_TMPDIR=$PWD
export SINGULARITY_CACHEDIR=$PWD
unset XDG_RUNTIME_DIR

# Activate  Nextflow on Puhti
module load nextflow
module load hyperqueue

# Create a per job directory

export HQ_SERVER_DIR=$PWD/.hq-server-$SLURM_JOB_ID
mkdir -p $HQ_SERVER_DIR

hq server start &
srun --cpu-bind=none --hint=nomultithread --mpi=none -N $SLURM_NNODES -n $SLURM_NNODES -c 4 hq worker start --cpus=$SLURM_CPUS_PER_TASK &

num_up=$(hq worker list | grep RUNNING | wc -l)
while true; do

    echo "Checking if workers have started"
    if [[ $num_up -eq $SLURM_NNODES ]];then
        echo "Workers started"
        break
    fi
    echo "$num_up/$SLURM_NNODES workers have started"
    sleep 1
    num_up=$(hq worker list | grep RUNNING | wc -l)

done

#run the nextflow 
#for testing (in profile params use test,singularity flag )

nextflow run nf-core/taxprofiler -r 1.1.5 -c /scratch/project_2010455/config/nextflow.config -resume -with-report /scratch/project_2010455/USERS/report_taxprofiler.html \
   -profile singularity \
   --input /scratch/project_2010455/config/samplesheet.csv \
   --databases /scratch/project_2010455/config/database.csv \
   --outdir ./RESULTS  \
   --perform_shortread_qc \
   --perform_shortread_complexityfilter  --shortread_complexityfilter_tool bbduk \
   --perform_shortread_hostremoval \
   --hostremoval_reference /projappl/project_2009677/DB/T2T-CHM13v2.0.zip \
   --shortread_hostremoval_index /projappl/project_2009677/DB/human_CHM13/ \
   --perform_runmerging \
   --save_analysis_ready_fastqs \
   --save_hostremoval_unmapped \
   --run_profile_standardisation \
   --run_metaphlan


# Make sure we exit cleanly once nextflow is done
hq worker stop all
hq server stop    