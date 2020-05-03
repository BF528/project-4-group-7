#!/bin/bash -l

# qsub options
#$ -P bf528 # project
#$ -pe mpi_16_tasks_per_node 16
#$ -l h_rt=48:00:00 # maximum run time
#$ -N index  # job name
#$ -j y # join stdout and stderr
#$ -o index.qlog  # log file name

# job info
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""


salmon index -i index -k 31 --gencode -p 4 -t gencode.v33.transcripts.fa.gz

echo "Analysis complete!"
