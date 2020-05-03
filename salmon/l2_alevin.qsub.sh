#!/bin/bash -l

# qsub options
#$ -P bf528 # project
#$ -pe mpi_16_tasks_per_node 16
#$ -l h_rt=48:00:00 # maximum run time
#$ -N l2_alevin # job name
#$ -j y # join stdout and stderr
#$ -o l2_alevin.qlog  # log file name

# job info
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""



salmon alevin -l ISR -1 SRR3879605_1_bc.fastq.gz -2 SRR3879605_2.fastq.gz --end 5 --barcodeLength 19 --umiLength 6  -i index -p 10 --whitelist library_2_whitelist.txt -o l2_alevin_output --tgMap tx2gene_ensemble.tsv --dumpMtx

echo "Analysis complete!"


