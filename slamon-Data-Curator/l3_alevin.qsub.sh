#!/bin/bash -l

# qsub options
#$ -P bf528 # project
#$ -pe mpi_16_tasks_per_node 16
#$ -l h_rt=48:00:00 # maximum run time
#$ -N l3_alevin # job name
#$ -j y # join stdout and stderr
#$ -o l3_alevin.qlog  # log file name

# job info
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""



salmon alevin -l ISR -1 SRR3879606_1_bc.fastq.gz -2 SRR3879606_2.fastq.gz --end 5 --barcodeLength 19 --umiLength 6  -i index -p 10 --whitelist library_3_whitelist.txt -o l3_alevin_output --tgMap tx2gene_ensemble.tsv --dumpMtx

echo "Analysis complete!"


