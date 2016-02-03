

#!/bin/bash
#
#$ -S /bin/bash
#$ -o /netapp/home/dreuxj/Rando                      
#$ -e /netapp/home/dreuxj/Rando                   
#$ -r y                                 
#$ -j y                                 
#$ -l mem_free=5G                    
#$ -l arch=linux-x64                    
#$ -l netapp=1G,scratch=1G              
#$ -l h_rt=24:00:00                     
#$ -t 1                



echo "Job ID is:" $JOB_ID
echo "SGE Task ID:" $SGE_TASK_ID
hostname
date

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA462/ERA462794/fastq/Quiescent_1-P1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA462/ERA462794/fastq/Quiescent_1-P2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA462/ERA462794/fastq/Quiescent_2-P1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA462/ERA462794/fastq/Quiescent_2-P2.fastq.gz

date

qstat -j $JOB_ID