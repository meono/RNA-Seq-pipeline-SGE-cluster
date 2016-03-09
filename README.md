#RNA-Seq pipeline on SGE cluster

A pipeline for RNA-Seq data processing using the UCSF qb3 cluster (SGE). 

#TOOLS
FASTQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

FASTx_toolkit - http://hannonlab.cshl.edu/fastx_toolkit/download.html

STAR - https://github.com/alexdobin/STAR

SAMtools - http://www.htslib.org/download/

HTSeq - http://www-huber.embl.de/HTSeq/doc/overview.html

Subread bioconductor package (featureCounts) - http://subread.sourceforge.net/ 

#STEPS

1. FASTQC - check read quality
2. FASTQ trimmer (fastx_toolkit) - trim low quality bases. Input is raw FASTQ files, output is (hopefully cleaner) FASTQ files
3. STAR aligner - use hg38 as a reference. Input is annotation file (I use GRCh38-Gencode24) and FASTQ files from step 1 or 2 Output is SAM file.
4. Convert the alignement sam file to a bam file, sort and index it using SAMTOOLS.
5. Count the reads with HTSeq or featureCounts (faster). Input is a sorted BAM file.
6. Locally analyze with EdgeR or DeSeq2 for differential expression


Note: the masterfile to qsub your jobs also supports the cufflinks suite (bowtie, tophat, cuffdiff...etc.) although they are no longer part of my pipeline.
