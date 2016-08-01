#RNA-Seq pipeline for iLoop on Computerome

Basic functionality is implemented.

### Operation:

1. Quality check: FASTQC run. Output still has to be handled manually. Does not impede rest of the process.

2. Map and link jobs: HISAT2, cufflinks, stringtie and htseq-count. All are run by default but can be limited to requirements.

3. Merge job: Cuffmerge run. Depends on successful completion of step 2.

4. Quantify job: Cuffquant run. Depends on successful completion of step 3.

5. DGE job: CUffdiff run. Depends on successful completion of step 4.


## Operational requirements:

Package includes templates for "RNAseq_pipeline_defaults.txt" and "RNAseq_pipeline_references.tsv".

"RNAseq_pipeline_defaults.txt" file contains preferred options for the RNAseq pipeline components. The template contains defaults that will be used unless a file with the same name isn't found under project path or home path. Parameters in file under project path has precedence over file under home path, which in turn has precedence over defaults.

"RNAseq_pipeline_references.tsv" file contains necessary paths for reference annotations, sequences, indexes and so on. Currently, these are done manually as well. Same precedence pattern for "RNAseq_pipeline_defaults.txt" applies.

Optional: An environmental variable "RNAseq_env" for the virtual environment name (if installed as such) where iLoop_RNAseq_pipeline package is installed. 

##Requirements:

HISAT2 and htseq-count depends on local installations due to version issues. The rest of the packages are used through Computerome module system.

A local installation of featureCounts since the module version is old.

R packages:
- edgeR
- DESeq
- optparse