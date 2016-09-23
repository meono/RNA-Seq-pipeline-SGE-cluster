#!/usr/bin/env Rscript

# Based on https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Module4_Part2_cummeRbund.R

library("optparse")
library("cummeRbund")

option_list = list(
                   make_option(c("-p", "--project-path"),
                               type="character",
                               default=NULL,
                               help="Project path",
                               metavar="character"),
                   make_option(c("-o", "--output"),
                               type="character",
                               default=NULL,
                               help="Output path",
                               metavar="character")
                  );

opt <- parse_args(OptionParser(option_list=option_list), convert_hyphens_to_underscores=TRUE);
# set project path
setwd(opt$project_path)
# create results path
if (!(dir.exists(opt$output))) {dir.create(opt$output, recursive=TRUE)}

refCuffdiff=paste(opt$project_path, "cdiff", "diff_out", sep="/")
gtfFilePath=paste(opt$project_path, "cmerge", "merged_asm", "merged.gtf", sep="/")
outfile=paste(opt$output, "cummeRbund_output.pdf", sep="/")

# read in Cufflinks output
cuff <- readCufflinks(dir=refCuffdiff, rebuild=TRUE)

#Set pdf device
pdf(file=outfile)

# Plot #1 - A dispersion of FPKM within samples
isoforms.disp<-dispersionPlot(isoforms(cuff))
isoforms.disp + theme(legend.position="none")

# Plot #2 - A density plot of FPKM across samples
isoforms.dens<-csDensity(isoforms(cuff))
isoforms.dens + xlim(c(-10, 10))

# Plot #3a - A box plot of FPKM across samples
isoforms.boxP<-csBoxplot(isoforms(cuff))
isoforms.boxP + theme(legend.position="none")

# Plot #3b - A box plot of FPKM across sample replicates
isoforms.boxP.rep<-csBoxplot(isoforms(cuff))
isoforms.boxP.rep + theme(legend.position="none")

# Plot #4 - A scatter-plot matrix for all samples
isoforms.scaPmat<-csScatterMatrix(isoforms(cuff))
isoforms.scaPmat

# Plot #5 - Volcano plot matrix across samples
isoforms.volP<-csVolcanoMatrix(isoforms(cuff))
isoforms.volP

# Plot #6a - Using k-means clustering a dendrogram of the distance between samples
isoforms.dend<-csDendro(isoforms(cuff))
isoforms.dend

# Plot #6b - Using k-means clustering a dendrogram of the distance between sample replicates
isoforms.dend.rep<-csDendro(isoforms(cuff), replicate=T)
isoforms.dend.rep

# Plot #7a - A heatmap of sample distace based on JS distance
isoforms.DistHeat<-csDistHeat(isoforms(cuff))
isoforms.DistHeat

# Plot #7b - A heatmap of sample replicate distace based on JS distance
isoforms.DistHeat.rep<-csDistHeat(isoforms(cuff), replicate=T)
isoforms.DistHeat.rep

# Plot #8a - Principal Component Analysis of all genes across each sample
isoforms.PCA<-PCAplot(isoforms(cuff),"PC1","PC2")
isoforms.PCA + theme(legend.position="none")

# Plot #8b - Principal Component Analysis of all genes across each sample replicate
isoforms.PCA.rep<-PCAplot(isoforms(cuff),"PC1","PC2", replicate=T, showPoints=F)
isoforms.PCA.rep + theme(legend.position="none")

# Plot #8c - Principal Component Analysis of all genes across each sample
isoforms.PCA_23<-PCAplot(isoforms(cuff),"PC2","PC3", showPoints=F)
isoforms.PCA_23  + theme(legend.position="none")

# Plot #8d - Principal Component Analysis of all genes across each sample replicate
isoforms.PCA_23.rep<-PCAplot(isoforms(cuff),"PC2","PC3", replicate=T, showPoints=F)
isoforms.PCA_23.rep + theme(legend.position="none")

# Plot #9a - MDS scaling of all genes across samples
isoforms.MDS<-MDSplot(isoforms(cuff))
isoforms.MDS + theme(legend.position="none")

# Plot #9b - MDS scaling of all genes across sample replicates
isoforms.MDS.rep<-MDSplot(isoforms(cuff),replicates=T)
isoforms.MDS.rep + theme(legend.position="none")

#Close pdf device - necessary before you can open it in your browser
dev.off()
