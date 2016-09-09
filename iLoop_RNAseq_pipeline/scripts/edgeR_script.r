#!/usr/bin/env Rscript
library("optparse")
library(splines)
library(limma)
library(edgeR)
library(xlsx)

option_list = list(
                   make_option(c("-p", "--project-path"),
                               type="character",
                               default=NULL,
                               help="Project path",
                               metavar="character"),
                   make_option(c("-c", "--counts"),
                               type="character",
                               default=NULL,
                               help="Collected counts file ready for R.",
                               metavar="character"),
                   make_option(c("-s", "--strains"),
                               type="character",
                               default=NULL,
                               help="File describing strains and/or conditions.",
                               metavar="character"),
                   make_option(c("-o", "--output"),
                               type="character",
                               default=NULL,
                               help="Output path",
                               metavar="character"),
                   make_option(c("-f", "--filter"),
                               type="double",
                               default=0.5,
                               help="Filter threshold for CPM",
                               metavar="double")
                  );

opt <- parse_args(OptionParser(option_list=option_list), convert_hyphens_to_underscores=TRUE);
# set project path
setwd(opt$project_path)
# create results path
if (!(dir.exists(opt$output))) {dir.create(opt$output, recursive=TRUE)}

# write edgeR plots on this file
pdf(paste(opt$output, "edgeR_plots.pdf", sep="/"))

counts <- read.delim(opt$counts, header=T, row.names=1, sep=",")
counts <- counts[ , order(names(counts))]

# Plot mean-difference accross all pairs
# It's probably not worth doing this "per replicate" bases. It might make sense to see it "per averaged replicates" bases - which I haven't figured out how.

# box plot for counts
boxplot(log2(counts+1), las=1)
title('Gene level counts', ylab='log2(count+1)')

# bar plot for total counts
par(mar=c(5,6,4,2)+0.1,mgp=c(5,1,0)) # to avoid overlaping axis label/title
total_counts <- colSums(counts)
barplot(total_counts, las=1)
title('Total mapped read counts', ylab='Read counts')

## Make design matrix
conditions <- read.csv(opt$strains)
group <- factor(paste(conditions$Strain,conditions$Treatment,sep="."))
cbind(conditions,group=group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

## Make new DGEList
y <- DGEList(counts=counts, group=group)

## This is a manual filter to keep only Features that have a CPM>0.5 in at least two sample. Features with less are
## considered not expressed and won't be included in comparisons.
## CPM bound depends very much on the read/genome size. A CPM corresponding to 10 counts should be aimed for. Ideally,
## this shouldn't not be a static value. Since we usually get generous amount of reads, it will be set to 0.5 for now.
## This may not apply in all conditions, so a switch as an argument might be useful later on.
keep <- rowSums(cpm(y)>opt$filter) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

## Calculate normalization factor to account for library size
y <- calcNormFactors(y, method='upperquartile')

# Plot for MDS
plotMDS(y)
title('Fold-change based multi-dimensional scaling (MDS) plot')


## normalize by library size, and estimate dispersion allowing possible trend with average count size
y <- estimateDisp(y, design, robust=TRUE)

# Plot for MDS
plotMDS(y, method='BCV')
title('Biological coefficient of variation based multi-dimensional scaling (MDS) plot')

# plot BCV for replicates
plotBCV(y)
title('Biological coefficient of variation')

## Run an Exact Test for Differences between Two Negative Binomial Groups
#et <- exactTest(y)

fit <- glmQLFit(y, design, robust=TRUE)
# generate comparison string based in available groups to make pairwise comparisons for all
constr <- paste(combn(levels(group), 2, simplify=TRUE, FUN=function(x) {sprintf("%s_vs_%s=%s-%s", x[1], x[2], x[1], x[2])}), sep=', ')
con <- makeContrasts(contrasts=constr, levels=design)

# Create the excel file for results and design
wb <- createWorkbook()
saveWorkbook(wb, paste(opt$output, "edgeR_pairwise_comparisons.xlsx", sep='/'))
write.xlsx(con, paste(opt$output, "edgeR_pairwise_comparisons.xlsx", sep='/'), sheetName ="Comparisons", append=TRUE)
write.xlsx(counts, paste(opt$output, "edgeR_pairwise_comparisons.xlsx", sep='/'), sheetName ="Counts", append=TRUE)

for (i in colnames(con)){
	lrt = glmLRT(fit, contrast=con[,i])
	assign(i, topTags(lrt, n=nrow(y))$table)
	write.xlsx(get(i), paste(opt$output, "edgeR_pairwise_comparisons.xlsx", sep='/'), sheetName = i, append = TRUE)
	dt <- decideTestsDGE(lrt)
	isDE <- as.logical(dt)
	DEnames <- rownames(y)[isDE]
	plotSmear(lrt, de.tags=DEnames)
	abline(h=c(-1,1), col="blue")
	title(con[1])
}
dev.off()
