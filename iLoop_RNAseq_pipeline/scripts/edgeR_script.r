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
                               help="Output path + prefix",
                               metavar="character")
                  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# set project path
setwd(opt$project_path)

# write edgeR plots on this file
pdf(paste(opt$output, ".pdf"), sep="")

counts <- read.delim(opt$counts, header=T, row.names=1)

## Make design matrix
conditions <- read.csv(opt$strains)
Group <- factor(paste(conditions$Strain,conditions$Treatment,sep="."))
cbind(conditions,Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
y <- DGEList(counts=counts)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

fit <- glmFit(y, design)
## As example, use the line below
# my.contrasts <- makeContrasts(comparison1 = condition1-condition2, comparison2 = condition3-condition4, crossc_comparison = (condition1-condition2)-(condition3-condition4), levels=design)
my.contrasts <- makeContrasts(COMPARISONS_STRING, levels=design)

wb <- createWorkbook()
saveWorkbook(wb, 'edgeR_results.xlsx')
write.xlsx(my.contrasts, 'edgeR_results.xlsx', sheetName ="my.contrasts")

for (i in colnames(my.contrasts)){
	lrt = glmLRT(fit, contrast=my.contrasts[,i])
	assign(i, topTags(lrt, n=nrow(y))$table)
	write.xlsx(get(i), 'edgeR_results.xlsx', sheetName = i, append = TRUE)
  }

## write.csv(etable, 'results.csv')

## to call the excel writer/wrapper function - which seems to have a problem, so ignore it
## source('saveXlsx.r')
## save.xlsx('edgeR_results.xlsx', my.contrasts, get(colnames(my.contrasts)[1]), get(colnames(my.contrasts)[2]), get(colnames(my.contrasts)[3]), get(colnames(my.contrasts)[4]))