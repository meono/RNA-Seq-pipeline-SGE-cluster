## R script to work out IL_fermentation data

## ... using edgeR

library(splines)
library(limma)
library(edgeR)
library(xlsx)

setwd("EDGER_PATH")

# write edgeR plots on this file
pdf("edgeR_plots.pdf")

counts <- read.delim("counts_collected_Rready.txt", header=T, row.names=1)

## Make design matrix
conditions <- read.csv('conditions_Rready.csv')
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