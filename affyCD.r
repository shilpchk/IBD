nam="GSE4183"
nam1 <- paste(nam,"_RAW.tar",sep="")
dy="/Users/shilpa/BIDMC/CD-Proj"
dy1 <- paste(dy,nam,sep="/")
setwd("/Users/shilpa/BIDMC/CD-Proj/Datasets/")

library(GEOquery)
getGEOSuppFiles("GSE20881")

untar(paste(dy1,nam1,sep="/"), exdir=paste(nam,"data",sep="/"))


cels <- list.files(paste(nam,"data",sep="/"), pattern = "[gz]")
sapply(paste(paste(nam,"data",sep="/"), cels, sep="/"), gunzip)
cels


setwd(paste(getwd(),nam,sep="/"))

library(simpleaffy)
celfiles <- read.affy(covdesc="phenodataCDN.txt", path="data")
celfiles.rma <- rma(celfiles)




require(arrayQualityMetrics)
arrayQualityMetrics(celfiles, outdir = "QAraw_CDN", force = FALSE, do.logtransform = TRUE, intgroup = "dis", reporttitle = paste("Quality metrics report for", "QCraw_CDN"))

arrayQualityMetrics(celfiles.rma, outdir = "QARMA_CDN", force = FALSE, do.logtransform = TRUE, intgroup = "dis", reporttitle = paste("Quality metrics report for", "QCRMA_CDN"))




data_rma_linear <- 2**exprs(celfiles.rma)
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(row.names(data_rma_linear), "hgu133plus2")
write.table(cbind(gene.symbols,data_rma_linear),file="Linear_RMA.xls",sep='\t')


celfiles.filtered <- nsFilter(celfiles.rma, require.entrez=FALSE, remove.dupEntrez=T)
samples <- celfiles.rma$dis
samples <- as.factor(samples)

design <- model.matrix(~0 + samples)
colnames(design)<- levels(samples)

library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(CD_Nor=CD-Nor,levels=design)
contrast.matrix

eset2 <- exprs(celfiles.filtered$eset)

huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
write.table(huvec_ebFit,file="ALL_Tested.xls",sep='\t')


top <- topTable(huvec_ebFit, coef=1, n=nrow(eset2), adjust="fdr")
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="CD_Nor_FDR_all.xls",sep='\t')

top <- topTable(huvec_ebFit, coef=1, n=nrow(eset2), adjust="none")
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="CD_Nor_all.xls",sep='\t')




top <- topTable(huvec_ebFit, coef=1, n="null", adjust="fdr",p.value=0.05,lfc=.584)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="4183_1.5-FDR05.xls",sep='\t')


top <- topTable(huvec_ebFit, coef=1, n="null", adjust="none",p.value=0.05,lfc=.584)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="4183_1.5-05.xls",sep='\t')



top <- topTable(huvec_ebFit, coef=1, n="null", adjust="none",p.value=0.01,lfc=.584)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="4183_1.5-01.xls",sep='\t')



top <- topTable(huvec_ebFit, coef=1, n="null", adjust="none",p.value=0.001,lfc=.584)
gene.symbols <- getSYMBOL(row.names(top), "hgu133plus2")
results <- cbind(top, gene.symbols)
results <- data.frame(results,eset2[row.names(top),])
write.table(results,file="4183_1.5-001.xls",sep='\t')


