library(gplots)

h1 <- read.table("Heatmap1.txt", header=T, row.names=1, sep="\t")
pdf("heatmap1.pdf")
heatmap.2(as.matrix(h1), dendrogram ='row', Colv=FALSE, col=greenred(60), 
                 key=FALSE,keysize=1.0, symkey=FLASE,density.info='none',
                 trace='none', colsep=1:11,
                 sepcolor='white', sepwidth=0.001,
                 scale="none",cexRow=0.5,cexCol=1,
                 labCol = colnames(h1), margins=c(5,25),
                 hclustfun=function(c){hclust(c, method='complete')},
                 lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.5, 21, 0.5 )
                 )
dev.off()

h2 <- read.table("Heatmap2.txt", header=T, row.names=1, sep="\t")
pdf("heatmap2.pdf")
heatmap.2(as.matrix(h2), dendrogram ='row', Colv=FALSE, col=greenred(60), 
                 key=FALSE, trace='none', colsep=1:10, rowsep=1:100,
                 sepcolor='white', sepwidth=0.001,
                 scale="none",cexRow=0.6,cexCol=1,
                 hclustfun=function(c){hclust(c, method='complete')},
		margins = c(11, 5),
                offsetRow = 0.5,
                offsetCol = 0.5, lhei=c(2,20),
                 )
dev.off()
