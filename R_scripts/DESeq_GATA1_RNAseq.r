setwd("/home/tommy/scripts")
library("DESeq")
countsTable<- read.delim("qiu_lab_with_header.txt", header=TRUE)
rownames(countsTable)<- countsTable$Gene
countsTable<- countsTable[,-1]
countsTable<- countsTable[which(rowSums(countsTable)>25),]
head(countsTable)
conds<- factor(c("GIE","GIE", "GATA1","V205M", "V205M","2RA"))

abDesign<- data.frame(row.names=colnames(countsTable),
                      condition = c("GIE","GIE","GATA1","V205M","V205M","2RA"),
                      libType = c("pair-end","pair-end","pair-end","pair-end","pair-end","pair-end"))


cds<- newCountDataSet(countsTable, abDesign$condition)
cds<- estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds))
head(counts(cds,normalized=TRUE))
cds<- estimateDispersions(cds)

plotDispEsts(cds)

res <- nbinomTest( cds, 'GIE', 'GATA1' )
head(res)
plotMA(res)
# if I use the command above, the red colored points (padj<0.1) are not at the edge of the cloud
# rather they spread around.

# there is a bug in the plotMA function that messes up with the order of the dataframe.
# see a post here https://www.biostars.org/p/103855/#106636
# my question post https://www.biostars.org/p/106588/#106627
# to fix it now:
plotMA(res[order(res$padj),]) 

# mannually creat the MA plot
# x=subset(res, res$baseMean!=0)
# col=ifelse(x$padj>=0.1, "gray32", "red3")
# plot(x$baseMean, x$log2FoldChange, col=col, log="x", ylim=c(-5,5), pch=19, cex=0.5)
# abline(h=0, col='red')

hist(res$padj, breaks=100, col='skyblue', border='slateblue', main='')

resSig <- res[ res$padj < 0.01 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]
resSig <- na.omit(resSig)
head(resSig)
write.table(resSig, "2RA_VS_V205M_2fold.txt", sep="\t", quote=F)

##### GAGE pathway analysis ###########
require(gage)
data(kegg.gs)
#use all the gene fold change 
deseq.fc<- res$log2FoldChange
names(deseq.fc)<- res$id
sum(is.infinite(deseq.fc))  # there are some infinite numbers, if use DESeq2, no such problem.
deseq.fc[deseq.fc>10]=10
deseq.fc[deseq.fc<-10]=-10
exp.fc<- deseq.fc

#kegg.gsets works with 3000 KEGG speicies
data(korg)
head(korg[,1:3], n=20)


#let's get the annotation files for mouse and convert the gene set to gene symbol format
kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
lapply(kegg.gs[1:3],head)

# egSymb is only for human data, so eg2sym and sym2eg functions are only for human data.
#data(egSymb)
#kegg.gs.sym<- lapply(kegg.gs, eg2sym)
#lapply(kegg.gs.sym[1:3],head)

# to convert IDs among gene/transcript ID to Entrez GeneID or reverse, use eg2id and id2eg in the pathview package written by the same person.
library(pathview)
data(bods)
bods

gene.symbol.eg<- id2eg(ids=names(exp.fc), category='SYMBOL', org='Mm') # convert the gene symbol to Entrez Gene ID
head(gene.symbol.eg, n=100)
head(gene.symbol.eg[,2], n=10)

names(exp.fc)<- gene.symbol.eg[,2]

fc.kegg.p<- gage(exp.fc, gsets= kegg.gs, ref=NULL, samp=NULL)
sel<- fc.kegg.p$greater[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel)

sel.l<- fc.kegg.p$less[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel.l)

##### PCA analysis ###################
=======
res = nbinomTest( cds, 'GATA1', 'GIE' )
head(res)
plotMA(res)

hist(res$padj, breaks=100, col='skyblue', border='slateblue', main='')

resSig = res[ res$padj < 0.01 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]
head(resSig)
write.table(resSig, "2RA_VS_V205M_2fold.txt", sep="\t", quote=F)


cdsBlind<- estimateDispersions(cds, method= "blind")
vsd<- varianceStabilizingTransformation(cdsBlind)
plotPCA(vsd, intgroup=c("condition"))

######## heatmap #####################


library("RColorBrewer")
library("gplots")
library("genefilter")


##library size normalized counts heatmap
new_counts_table<- counts(cds, normalize=T)
new_counts_table<- log2(new_counts_table)
rv<- rowVars(new_counts_table)
idx<- order(-rv)[1:500]
hmcols<- colorRampPalette(c("green","green4","red","red4","yellow"))(256)
heatmap.2(new_counts_table[idx,],col=hmcols,trace="none",Colv=T,density.info="none")

###variance stablized counts heatmap
rv<-rowVars(exprs(vsd))  #vsd is an expression set object after variance stablizing transformation
idx<- order(-rv)[1:500]
hmcols<- colorRampPalette(c("green","green4","red","red4"))(256)

heatmap.2(exprs(vsd)[idx,], col=hmcols, hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(10,6), scale='row', density.info="none",labRow=NA)
# ?cor:If x and y are matrices then the covariances (or correlations) between the columns of x and the columns of y are computed.
# ?dist: This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
# so dist computes the distance between rows while cor computes the correlation between columns

#####sample distance  heatmap,  note the transposed the matrix #########
dists<- dist(t(exprs(vsd)))
mat<- as.matrix(dists)
rownames(mat)=colnames(mat)= with(pData(cdsBlind), paste(condition))
heatmap.2(mat, trace="none", col=rev(hmcols), margin=c(13,13))

#multi-dimension scaling  MDS plot
mds<- cmdscale(dists)
km<- kmeans(t(exprs(vsd)), centers=4)
plot(mds, col=km$cluster, pch=16)


##### make a vacano plot ##################

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col="orange"))
with(subset(res, padj <.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot

library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)> 4), textxy(log2FoldChange, -log10(pval), labs=id, cex=.8))
