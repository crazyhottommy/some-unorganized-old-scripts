
setwd("/home/tommy/MCF7_hypoxia_datasets/RNA-seq")
library("DESeq")
countsTable<- read.delim("RNA-seq_counts_final.txt", header=TRUE)
rownames(countsTable)<- countsTable$Gene
countsTable<- countsTable[,-1]
countsTable<- countsTable[which(rowSums(countsTable)>30),]
head(countsTable)

# if no other factors, just use the conds vector when constructing the cds object
# conds<- factor(c("hypoixa","normoxia", "normoxia","hypoxia","hypoxia","hypoxia"))

abDesign<- data.frame(row.names=colnames(countsTable),
                      condition = c("hypoxia","normoxia","normoxia","hypoxia","hypoxia","hypoxia"),
                      libType = c("single-end","single-end","single-end","single-end","paired-end","paired-end"))

# read the DESeq mannual, using abDesign$condition makes DESeq not be aware of the multi-factor design 
# cds<- newCountDataSet(countsTable, abDesign$condition)
############### 
# make DESeq aware of the libType factor in addtion to the treatment factor, put abDesign dataframe into cds
cds<- newCountDataSet(countsTable, abDesign)
###############
cds<- estimateSizeFactors(cds)

sizeFactors(cds)

head(counts(cds))
head(counts(cds,normalized=TRUE))
cds<- estimateDispersions(cds)

plotDispEsts(cds)

############################# This is for multi-factor design comparasion 
###### fit full model and reduced model 


fit1<- fitNbinomGLMs(cds, count ~ condition + libType)
fit0<- fitNbinomGLMs(cds, count ~ libType )
# fit_null<- fitNbinomGLMs(cds, count ~1)  # null hypothesis
# or including interaction term
# fit_full<- fitNbinomGLMs(cds, count~ condition + libType + condition:libType)

fit1<- fitNbinomGLMs(cds, count ~ libType + condition)
fit0<- fitNbinomGLMs(cds, count ~ libType )


pvalsGLM<- nbinomGLMTest (fit1, fit0)
padjGLM<- p.adjust ( pvalsGLM, method= "BH")

df<- data.frame(geneID =row.names(counts(cds)), pval=pvalsGLM, padj=padjGLM, fit1)
df<- df[order(df$padj, decreasing =F),]

head(df)

# str(fit1)
# fit1[(padjGLM <=0.05) & !is.na(padjGLM),]
# plot(density(pvalsGLM))


head(df)


hist(df$pval, breaks=100, col='skyblue', border='slateblue', main='')
write.table(df, file="results.txt", sep="\t", quote=F)

# the colname with name conditionnormoxia is the log2 fold change (however it is normoxia/hypoxia ratio, one needs to change the minus sign to plus)


################################## for single factor design

res<- nbinomTest( cds, 'normoxia', 'hypoxia' )
head(res)
plotMA(res)

# identify()


hist(res$pval, breaks=100, col='skyblue', border='slateblue', main='')
head(res)
head(res[order(res$padj, -res$foldChange),], n=20)
resSig = res[ res$padj < 0.05 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]
resSig = res[ res$padj < 0.05 & (res$log2FoldChange >1), ]
head(resSig[order(resSig$padj, -resSig$foldChange),], n=20)
write.table(res, "hypoxia_genes_all_DESeq.txt", sep="\t", quote=F)


cdsBlind<- estimateDispersions(cds, method= "blind")
vsd<- varianceStabilizingTransformation(cdsBlind)
plotPCA(vsd, intgroup=c("condition","libType"))

library("RColorBrewer")
library("gplots")
library("genefilter")


rv<-rowVars(exprs(vsd))  #vsd is an expression set object after variance stablizing transformation
idx<- order(-rv)[1:500]  # idx by the variance


index<- res$padj<0.05    # significantly changed genes
hmcols<- colorRampPalette(c("green","green4","red","red4"))(256)

#set scale='row'  get a standardized z-score to remove the mean gene expression 

heatmap.2(exprs(vsd)[index,],col=hmcols, trace="none",margin=c(10,6), scale='row', density.info="none")

# the defaut distance calculating method in hclust is Euclidian, but it is problemetic. genes are upregulated
# did not cluster together. see a post here https://www.biostars.org/p/91978/ and 
# here http://liorpachter.wordpress.com/2014/01/19/why-do-you-look-at-the-speck-in-your-sisters-quilt-plot-and-pay-no-attention-to-the-plank-in-your-own-heat-map/
# let's try dist but using the method single
mat<- exprs(vsd)[index,]
dist<- dist(mat)
hc.rows<- hclust(dist, method='single')
plot(hc.rows)
rowDend<- as.dendrogram(hc.rows)
heatmap.2(mat[,c(2,3,1,4,5,6)], col=hmcols, Rowv=rowDend, dendrogram = "row", Colv=F, trace="none",margin=c(10,6), scale='row', density.info="none")

# or simplely
heatmap.2(mat[,c(2,3,1,4,5,6)], col=hmcols, hclust=function(x) hclust(x, method='single'), distfun=function(x) dist(x, method='euclidean'), dendrogram = "both", Colv=T, trace="none",margin=c(10,6), scale='row', density.info="none")


#####################################################
########### this is many times the best way to cluster using person correlation distance
# let's try to compute the pearson distances instead
# do something like this http://stackoverflow.com/questions/6719747/heatmap-of-microarray-data-using-pearson-distance
# both row and colum are clustered by pearson distance
heatmap.2(mat[,c(2,3,5,6)], col=hmcols, hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(10,6), scale='row', density.info="none")

#####################################################
dists<- dist(t(exprs(vsd)))
mat<- as.matrix(dists)
rownames(mat)=colnames(mat)= with(pData(cdsBlind), paste(condition))
heatmap.2(mat, trace="none", col=rev(hmcols), margin=c(13,13))


################################

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
