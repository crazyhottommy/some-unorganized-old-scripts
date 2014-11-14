library(gplots)
getwd()
setwd("/home/tommy/Tet1/shDnmt3L")
d<- read.table("co_up_or_down_uniq.txt", header=T)

# heatmap.2 works only with matrix, convert the dataframe to matrix
m<-as.matrix(d[,2:3])
rownames(m)<- d$genes # add the gene names as the row lable

png(filename = "co_regulated1.png", width=400, height = 800) #save the heatmap to a png or a pdf by pdf(filename=...)

bk = unique(c(seq(-5.73,-0.5, length=100),seq(-0.5,0, length=100), seq(0,5,length=100)))
hmcols<- colorRampPalette(c("green","black", "red"))(length(bk)-1)

heatmap.2(m, breaks=bk, col=hmcols,trace="none",Colv=FALSE, dendrogram = "row",density.info="none",cexCol=1.2, labRow=NA, symm=F,symkey=F,symbreaks=T, scale="none",keysize=1.5)
dev.off()

# to get the the matrix after clustering
hm<- heatmap.2(m, breaks=bk, col=hmcols,trace="none",Colv=FALSE, dendrogram = "row",density.info="none",cexCol=1.2, labRow=NA, symm=F,symkey=F,symbreaks=T, scale="none")

names(hm)
# return the maxtrix returned after clustering as in the heatmap
m.afterclust<- m[rev(hm$rowInd),rev(hm$colInd)]

# to extract subgroups that are clustered together
# rowDendrogram is a list object 
# convert the rowDendrogram to a hclust object
hc<- as.hclust(hm$rowDendrogram)

names(hc)
plot(hc)  # rotate the dendrogram 90 degree, it is the same as in the heatmap
rect.hclust(hc,h=8) # based on the height of the tree, you can specify h

ct<- cutree(hc,h=8)

# get the members of each subgroup in the order of the cluster(left--->right), the row order will
# it is reversed compared to the heatmap.
table(ct)
ct[hc$order]


# get the matrix after clustering in the order of the heatmap (up--->down)

tableclustn<-  data.frame(m.afterclust, rev(ct[hc$order]))
head(tableclustn)
write.table(tableclustn, file="tableclustn.xls", row.names=T, sep="\t")

# remake the heatmap adding the RowSide bar based on the subgroups

png("Rheatmap2.png", width=400, height=800)
mycolhc<- sample(rainbow(256))
mycolhc<-mycolhc[as.vector(ct)]
rowDend<- as.dendrogram(hc)

heatmap.2(m, breaks=bk, Rowv=rowDend, Colv = FALSE, dendrogram = "row", col=hmcols, RowSideColors=mycolhc, labRow=NA, cexCol=1.2,trace="none", density.info="none", keysize=1.5)

dev.off()
