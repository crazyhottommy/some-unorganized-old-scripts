# read GEO data sets from NCBI by GEOquery

setwd("/home/tommy/Tet1")# set the working directory
library(Biobase)
library(GEOquery)

# only set the GSEMatrix to FALSE can it be parsed for later use of function like
# Meta(gse) 
gse<- getGEO('GSE26830', GSEMatrix=FALSE,  destdir=".")

Meta(gse)
names(GSMList(gse))

# I want the expression matrix, so this time I set the GSEMatrix to TRUE.
gse<- getGEO('GSE26830', GSEMatrix=TRUE,  destdir=".") #shTet1 in mESC

class(gse)
length(gse)
gse

#gse is a list of length 2
# use double brackets to access the elements of the list
class(gse[1]]) # the ExpressionSet object 
str(gse[[1]]) # check the data structure of gse

show(pData(phenoData(gse[[1]]))[,c(1,6,8)])

eset<- gse[[1]]
eset
featureNames(eset)[1:10]
sampleNames(eset)

#anotation file, packages from Bioconductor
source("http://www.bioconductor.org/biocLite.R")
#biocLite("mouse4302.db")
library("mouse4302.db")
library(annotate)
ID <- featureNames(eset)
Symbol<- getSYMBOL(ID, "mouse4302.db")
fData(eset)<- data.frame(ID, Symbol)


#differentially expressed genes by limma
library(limma)
design<- model.matrix(~0+factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design)<- c("group1", "group2","group3")
fit<- lmFit(eset, design)
contrast.matrix<- makeContrasts(group2-group1,group3-group2,group3-group1, levels=design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)

topgenes<- topTable(fit2, coef=1,number=2000, p.value=0.21, lfc=0.6, adjust="BH")
write.table(topgenes,"topgenes_Fdr_0.21.txt", sep="\t")
