setwd("/home/tommy/HIF1_ChIP-seq")
library(GenomicRanges)
library(genomation)
HREs<- readGeneric("total_HRE_sites.bed",header=TRUE, keep.all.metadata=TRUE, zero.based=TRUE)
head(HREs)
promoters<- readGeneric("3kb_up_and_downstream_TSS.hg19.bed", strand=6, meta.col=list(transcript=4,gene_name=5))

# bam.files<- list.files(full.names=T, pattern="bam$")
bam.files<-c("./HIF1_hg19.sorted.bam","./Mcf7Jund.sorted.bam", "./DnaseMcf7contrl.sorted.bam")

HREs<- resize(HREs, width=6000, fix="center")

################## single bam file #############################
#sm<-ScoreMatrixBin(bam.files, HREs, bin.num=600, type="bam")
#sm.scaled<-scaleScoreMatrix(sm, scalefun=function(x) log2(x+1)) # or just use sm<- log2(sm+1) to scale it
#heatMatrix(sm.scaled,xcoords=c(-3000,3000),col=c("white","red"))
################## k-means clustering #############################################
sml<- ScoreMatrixList(bam.files, HREs, bin.num=600, type="bam")
sml.scaled<- scaleScoreMatrixList(sml, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sml.scaled)<- c("HIF1","Jund","Dnase")
set.seed(1000)
multiHeatMatrix(sml.scaled,xcoords=c(-3000,3000), col=c("white","red"),kmeans=TRUE, k=3, winsorize=c(0,99))

##################divide the HREs into looping and non-looping clusters#########################################
HREs<- readGeneric("HRE_total.bed",header=TRUE, keep.all.metadata=TRUE, zero.based=TRUE)
head(HREs)
looping.HREs<- which(HREs$looping_event==1)
non_looping.HREs<- which(HREs$looping_event==0)
HREs<- resize(HREs, width=6000, fix="center")

# bam.files<- list.files(full.names=T, pattern="bam$")
bam.files<-c("./HIF1_hg19.sorted.bam","./Mcf7Jund.sorted.bam", "./DnaseMcf7contrl.sorted.bam")

sml<- ScoreMatrixList(bam.files, HREs, bin.num=600, type="bam")
sml.scaled<- scaleScoreMatrixList(sml, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sml.scaled)<- c("HIF1","Jund","Dnase")
set.seed(1000)

multiHeatMatrix(sml.scaled,xcoords=c(-3000,3000), group=list(looping=looping.HREs,non_looping=non_looping.HREs), order=TRUE, col=c("white","red"), winsorize=c(0,99))
################### further divide into four groups promoter-looping, distal-looping etc #######
HREs<- readGeneric("HRE_total.bed",header=TRUE, keep.all.metadata=TRUE, zero.based=TRUE)
head(HREs)
promoters<- readGeneric("3kb_up_and_downstream_TSS.hg19.bed", strand=6, meta.col=list(transcript=4,gene_name=5))
head(promoters)
looping_distal_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)==0)
looping_promoter_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)>0)
non_looping_distal_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)==0)
non_looping_promoter_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)>0)

HREs<- resize(HREs, width=10000, fix="center")

bam.files<-c("./HIF1_hg19.sorted.bam", "./HIF2_hg19.sorted.bam","SRR358668_ctrl_H3k4me2.sorted.bam","SRR816996_ctrl_H3k27ac.sorted.bam")

sml<- ScoreMatrixList(bam.files, HREs, bin.num=1000, type="bam", rpm= T)
sml.scaled<- scaleScoreMatrixList(sml, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sml.scaled)<- c("HIF1","HIF2","H3k4me2","H3k27ac")
set.seed(1000)

multiHeatMatrix(sml.scaled,xcoords=c(-3000,3000), group=list(looping_distal=looping_distal_HRE, looping_promoter=looping_promoter_HRE, non_looping_distal=non_looping_distal_HRE, non_looping_promoter=non_looping_promoter_HRE),order=TRUE, col=c("white","red"), winsorize=c(0,99))
plotMeta(sml.scaled, xcoords=c(-5000,5000),profile.names=names(sml.scaled), line.col=c("blue", "red"))
heatMeta(sml.scaled, xcoords=c(-3000,3000),xlab="bp around peaks")
################################################ metaplot subgroups of HREs#################

looping_distal_HRE<- readGeneric("looping_distal_HRE_tag2.bed",header=FALSE, keep.all.metadata=TRUE, zero.based=TRUE)
non_looping_distal_HRE<- readGeneric("non_looping_distal_HRE_tag2.bed",header=FALSE, keep.all.metadata=TRUE, zero.based=TRUE)


looping_distal_HRE<- resize(looping_distal_HRE, width=10000, fix="center")

non_looping_distal_HRE<- resize(non_looping_distal_HRE, width=10000, fix="center")

bam.file<- list.files(".", full.names=T, pattern="HIF1_hg19.sorted.bam$")

sm1<- ScoreMatrixBin(target=bam.file, windows=looping_distal_HRE, bin.num=1000, strand.aware=FALSE, type="bam", rpm=T)
sm2<- ScoreMatrixBin(target=bam.file, windows=non_looping_distal_HRE, bin.num=1000, strand.aware=FALSE, type="bam",rpm=T)
sm.list<- new("ScoreMatrixList", list(LD=sm1,ND=sm2)) # create a ScoreMatrixList object 

sm.list.scaled<- scaleScoreMatrixList(sm.list, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sm.list.scaled)<- names(sm.list)
pdf("p300_meta.pdf", width= 4, height= 4)
plotMeta(sm.list.scaled, xcoords=c(-5000,5000), profile.names=names(sm.list), line.col=c("red","blue"))
dev.off()

####


bam.file2<- "Mcf7Jund.sorted.bam"
sm3<- ScoreMatrixBin(target=bam.file2, windows=looping_distal_HRE, bin.num=1000, strand.aware=FALSE, type="bam")
sm4<- ScoreMatrixBin(target=bam.file2, windows=non_looping_distal_HRE, bin.num=1000, strand.aware=FALSE, type="bam")
sm.list2<- new("ScoreMatrixList", list(looping_distal=sm3,non_looping_distal=sm4)) # create a ScoreMatrixList object 

sm.list.scaled2<- scaleScoreMatrixList(sm.list2, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sm.list.scaled2)<- names(sm.list2)
plotMeta(sm.list.scaled2, xcoords=c(-5000,5000), profile.names=names(sm.list2))

#========================order ScoreMatrixList from high to low in subgroups according to one of the matrix (HIF1) in the matrice=======================

HREs<- readGeneric("HRE_total.bed",header=TRUE, keep.all.metadata=TRUE, zero.based=TRUE)
head(HREs)
promoters<- readGeneric("3kb_up_and_downstream_TSS.hg19.bed", strand=6, meta.col=list(transcript=4,gene_name=5))
head(promoters)
looping_distal_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)==0)
looping_promoter_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)>0)
non_looping_distal_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)==0)
non_looping_promoter_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)>0)

HREs<- resize(HREs, width=6000, fix="center")

bam.files<-c("./HIF1_hg19.sorted.bam", "./HIF2_hg19.sorted.bam","SRR358668_ctrl_H3k4me2.sorted.bam","SRR816996_ctrl_H3k27ac.sorted.bam")

sml<- ScoreMatrixList(bam.files, HREs, bin.num=600, type="bam", rpm=T)
sml.scaled<- scaleScoreMatrixList(sml, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sml.scaled)<- c("HIF1","HIF2","H3k4me2","H3k27ac")

HIF1.matrix<- sml.scaled$HIF1
HIF1_order<- order(rowSums(HIF1.matrix), decreasing=TRUE)
sml.scaled.ordered<- orderBy(sml.scaled, HIF1_order))

set.seed(1000)

multiHeatMatrix(sml.scaled.ordered ,xcoords=c(-3000,3000), group=list(looping_distal=looping_distal_HRE, looping_promoter=looping_promoter_HRE, non_looping_distal=non_looping_distal_HRE, non_looping_promoter=non_looping_promoter_HRE),order=TRUE, col=c("white","red"), winsorize=c(0,99))
plotMeta(sml.scaled, xcoords=c(-3000,3000),profile.names=names(sml.scaled))
heatMeta(sml.scaled, xcoords=c(-3000,3000),xlab="bp around peaks")

##=============k-means clustering in sub-groups============######

HREs<- readGeneric("HRE_total.bed",header=TRUE, keep.all.metadata=TRUE, zero.based=TRUE)
head(HREs)
promoters<- readGeneric("3kb_up_and_downstream_TSS.hg19.bed", strand=6, meta.col=list(transcript=4,gene_name=5))
head(promoters)
looping_distal_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)==0)
looping_promoter_HRE<- which(HREs$looping_event==1 & countOverlaps(HREs,promoters)>0)
non_looping_distal_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)==0)
non_looping_promoter_HRE<-which(HREs$looping_event==0 & countOverlaps(HREs,promoters)>0)

HREs<- resize(HREs, width=10000, fix="center")

bam.files<-c("./HIF1_hg19.sorted.bam", "./HIF2_hg19.sorted.bam","SRR358668_ctrl_H3k4me2.sorted.bam","SRR816996_ctrl_H3k27ac.sorted.bam")

sml<- ScoreMatrixList(bam.files, HREs, bin.num=1000, type="bam", rpm=T)
sml.scaled<- scaleScoreMatrixList(sml, strand.aware=FALSE, scalefun=function(x) log2(x+1))
names(sml.scaled)<- c("HIF1","HIF2","H3k4me2","H3k27ac")

# creat sub-scoreMatrixList
sub1<- new("ScoreMatrixList", lapply(sml.scaled, function(x) x[looping_distal_HRE,]))
sub2<- new("ScoreMatrixList", lapply(sml.scaled, function(x) x[looping_promoter_HRE,]))
sub3<- new("ScoreMatrixList", lapply(sml.scaled, function(x) x[non_looping_distal_HRE,]))
sub4<- new("ScoreMatrixList", lapply(sml.scaled, function(x) x[non_looping_promoter_HRE,]))

# k-means clustering for each sub-scoreMatrixList, multiheatMatrix function returns the order of rows
set.seed(1000)

sub1.k.means<- multiHeatMatrix(sub1,xcoords=c(-5000,5000), kmeans=T, k=3, order=TRUE, col=c("white","red"))
sub2.k.means<- multiHeatMatrix(sub2,xcoords=c(-5000,5000), kmeans=T, k=3, order=TRUE, col=c("white","red"))
sub3.k.means<- multiHeatMatrix(sub3,xcoords=c(-5000,5000), kmeans=T, k=3, order=TRUE, col=c("white","red"))
sub4.k.means<- multiHeatMatrix(sub4,xcoords=c(-5000,5000), kmeans=T, k=3, order=TRUE, col=c("white","red"))

# increment each subgroup kmeans index by a k, so when plot together, they are not from the same k-means clustering
sub2.k.means<- sub2.k.means +3
sub3.k.means<- sub3.k.means +6
sub4.k.means<- sub4.k.means +9

sub.k.means<- c(sub1.k.means,sub2.k.means,sub3.k.means,sub4.k.means)
png("k.means_subgroup1.png", width=1600, height=2000)
multiHeatMatrix(sml.scaled, xcoords=c(-5000,5000), group=factor(sub.k.means), order=T,col=c("white","red"), common.scale=T, cex.lab=2, cex.main=2, cex.axis=2, cex.legend=2, winsorize=c(0,99))
dev.off()
