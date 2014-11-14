# VennDiagram
library("VennDiagram")

venn.plot <- draw.pairwise.venn(area1=42381, area2=3699, cross.area= 866, category=c("GATA3 binding sites","LSD1 binding sites"), fill=c("red","blue"), cex=1, cat.cex=1, cat.prompts=TRUE, cat.dist=c(-0.06,-0.06))
grid.draw(venn.plot)
