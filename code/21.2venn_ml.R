?Sys.setenv
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
##install.packages("VennDiagram")
library(VennDiagram)
library(grid)
setwd("Your path")
rt1="importanceGene.RF.txt"
rt2="importanceGene.SVM.txt"
rt3="importanceGene.XGB.txt"
rt4="importanceGene.LASSO.txt"
rt5="importanceGene.GLM.txt"
geneList=list()

rt1=read.table(rt1, header=F, sep="\t", check.names=F)
rt2=read.table(rt2, header=F, sep="\t", check.names=F)
rt3=read.table(rt3, header=F, sep="\t", check.names=F)
rt4=read.table(rt4, header=F, sep="\t", check.names=F)
rt5=read.table(rt5, header=F, sep="\t", check.names=F)


rt1 <- rt1[-1, ]
rt1 <- c(rt1$V1)
rt1 <- data.frame(Gene = rt1)
geneNames=as.vector(rt1[,1])  
geneNames=gsub("^ | $","",geneNames) 
uniqGene=unique(geneNames)
geneList[["RF"]]=uniqGene

rt2 <- rt2[-1, ]
rt2 <- c(rt2$V1)
rt2 <- data.frame(Gene = rt2)
geneNames=as.vector(rt2[,1])  
geneNames=gsub("^ | $","",geneNames) 
uniqGene=unique(geneNames)
geneList[["SVM"]]=uniqGene

rt3 <- rt3[-1, ]
rt3 <- c(rt3$V1)
rt3 <- data.frame(Gene = rt3)
geneNames=as.vector(rt3[,1])  
geneNames=gsub("^ | $","",geneNames) 
uniqGene=unique(geneNames)
geneList[["XGB"]]=uniqGene

rt4 <- rt4[-1, ]
rt4 <- c(rt4$V1)
rt4 <- data.frame(Gene = rt4)
geneNames=as.vector(rt4[,1])  
geneNames=gsub("^ | $","",geneNames) 
uniqGene=unique(geneNames)
geneList[["LASSO"]]=uniqGene

rt5 <- rt5[-1, ]
rt5 <- c(rt5$V1)
rt5 <- data.frame(Gene = rt5)
geneNames=as.vector(rt5[,1])  
geneNames=gsub("^ | $","",geneNames) 
uniqGene=unique(geneNames)
geneList[["GLM"]]=uniqGene

venn_plot <- venn.diagram(
  x <- geneList,
  scaled = F, 
  filename = NULL,
  fill=c("#ff8181","#8181e5","#a5c5a2","blue","yellow"),
  alpha=0.3,
  height = 4000,
  width = 4000,
  cex = 2,
  cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), 
  cat.pos = c(0, -10, 240, 120, 20),
  cat.cex = 1.5)

pdf(file = "vennplot.pdf", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()

interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)
