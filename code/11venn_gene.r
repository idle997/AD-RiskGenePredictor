?Sys.setenv
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
##install.packages("VennDiagram")


library(VennDiagram) 
diffFile="diff.txt"       #差异分析的结果文件
wgcnaFile1="module_red.txt"     #共表达的结果文件
wgcnaFile2="module_green.txt"
wgcnaFile3="module_greenyellow.txt"
wgcnaFile4="module_magenta.txt"
setwd("Your path")
geneList=list()

#整合wgcna分析结果
#读取模块基因的列表文件
rt1=read.table(wgcnaFile1, header=F, sep="\t", check.names=F)
rt2=read.table(wgcnaFile2, header=F, sep="\t", check.names=F)
rt3=read.table(wgcnaFile3, header=F, sep="\t", check.names=F)
rt4=read.table(wgcnaFile4, header=F, sep="\t", check.names=F)

# 将四个数据框的基因名列合并
rt <- c(rt1$V1, rt2$V1, rt3$V1,rt4$V1)
# 去除重复的基因名
unique_genes <- unique(rt)
# 将结果转换为数据框
rt <- data.frame(Gene = unique_genes)
geneNames=as.vector(rt[,1])              #提取模块基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)
geneList[["WGCNA"]]=uniqGene

#读取差异分析的结果文件
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取差异基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对差异基因取unique
geneList[["DEG"]]=uniqGene

#读取孟德尔随机化分析结果文件
rt=read.csv("IVW.filter.csv", header=T, sep=",", check.names=F)
rt=rt[rt$method=="Inverse variance weighted",]
upRT=rt[rt$or>1,]
upGenes=unique(upRT[,"exposure"])
geneList[["MR_or>1"]]=upGenes

#绘制venn图
library(VennDiagram)
library(grid)   # 绘制 grid 对象需要

# 绘制 Venn 图，返回一个 gList 对象
venn_plot <- venn.diagram(
  x = geneList,
  filename = NULL,   # 不直接输出文件
  fill = c("red", "blue", "green"),
  alpha = 0.3,
  cex = 2,
  cat.cex = 1.5
)

# 打开 PDF 设备
pdf(file = "vennplot5.pdf", width = 8, height = 8)

# 使用 grid.draw 绘制到 PDF
grid.draw(venn_plot)

# 关闭设备
dev.off()





#输出交集基因的列表
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)







