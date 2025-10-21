?Sys.setenv
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


##install.packages("tidyverse")
##BiocManager::install('GEOquery')
##install.packages('GEOquery')
##BiocManager::install("GEOquery")



library(Biobase)
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("Your path")
gset = getGEO(GEO='GSE132903', destdir=".",getGPL = F)


e2=gset[[1]]

exp2=e2@assayData[["exprs"]]

##提取数据临床信息
phe=e2@phenoData@data

library(data.table)
anno=fread("GPL10558-50081.txt",header =T,data.table = F)
a3=data.frame(anno$ID,anno$`Symbol`)
colnames(a3) <- c("ID", "gene.all")
a3$gene.all <- sub(" ///.*", "", a3$gene.all)


exp=as.data.frame(exp2)###讲矩阵转化为数据框的目的是为了好对表达矩阵进行操作

exp1=merge(x=a3,y=exp,by.x =1 ,by.y =0 )


rownames(exp1)=exp1$gene.all

exp2=distinct(exp1,gene.all,.keep_all = T)

rownames(exp2)=exp2$gene.all

exp3=na.omit(exp2)
rownames(exp3)=exp3$gene.all
exp4 <- exp3[!apply(exp3, 1, function(row) any(grepl("^\\s*$", row))), ]
exp4=exp4[,-c(1,2)]

exp5 <- cbind(geneNames = rownames(exp4), exp4)

write.table(exp5, file = "geneMatrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(phe, file = "clinical.csv", row.names = FALSE)




