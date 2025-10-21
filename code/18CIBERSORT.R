?Sys.setenv
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      #表达数据文件
refFile="humanREF.txt"         #参考文件
setwd("Your path")      #设置工作目录
source("CIBERSORT.R")       #引用包

#免疫细胞浸润分析
outTab=CIBERSORT(refFile, inputFile, perm=1000)

#对免疫浸润结果过滤，并且保存免疫细胞浸润的结果
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


