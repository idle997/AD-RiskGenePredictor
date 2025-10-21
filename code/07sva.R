

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)
setwd("Your path")      #设置工作目录

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
  if(file=="merge.preNorm.txt"){next}
  if(file=="merge.normalize.txt"){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
  geneNames=as.vector(rt[,1])      #提取基因名称
  uniqGene=unique(geneNames)       #基因取unique
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}
# 获取交集基因
interGenes = Reduce(intersect, geneList)

# 数据合并
allTab = data.frame()
batchType = c()

for (i in 1:length(files)) {
  inputFile = files[i]
  if (inputFile == "merge.preNorm.txt") { next }
  if (inputFile == "merge.normalize.txt") { next }
  
  header = unlist(strsplit(inputFile, "\\.|\\-"))
  
  # 读取输入文件，并对输入文件进行整理
  rt = read.table(inputFile, header = T, sep = "\t", check.names = F)
  rt = as.matrix(rt)
  rownames(rt) = rt[, 1]
  exp = rt[, 2:ncol(rt)]
  dimnames = list(rownames(exp), colnames(exp))
  data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  rt = avereps(data)
  
  # 确保基因名称一致
  rownames(rt) = toupper(rownames(rt))  # 统一为大写
  interGenes = toupper(interGenes)  # 统一为大写
  
  # 检查是否有缺失基因
  missingGenes = setdiff(interGenes, rownames(rt))
  if (length(missingGenes) > 0) {
    warning(paste("Missing genes in", inputFile, ":", paste(missingGenes, collapse = ", ")))
  }
  
  # 只保留交集基因
  rt = rt[intersect(interGenes, rownames(rt)), ]
  
  colnames(rt) = paste0(header[1], "_", colnames(rt))
  
  # 数据合并
  if (i == 1) {
    allTab = rt
  } else {
    allTab = cbind(allTab, rt[rownames(allTab), ])
  }
  
  batchType = c(batchType, rep(i, ncol(rt)))
}

#输出合并后的表达数据
outTab=rbind(geneNames=colnames(allTab), allTab)
write.table(outTab, file="merge.preNorm.txt", sep="\t", quote=F, col.names=F)

#对合并后数据进行批次矫正，输出批次矫正后的表达数据
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)

