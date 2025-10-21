

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="GSE138260.normalize.txt"      #表达数据文件
geneFile="gene.txt"        #基因列表文件
setwd("Your path")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=t(data[as.vector(geneRT[,1]),])

#获取样品的分组信息(对照组和实验组)
Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

#将数据转换为箱线图的输入文件
data <- reshape2::melt(rt, id.vars = c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", fill = "Type",
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("#0088FF", "#FF5555"), width=0.75)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")
p1 = p1 + ggtitle("GSE138260") +
  theme(plot.title = element_text(hjust = 0.5))

#输出图形
pdf(file="boxplotGSE138260.pdf", width=6, height=4.5)
print(p1)
dev.off()


