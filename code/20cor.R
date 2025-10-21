# 安装必要的包（如果尚未安装）
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("reshape2")

# 加载必要的R包
library(limma)      # 用于数据预处理的包
library(reshape2)   # 用于数据重构的包
library(tidyverse)  # 数据清理与可视化工具包
library(ggplot2)    # 高级绘图工具包


# 设置工作目录
expFile = "merge.normalize.txt"
geneFile = "interGenes.txt"
immFile = "CIBERSORT-Results.txt"  # 免疫细胞浸润结果文件
setwd("Your path")  # 设置工作目录


# 读取表达数据文件，并进行数据处理
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)  # 读取带有标题的表达数据文件，禁用列名检查
rt = as.matrix(rt)  # 将数据框转换为矩阵形式
rownames(rt) = rt[,1]  # 将第一列设置为行名（基因名）
exp = rt[, 2:ncol(rt)]  # 去除第一列，只保留表达值数据
dimnames = list(rownames(exp), colnames(exp))  # 设置矩阵的行名和列名
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  # 将数据转换为数值矩阵，确保数据为数值类型
data = avereps(data)  # 使用avereps函数去除重复基因，合并重复的基因
data = data[rowMeans(data) > 0, ]  # 过滤掉平均表达量为0或以下的基因

# 读取基因列表文件，并筛选出在表达数据中同时存在的基因
gene = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)  # 读取基因列表文件
sameGene = intersect(as.vector(gene[, 1]), rownames(data))  # 找到基因列表与表达数据中的共同基因
geneExp = data[sameGene, ]  # 提取共同基因的表达数据



dimnames = list(rownames(geneExp), colnames(geneExp))  # 创建行名和列名的维度名称
data = matrix(as.numeric(as.matrix(geneExp)), nrow = nrow(geneExp), dimnames = dimnames)  # 将数据转换为数值矩阵
data = avereps(data)  # 对重复数据进行平均

# 筛选处理组样本数据
group = gsub("(.*)\\_(.*)", "\\2", colnames(data))  # 提取样本组别信息
data = data[, group == "Treat", drop = FALSE]  # 仅保留处理组的样本
data = t(data)  # 转置矩阵，使样本作为行，基因作为列

# 读取免疫细胞浸润数据文件，并对样本进行匹配
immune = read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
sameSample = intersect(row.names(data), row.names(immune))  # 获取两个数据集中共有的样本
data = data[sameSample, , drop = FALSE]  # 保留相同样本的基因表达数据
immune = immune[sameSample, , drop = FALSE]  # 保留相同样本的免疫细胞数据

# 计算每个基因与每个免疫细胞类型的相关性
outTab = data.frame()  # 初始化空数据框
for (cell in colnames(immune)) {  # 遍历免疫细胞类型
  if (sd(immune[, cell]) == 0) { next }  # 如果免疫细胞数据的标准差为0，跳过该细胞
  for (gene in colnames(data)) {  # 遍历每个基因
    x = as.numeric(immune[, cell])  # 提取免疫细胞数据
    y = as.numeric(data[, gene])  # 提取基因表达数据
    corT = cor.test(x, y, method = "spearman")  # 进行Spearman相关性检验
    cor = corT$estimate  # 提取相关系数
    pvalue = corT$p.value  # 提取p值
    text = ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))  # 根据p值标注显著性符号
    outTab = rbind(outTab, cbind(Gene = gene, Immune = cell, cor, text, pvalue))  # 将结果保存到outTab中
  }
}

# 生成相关性热图
outTab$cor = as.numeric(outTab$cor)  # 将相关系数转换为数值型
pdf(file = "cor.pdf", width = 7, height = 5)  # 保存为PDF
ggplot(outTab, aes(Immune, Gene)) +  # 绘制免疫细胞与基因的相关性图
  geom_tile(aes(fill = cor), colour = "grey", size = 1) +  # 使用颜色表示相关系数
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +  # 设置颜色梯度（负相关为蓝色，正相关为红色）
  geom_text(aes(label = text), col = "black", size = 3) +  # 添加显著性符号
  theme_minimal() +  # 使用简洁主题
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(),  # 隐藏轴标题和刻度
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),  # 设置x轴标签格式
        axis.text.y = element_text(size = 8, face = "bold")) +  # 设置y轴标签格式
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", " *  p<0.05", "\n", "\n", "Correlation")) +  # 添加图例标签
  scale_x_discrete(position = "bottom")  # 将x轴标签显示在下方
dev.off()  # 关闭PDF输出
