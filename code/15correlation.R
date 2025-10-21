#install.packages("corrplot")
#install.packages("circlize")
# 如果没有安装这些包，请先安装它们

# 加载所需的R包
library(corrplot)  # 用于绘制相关性图
library(circlize)  # 用于绘制环形图
library(limma)


expFile = "merge.normalize.txt"
geneFile = "interGenes.txt"
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


# 选择"Treat"组的数据
group = gsub("(.*)\\_(.*)", "\\2", colnames(geneExp))  # 提取样本分组信息
geneExp = geneExp[, group == "Treat", drop = FALSE]  # 选择分组为"Treat"的列数据
rt = t(geneExp)  # 转置数据矩阵，使基因为列，样本为行

# 计算相关性矩阵
cor1 = cor(rt)  # 计算基因之间的相关性矩阵

# 设置相关性矩阵的颜色
col = c(rgb(1, 0, 0, seq(1, 0, length = 32)), rgb(0, 1, 0, seq(0, 1, length = 32)))  # 颜色渐变从红色到绿色
cor1[cor1 == 1] = 0  # 将完全相关（值为1）的部分设为0，以便在绘图时不显示
c1 = ifelse(c(cor1) >= 0, rgb(1, 0, 0, abs(cor1)), rgb(0, 1, 0, abs(cor1)))  # 根据相关性值设置颜色
col1 = matrix(c1, nc = ncol(rt))  # 将颜色矩阵按数据矩阵的列数重新排列

# 绘制环形相关性图
pdf(file = "circos.pdf", width = 7, height = 7)  # 创建PDF文件，指定宽度和高度
par(mar = c(2, 2, 2, 4))  # 设置边距
circos.par(gap.degree = c(3, rep(2, nrow(cor1) - 1)), start.degree = 180)  # 设置环形图的参数
chordDiagram(cor1, grid.col = rainbow(ncol(rt)), col = col1, transparency = 0.5, symmetric = TRUE)  # 绘制环形图
par(xpd = TRUE)  # 允许图形超出绘图区域
colorlegend(col, vertical = TRUE, labels = c(1, 0, -1), xlim = c(1.1, 1.3), ylim = c(-0.4, 0.4))  # 添加颜色图例
dev.off()  # 关闭PDF设备
circos.clear()  # 清除circos设置

# 绘制相关性矩阵的热图
pdf(file = "corrplot.pdf", width = 7, height = 7)  # 创建PDF文件，指定宽度和高度
corrplot(cor1,  # 绘制相关性图
         method = "pie",  # 使用饼图表示相关性
         order = "hclust",  # 按照层次聚类的顺序排列
         type = "upper",  # 只显示上三角矩阵
         col = colorRampPalette(c("green", "white", "red"))(50)  # 颜色渐变从绿色到白色再到红色
)
dev.off()  # 关闭PDF设备
