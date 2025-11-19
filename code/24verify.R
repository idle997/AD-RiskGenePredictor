library(randomForest)
library(pROC)

# 读取训练好的 RF 模型（来自训练阶段）
setwd("Your path")
rf_model <- readRDS("RF_final_model.rds")   # 你的训练模型
features <- read.table("RF_features.txt", header=FALSE, stringsAsFactors=FALSE)[,1]  
features <- gsub("-", "_", features)  # 特征名标准化（与你训练代码一致）

# 读取外部表达数据
expFile="GSE28146.normalize.txt"
rt = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# 构造标签 y
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y == "Control", 0, 1)

# 仅保留模型训练时的特征
#    --- 外部表达矩阵必须匹配模型特征顺序 ---
rt = rt[features, ]    # 保留训练时使用的基因
rt = as.data.frame(t(rt))   # 行：样本，列：基因

# 使用训练时的模型进行预测
pred = predict(rf_model, newdata = rt, type = "prob")[,2]  # AD 的概率

#ROC 曲线
roc1 = roc(y, pred)
ci1 = ci.auc(roc1, method="bootstrap")

pdf(file="ROC_RF_External_GSE28146.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=TRUE,
     main="External Validation: Random Forest")
polygon(c(roc1$specificities, 1),
        c(roc1$sensitivities, 0),
        col=rgb(1,0,0,0.2), border=NA)

text(0.38, 0.43,
     paste0("95% CI: ", sprintf("%.3f", ci1[1]), "-", sprintf("%.3f", ci1[3])),
     col="red")
dev.off()
