# 安装所需的 R 包
#install.packages("caret")
#install.packages("DALEX")
#install.packages("ggplot2")
#install.packages("randomForest")
#install.packages("kernlab")
#install.packages("pROC")
#install.packages("xgboost")

# 引用已安装的 R 包
library(caret)         # 用于创建预测模型
library(DALEX)         # 用于解释机器学习模型
library(ggplot2)       # 用于数据可视化
library(randomForest)  # 用于随机森林模型
library(kernlab)       # 用于支持向量机模型
library(xgboost)       # 用于梯度提升模型
library(pROC)          # 用于绘制 ROC 曲线

# 设置随机种子以确保结果可重复
set.seed(123)

# 定义输入文件路径
inputFile="merge.normalize.txt"      # 基因表达数据文件
geneFile="interGenes.txt"      # 基因列表文件

# 设置工作目录
setwd("Your path")

# 读取表达数据文件，行名为基因名称，列名为样本名称
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# 读取基因列表文件，提取核心基因的表达数据
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]     # 根据基因列表提取表达矩阵中的对应基因
row.names(data)=gsub("-", "_", row.names(data))  # 将基因名称中的 "-" 替换为 "_"

# 获取样本分组信息，并将其作为新列添加到数据框中
data=t(data)   # 转置数据矩阵，使样本为行，基因为列
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))  # 提取样本分组信息
data=as.data.frame(data)  # 转换为数据框格式
data$Type=group           # 添加样本分组信息作为数据框中的一列

# 将数据分为训练集和测试集，训练集占总样本的 70%
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]  # 训练集数据
test<-data[-inTrain,]  # 测试集数据

# 创建 RF 随机森林模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)  # 定义交叉验证方法
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

# 创建 SVM 机器学习模型
mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

# 创建 XGB 模型
mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

# 创建 GLM 模型
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)

# 创建 GBM 模型
mod_gbm=train(Type ~., data = train, method = "gbm", trControl=control)

# 创建 KNN 模型
mod_knn=train(Type ~., data = train, method = "knn", trControl=control)

# 创建 NNET 模型
mod_nnet=train(Type ~., data = train, method = "nnet", trControl=control)

# 创建 LASSO 模型
mod_lasso=train(Type ~., data = train, method = "glmnet", trControl=control)

# 创建 DT 决策树模型
mod_dt=train(Type ~., data = train, method = "rpart", trControl=control)

# 定义预测函数，用于获取模型预测的概率
p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)  # 将测试集的分类结果转换为二进制

# RF 随机森林模型预测结果
explainer_rf=explain(mod_rf, label = "RF",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)

# SVM 机器学习模型预测结果
explainer_svm=explain(mod_svm, label = "SVM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)

# XGB 模型预测结果
explainer_xgb=explain(mod_xgb, label = "XGB",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)

# GLM 模型预测结果
explainer_glm=explain(mod_glm, label = "GLM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)

# GBM 模型预测结果
explainer_gbm=explain(mod_gbm, label = "GBM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_gbm=model_performance(explainer_gbm)

# KNN 模型预测结果
explainer_knn=explain(mod_knn, label = "KNN",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_knn=model_performance(explainer_knn)

# NNET 模型预测结果
explainer_nnet=explain(mod_nnet, label = "NNET",
                       data = test, y = yTest,
                       predict_function = p_fun,
                       verbose = FALSE)
mp_nnet=model_performance(explainer_nnet)

# LASSO 模型预测结果
explainer_lasso=explain(mod_lasso, label = "LASSO",
                        data = test, y = yTest,
                        predict_function = p_fun,
                        verbose = FALSE)
mp_lasso=model_performance(explainer_lasso)

# DT 决策树模型预测结果
explainer_dt=explain(mod_dt, label = "DT",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_dt=model_performance(explainer_dt)

# 绘制残差反向累计分布图
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, mp_gbm, mp_knn, mp_nnet,mp_lasso)
print(p1)
dev.off()

# 绘制残差箱线图
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, mp_gbm, mp_knn , mp_nnet, mp_lasso, geom = "boxplot")
print(p2)
dev.off()

# 绘制 ROC 曲线并保存到 PDF 文件
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
pred5=predict(mod_gbm, newdata=test, type="prob")
pred6=predict(mod_knn, newdata=test, type="prob")
pred7=predict(mod_nnet, newdata=test, type="prob")
pred8=predict(mod_lasso, newdata=test, type="prob")
pred9=predict(mod_dt, newdata=test, type="prob")

roc1=roc(yTest, as.numeric(pred1[,2]))  # 计算 ROC 曲线
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
roc5=roc(yTest, as.numeric(pred5[,2]))
roc6=roc(yTest, as.numeric(pred6[,2]))
roc7=roc(yTest, as.numeric(pred7[,2]))
roc8=roc(yTest, as.numeric(pred8[,2]))
roc9=roc(yTest, as.numeric(pred9[,2]))

pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")  # 绘制 ROC 曲线
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
plot(roc5, print.auc=F, legacy.axes=T, main="", col="orange", add=T)
plot(roc6, print.auc=F, legacy.axes=T, main="", col="purple", add=T)
plot(roc7, print.auc=F, legacy.axes=T, main="", col="black", add=T)
plot(roc8, print.auc=F, legacy.axes=T, main="", col="red", add=T)

# 添加图例
legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('XGB: ',sprintf("%.03f",roc3$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc)),
         paste0('GBM: ',sprintf("%.03f",roc5$auc)),
         paste0('KNN: ',sprintf("%.03f",roc6$auc)),
         paste0('NNET: ',sprintf("%.03f",roc7$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc8$auc))),
       col=c("red","blue","green","yellow","orange","purple","black","blue","red"), lwd=2, bty = 'n')
dev.off()

# 计算各模型的基因重要性评分
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
importance_knn<-variable_importance(
  explainer_knn,
  loss_function = loss_root_mean_square
)
importance_gbm<-variable_importance(
  explainer_gbm,
  loss_function = loss_root_mean_square
)
importance_nnet<-variable_importance(
  explainer_nnet,
  loss_function = loss_root_mean_square
)
importance_lasso<-variable_importance(
  explainer_lasso,
  loss_function = loss_root_mean_square
)
importance_dt<-variable_importance(
  explainer_dt,
  loss_function = loss_root_mean_square
)

# 绘制基因重要性图形并保存为 PDF 文件
pdf(file="importance.pdf", width=7, height=14)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_gbm[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_knn[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_nnet[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_lasso[c(1,(ncol(data)-8):(ncol(data)+1)),],
     importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()

# 输出重要性评分最高的基因到文本文件中
geneNum=5  # 设置输出基因的数量
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_gbm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GBM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_knn[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.KNN.txt", sep="\t", quote=F, row.names=F)
write.table(importance_nnet[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.NNET.txt", sep="\t", quote=F, row.names=F)
write.table(importance_lasso[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.LASSO.txt", sep="\t", quote=F, row.names=F)


