inputFile="exposure_data.csv"      #输入文件
setwd("Your path")     #设置工作目录

#读取输入文件
dat=read.csv(inputFile, header=T, sep=",", check.names=F)

#计算F检验值
dat$R2<-(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)/(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)+2*dat$se.exposure*dat$se.exposure*dat$samplesize.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)))     #计算R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)     #计算F检验值

#根据F值>10对数据进行过滤, 删除弱工具变量
outTab=dat[as.numeric(dat$F)>10,]
write.csv(outTab, file="exposure.F.csv", row.names=F)