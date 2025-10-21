
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")


#引用包
library(TwoSampleMR)

inputFile="exposureID.txt"      #暴露数据id文件
setwd("Your path")     #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

#对暴露数据ID进行循环
outTab=data.frame()
for(id in rt$ID) {
	expoData=extract_instruments(id,
                                 p1 = 5e-08, p2 = 5e-08,
                                 clump = T,
                                 kb = 10000, r2 = 0.001)
    outTab=rbind(outTab, expoData)
}
write.csv(outTab, file="exposure_data.csv", row.names=F)



