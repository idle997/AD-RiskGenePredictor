
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(TwoSampleMR)

exposureFile="exposure.F.csv"        #暴露数据文件
outcomeFile="outcome.csv"            #结局数据文件
geneFile="interGenes.txt"       #交集基因列表文件
outcomeName="Alzheimer disease"      #结局中展示疾病的名字
setwd("Your path")     #设置工作目录

#读取暴露数据输入文件
rt=read.csv(exposureFile, header=T, sep=",", check.names=F)

#读取交集基因的列表文件
sigGene=read.table(geneFile, header=F, sep="\t", check.names=F)

#对基因进行循环
for(geneName in as.vector(sigGene[,1])){
	#提取这个基因的暴露数据
	i=gsub("\\%|\\/|\\:", "_", geneName)
	singleExposureFile=paste0(i, ".exposure.csv")
	exposure_set=rt[rt$exposure==geneName,]
	write.csv(exposure_set, file=singleExposureFile, row.names=F)
	#读取这个基因的暴露数据
	exposure_dat=read_exposure_data(filename=singleExposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)
	
	#读取结局数据
	outcome_data=read_outcome_data(snps=exposure_dat$SNP,
		             filename="outcome.csv", sep = ",",
		             snp_col = "SNP",
		             beta_col = "beta.outcome",
		             se_col = "se.outcome",
		             effect_allele_col = "effect_allele.outcome",
		             other_allele_col = "other_allele.outcome",
		             pval_col = "pval.outcome",
		             eaf_col = "eaf.outcome")
	
	#将暴露数据和结局数据合并
	outcome_data$outcome=outcomeName
	dat=harmonise_data(exposure_dat, outcome_data)
	
	#输出用于孟德尔随机化的工具变量
	outTab=dat[dat$mr_keep=="TRUE",]
	write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
	
	#MR-PRESSO异常值检测(偏倚的SNP)
	#presso=run_mr_presso(dat)
	#write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
	#write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
	
	#孟德尔随机化分析
	mrResult=mr(dat)
	
	#对结果进行OR值的计算
	mrTab=generate_odds_ratios(mrResult)
	write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)
	
	#异质性检验
	heterTab=mr_heterogeneity(dat)
	write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)
	
	#多效性检验
	pleioTab=mr_pleiotropy_test(dat)
	write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)
	
	#绘制散点图
	pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
	p1=mr_scatter_plot(mrResult, dat)
	print(p1)
	dev.off()
	
	#森林图
	res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
	pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
	p2=mr_forest_plot(res_single)
	print(p2)
	dev.off()
	
	#漏斗图
	pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
	p3=mr_funnel_plot(singlesnp_results = res_single)
	print(p3)
	dev.off()
	
	#留一法敏感性分析
	pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
	p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
	print(p4)
	dev.off()
}

