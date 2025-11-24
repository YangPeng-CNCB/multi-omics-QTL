

rm(list=ls())

library(mediation)
library(data.table)
library(bda)
##蛋白和代谢

#读基因型
geno=fread("中介分析/input/readable_geno.raw")
colnames(geno) <- gsub("_.*", "", colnames(geno))
geno=geno[,-3:-6]

#读输入数据
pro=fread("中介分析/input/pro_input.txt")
metab=fread("中介分析/input/metab_input.txt")
me=fread("中介分析/input/me_input.txt")

#读取要做中介分析的pair,并且对齐snp
map=fread("中介分析/pair_infomation/readable_geno.map")
map=map[,c(1,2,4)]
colnames(map)=c("CHR","snp","BP")
pairs=fread("中介分析/pair_infomation/all_pairs_new.txt")

pairs=merge(map,pairs,by=c("CHR","BP"),suffixes = c("","_pair"))

#合并输入数据
input1=merge(pro,metab,by=c("FID","IID"))
input2=merge(input1,me,by=c("FID","IID"))
input=merge(input2,geno,by=c("FID","IID"))
input=as.data.frame(input)


result=data.frame()


for (i in 1:nrow(pairs)) {
  tryCatch({
  # 获取每个pair的snp, m_id, CpG
  snp <- as.character(pairs$snp[i])
  m_id <- as.character(pairs$m_id[i])
  cpg <- as.character(pairs$CpG[i])
  tag=paste(snp,m_id,cpg,sep = ";")
  expr <- input[, m_id]  
  methy <- input[, cpg]  
  geno <- input[, snp]  
  
  valid_data <- complete.cases(expr, methy, geno)  
  expr <- expr[valid_data]
  methy <- methy[valid_data]
  geno <- geno[valid_data]
  
  # 第一步：拟合暴露（SNP）对中介的回归
  mod_mediator <- lm(methy ~ geno)  
  
  # 第二步：拟合暴露和中介对结果的联合影响的回归
  mod_outcome <- lm(expr ~ geno + methy)  
  
  #
  contcont <- mediate(mod_mediator, mod_outcome, sims=1000, treat="geno", mediator="methy")

  sobel_p = mediation.test(methy,geno,expr)[2,1]
  
  mediation_each_SME = summary(contcont)
  aa=data.frame(
    exposure = snp,
    mediator = cpg,
    outcome = m_id,
    ACME_Estimate = mediation_each_SME$d1,
    ACME_pvalue = mediation_each_SME$d1.p,
    ADE_Estimate = mediation_each_SME$z1,
    ADE_pvalue = mediation_each_SME$z1.p,
    Total_Effect_Estimate = mediation_each_SME$tau.coef,
    Total_Effect_pvalue = mediation_each_SME$tau.p,
    Proportion_Mediated_Estimate = mediation_each_SME$n1,
    Proportion_Mediated_pvalue = mediation_each_SME$n1.p,
    sobel_p=sobel_p,
    tag=tag
  )
  result=rbind(result,aa)
  
  # 第一步：拟合暴露（SNP）对中介的回归
  mod_mediator <- lm(expr ~ geno)  
  
  # 第二步：拟合暴露和中介对结果的联合影响的回归
  mod_outcome <- lm(methy ~ geno + expr) 
  
  #
  contcont <- mediate(mod_mediator, mod_outcome, sims=1000, treat="geno", mediator="expr")

  sobel_p = mediation.test(expr,geno,methy)[2,1]
  
  mediation_each_SME = summary(contcont)
  aa=data.frame(
    exposure = snp,
    mediator = m_id,
    outcome = cpg,
    ACME_Estimate = mediation_each_SME$d1,
    ACME_pvalue = mediation_each_SME$d1.p,
    ADE_Estimate = mediation_each_SME$z1,
    ADE_pvalue = mediation_each_SME$z1.p,
    Total_Effect_Estimate = mediation_each_SME$tau.coef,
    Total_Effect_pvalue = mediation_each_SME$tau.p,
    Proportion_Mediated_Estimate = mediation_each_SME$n1,
    Proportion_Mediated_pvalue = mediation_each_SME$n1.p,
    sobel_p=sobel_p,
    tag=tag
  )
  result=rbind(result,aa)
  }, error = function(e) {
    message(e$message)
    return(NULL)  
  })

}

write.table(result,file = "中介分析/result/all_result.txt",
            append = FALSE,quote = FALSE ,sep = "\t",col.names = T, row.names = FALSE)


















