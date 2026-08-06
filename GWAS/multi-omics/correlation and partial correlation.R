

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
pairs=fread("中介分析/pair_infomation/all_pairs.txt")

pairs=merge(map,pairs,by=c("CHR","BP"),suffixes = c("","_pair"))

pairs=pairs[pairs$PP.H4.abf>0.75,]

#合并输入数据
input1=merge(pro,metab,by=c("FID","IID"))
input2=merge(input1,me,by=c("FID","IID"))
input=merge(input2,geno,by=c("FID","IID"))
input=as.data.frame(input)

#
result=data.frame()

# 
for (i in 1:nrow(pairs)) {
  tryCatch({
    aa=data.frame()

    snp <- as.character(pairs$snp[i])
    m_id <- as.character(pairs$m_id[i])
    cpg <- as.character(pairs$CpG[i])
    tag=paste(snp,m_id,cpg,sep = ";")
    expr <- input[, m_id]  
    methy <- input[, cpg] 
    geno <- input[, snp]  
    data=data.frame(expr,methy,geno)
    # Pearson 相关性分析
    cor_result <- cor.test(data$expr, data$methy, method = "pearson")
    pearson_corr <- cor_result$estimate  
    p_value <- cor_result$p.value        
    
    # 偏相关分析
    residuals_data <- data.frame(
      expr_resid = residuals(lm(data$expr ~ data$geno,na.action = na.exclude)),
      methy_resid = residuals(lm(data$methy ~ data$geno,na.action = na.exclude))
    )
    
    
    partial_cor_result <- cor.test(residuals_data$expr_resid, residuals_data$methy_resid, method = "pearson")
    partial_corr <- partial_cor_result$estimate
    partial_p_value <- partial_cor_result$p.value
    
    # 保存结果
    aa <- data.frame(
      m_id = m_id,
      cpg = cpg,
      corr = pearson_corr,
      corr.p = p_value,
      partial_corr = partial_corr,
      partial_corr.p = partial_p_value,
      tag=tag
    )
    
    result=rbind(result,aa)
  }, error = function(e) {
    message(e$message)
    return(NULL)  
  })
  
}

write.table(result,file = "中介分析/result/相关性和偏相关性.txt",
            append = FALSE,quote = FALSE ,sep = "\t",col.names = T, row.names = FALSE)








