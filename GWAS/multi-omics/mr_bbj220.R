
rm(list=ls())

# 目标目录
directory <- "mr_bbj220/metab_input/"
#directory <- "mr_bbj220/protein_input/"
# 获取目录下所有文件路径
all_files <- list.files(directory, recursive = TRUE, full.names = TRUE)

# 加载必要的包
library(foreach)
library(doParallel)

# 设置并行处理的核心数，这里设置为4，你可以根据需要调整
num_cores <- 30
registerDoParallel(cores = num_cores)


# 循环读入文件并处理
results_list <- foreach(file_path = all_files, .packages = c("TwoSampleMR")) %dopar% {
  # 读取数据
  exp_dat <- try(read_exposure_data(
    filename = file_path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "MAF",
    pval_col = "P",
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "BP",
    phenotype_col = "exposure"
  ), silent = TRUE)
  
  # 尝试读取 outcome 数据，如果失败返回空的数据框
  outcome_dat <- try(read_outcome_data(
    snps = exp_dat$SNP,
    filename = file_path,
    sep = "\t",
    snp_col = "SNP_out",
    beta_col = "BETA_out",
    se_col = "SE_out",
    effect_allele_col = "A1_out",
    other_allele_col = "A2_out",
    eaf_col = "FRQ",
    pval_col = "P_out",
    samplesize_col = "N_out",
    phenotype_col = "outcome"
  ), silent = TRUE)
  
  # 检查是否成功读取 outcome 数据
  if (inherits(outcome_dat, "try-error")) {
    cat("#######Error in processing", basename(file_path), "########\n")
    return(NULL)
  }
  
  # 检查 outcome 数据是否为空
  if (nrow(outcome_dat) == 0) {
    cat("No data found in", basename(file_path), "\n")
    return(NULL)
  }
  outcome_dat$outcome <- sub("^.*/([^/]+)/[^/]+$", "\\1", file_path)
  # 数据协调和 MR 分析
  dat <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = outcome_dat
  )
  dat=dat[dat$pval.outcome>5e-8,]
    # 检查 outcome 数据是否为空
  if (nrow(dat) == 0) {
    cat("No merge data found in", basename(file_path), "\n")
    return(NULL)
  }
  
  # 进行 MR 分析，获取结果
  results_mr <- mr(dat)
  
  # 进行异质性检验，获取结果
  heterogeneity <- mr_heterogeneity(dat)
  
  # 进行 Pleiotropy 测试，获取结果
  pleiotropy <- mr_pleiotropy_test(dat)
  
  # 返回三个数据框
  return(list(results_mr = results_mr, heterogeneity = heterogeneity, pleiotropy = pleiotropy, dat = dat))
}

# 关闭并行计算
stopImplicitCluster()


#####合并结果####
library(purrr)
library(dplyr)
library(data.table)
# 提取每个结果列表中的results_mr
results_mr_list <- map(results_list, pluck, "results_mr")
heterogeneity_list <- map(results_list, pluck, "heterogeneity")
pleiotropy_list <- map(results_list, pluck, "pleiotropy")
dat_list <- map(results_list, pluck, "dat")
# 合并所有的results_mr
combined_results_mr <- bind_rows(results_mr_list)
combined_heterogeneity <- bind_rows(heterogeneity_list)
combined_pleiotropy <- bind_rows(pleiotropy_list)
combined_dat <- bind_rows(dat_list)

# 输出最终结果
write.table(combined_results_mr, file = "mr_bbj220/metab_result/mrresult_results_all.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(combined_heterogeneity, file = "mr_bbj220/metab_result/mrresult_heterogeneity_all.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(combined_pleiotropy, file = "mr_bbj220/metab_result/mrresult_pleiotropy_all.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(combined_dat, file = "mr_bbj220/metab_result/mrresult_dat_all.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

####将结果合并####

a=combined_results_mr
b=combined_heterogeneity
c=combined_pleiotropy
colnames(c)[6]="se_pleiotropy"
colnames(c)[7]="p_pleiotropy"

ab=merge(a,b,by = c("id.exposure","id.outcome","outcome","exposure","method"),all.x = T)
abc=merge(ab,c,by = c("id.exposure","id.outcome","outcome","exposure"),all.x = T)
write.table(abc,file = "mr_bbj220/metab_result/mrresult_合并.txt",
            append = FALSE,quote = FALSE ,sep = "\t",col.names = T, row.names = FALSE)

abc <- abc[abc$method %in% c("Wald ratio", "Inverse variance weighted"), ]


write.table(abc,file = "mr_bbj220/metab_result/mrresult_合并_两种方法.txt",
            append = FALSE,quote = FALSE ,sep = "\t",col.names = T, row.names = FALSE)













