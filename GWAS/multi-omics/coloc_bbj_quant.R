
library(data.table)

# 加载必要的包
library(foreach)
library(doParallel)
library("coloc")
library(dplyr)

# 设置并行处理的核心数
num_cores <- 30
registerDoParallel(cores = num_cores)


# 读取路径文件
path <- fread("mr_bbj220/coloc/1coloc_path.txt")
path <- filter(path, type == "Quantitative")

# 定义两个目录
eqtl_input <- path$path
metab_input <- path$file_path

######################
# 多核并行处理
results <- foreach(i = 1:nrow(path), .combine = 'rbind', .packages = c("data.table", "coloc", "dplyr")) %dopar% {
  tryCatch({
    # 读取文件路径
    eqtl_file <- eqtl_input[i]
    metab_file <- metab_input[i]
    
    # 读取文件
    eqtl <- fread(file = eqtl_file)
    metab <- fread(file = metab_file)
    
    # 从 path 文件获取的 ss 值和 outcome 名称
    nn <- as.numeric(path$N_samples[i])
    bbj_name <- path$outcome[i] 
    
    # 设置 eqtl 数据的列名
    colnames(eqtl)[6:8] <- c( "A1", "A2", "FRQ")
    colnames(eqtl)[15] <- "P"
    eqtl$CHR <- as.numeric(eqtl$CHR)
    
    # 合并 eqtl 和 metab 数据
    input <- merge(eqtl, metab, by = c("CHR", "BP"), all = FALSE, suffixes = c("_eqtl", "_metab"))
    
    if (nrow(input) > 0) {
      # 去除重复的 SNPs
      input <- input[!duplicated(input$SNP_eqtl), ]
      input <- input[!duplicated(input$SNP_metab), ]
      input$SNP_eqtl <- paste(input$SNP_eqtl, bbj_name, input$CHR, input$BP, "P_bbj", input$P_eqtl, "P_metab", input$P_metab, sep = ";")
      input$SNP_metab <- input$SNP_eqtl
      
      # 运行 coloc.abf 分析
      result <- coloc.abf(dataset1=list(beta = input$BETA_eqtl, varbeta = (input$SE_eqtl)^2, 
                                        pvalues = input$P_eqtl, N = nn, 
                                        MAF = input$FRQ,type="quant",
                                        snp = input$SNP_eqtl,position=input$BP),
                          dataset2=list(beta = input$BETA_metab, varbeta = (input$SE_metab)^2, 
                                        pvalues = input$P_metab, N = input$N, 
                                        MAF = input$MAF, type="quant",
                                        snp = input$SNP_metab,position=input$BP),p12=1e-5)
      
      # 提取 PP4
      pp4_value <- max(result$results$SNP.PP.H4)
      
      # 添加 EqtlFile 和 MetabFile 列到 ColocData
      colocdata <- result$results[which.max(result$results$SNP.PP.H4), ]
      colocdata$EqtlFile <- basename(eqtl_file)
      colocdata$MetabFile <- basename(metab_file)
      colocdata$noverlap <- nrow(result$results)
      sum1=result$summary
      # 将 sum1 转换为单行数据框
      sum1 <- as.data.frame(t(sum1))
      
      # 直接返回合并后的 data.frame
      return(cbind(
        sum1,
        data.frame(EqtlFile = basename(eqtl_file), MetabFile = basename(metab_file), PP4 = pp4_value),
        colocdata
      ))
      
    } else {
      # 空结果处理
      return(data.frame(
        nsnps=numeric(0),
        PP.H0.abf=numeric(0),
        PP.H1.abf=numeric(0),
        PP.H2.abf=numeric(0),
        PP.H3.abf=numeric(0),
        PP.H4.abf=numeric(0),
        EqtlFile = basename(eqtl_file),
        MetabFile = basename(metab_file),
        PP4 = 0,
        SNP = character(0),  # 空的 SNP 列
        position = integer(0),  # 空的 position 列
        V.df1 = numeric(0),
        z.df1 = numeric(0),
        r.df1 = numeric(0),
        lABF.df1 = numeric(0),
        V.df2 = numeric(0),
        z.df2 = numeric(0),
        r.df2 = numeric(0),
        lABF.df2 = numeric(0),
        internal.sum.lABF = numeric(0),
        SNP.PP.H4 = numeric(0),
        noverlap = integer(0)
      ))
    }
  }, error = function(e) {
    # 输出错误信息
    cat("Error in processing files:", eqtl_file, ";", metab_file, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    NULL
  })
}
# 停止并行后台
stopImplicitCluster()
# 输出结果
write.table(results, file = "mr_bbj220/coloc/result/summary_bbj_quan_group.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

