library(data.table)

# 加载必要的包
library(foreach)
library(doParallel)
library("coloc")
library(dplyr)

# 设置并行处理核心数
num_cores <- 30
registerDoParallel(cores = num_cores)

# 定义文件目录和文件列表
eqtl_input <- "eqtlgen/eqtl/final/"
metab_input <- "coloc/metabinput/"
eqtl_files <- list.files(eqtl_input, pattern = "\\.txt$", full.names = TRUE)
metab_files <- list.files(metab_input, pattern = "\\.txt$", full.names = TRUE)

# 并行处理
colocdata_list <- foreach(eqtl_file = eqtl_files, .combine = rbind, .packages = c("coloc", "data.table", "dplyr")) %:% 
  foreach(metab_file = metab_files, .combine = rbind) %dopar% {
    tryCatch({
      # 读取文件
      eqtl <- fread(file = eqtl_file)
      metab <- fread(file = metab_file)
      
      # 合并数据
      input <- merge(eqtl, metab, by = c("CHR", "BP"), all = FALSE, suffixes = c("_eqtl", "_metab"))
      
      # 如果没有数据重叠，跳过此循环
      if (nrow(input) == 0) return(NULL)
      
      # 去除重复的SNPs
      input <- input[!duplicated(input$SNP_eqtl), ]
      input <- input[!duplicated(input$SNP_metab), ]
      input$SNP_eqtl <- paste(input$SNP_eqtl, input$CHR, input$BP, "P_eqtl", input$P_eqtl, "P_metab", input$P_metab, sep = ";")
      input$SNP_metab <- input$SNP_eqtl

      # 运行Coloc分析
      result <- coloc.abf(dataset1 = list(beta = input$BETA_eqtl, varbeta = (input$SE_eqtl)^2, 
                                          pvalues = input$P_eqtl, N = input$N_eqtl, 
                                          MAF = input$MAF_eqtl, type = "quant", snp = input$SNP_eqtl, position = input$BP),
                          dataset2 = list(beta = input$BETA_metab, varbeta = (input$SE_metab)^2, 
                                          pvalues = input$P_metab, N = input$N_metab, 
                                          MAF = input$MAF_metab, type = "quant", snp = input$SNP_metab, position = input$BP),
                          p12 = 1e-5)
      
      # 筛选并添加文件信息
      colocdata <- result$results[which.max(result$results$SNP.PP.H4), ]
      colocdata$EqtlFile <- basename(eqtl_file)
      colocdata$MetabFile <- basename(metab_file)
      ##
      sum1=result$summary
      # 将 sum1 转换为单行数据框
      sum1 <- as.data.frame(t(sum1))
      if(sum1$PP.H4.abf>=0.75){
      # 直接返回合并后的 data.frame
      return(cbind(sum1,colocdata))
      }else{return(NULL)}
      
    }, error = function(e) {
      message("Error in file combination ", basename(eqtl_file), " and ", basename(metab_file), ": ", e$message)
      return(NULL)  # 出现错误时返回 NULL，以跳过该组合
    })
  }

# 关闭并行计算
stopImplicitCluster()

# 输出最终的ColocData结果
colocdata <- colocdata_list

# 输出 colocdata
write.table(colocdata, file = "coloc/result/metab_eqtlgen.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

