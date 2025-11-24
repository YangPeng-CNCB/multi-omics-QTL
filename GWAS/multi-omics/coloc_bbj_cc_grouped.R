library(data.table)


library(foreach)
library(doParallel)
library("coloc")
library(dplyr)

num_cores <- 30
registerDoParallel(cores = num_cores)


path <- fread("mr_bbj220/coloc/1coloc_path.txt")
path <- filter(path, type == "Binary")


eqtl_input <- path$path
metab_input <- path$file_path


######################

results <- foreach(i = 1:nrow(path), .combine = 'rbind', .packages = c("data.table", "coloc", "dplyr")) %dopar% {
  tryCatch({
    #
    eqtl_file <- eqtl_input[i]
    metab_file <- metab_input[i]
    
    #
    eqtl <- fread(file = eqtl_file)
    metab <- fread(file = metab_file)
    
    #
    ss <- as.numeric(path$s[i])
    bbj_name <- path$outcome[i] 
    
    #
    colnames(eqtl)[3:8] <- c("BP", "SNP", "A2", "A1", "A1C", "FRQ")
    colnames(eqtl)[14] <- "P"
    eqtl$CHR <- as.numeric(eqtl$CHR)
    
    #
    input <- merge(eqtl, metab, by = c("CHR", "BP"), all = FALSE, suffixes = c("_eqtl", "_metab"))
    
    if (nrow(input) > 0) {
      #
      input <- input[!duplicated(input$SNP_eqtl), ]
      input <- input[!duplicated(input$SNP_metab), ]
      input$SNP_eqtl <- paste(input$SNP_eqtl, bbj_name, input$CHR, input$BP, "P_bbj", input$P_eqtl, "P_metab", input$P_metab, sep = ";")
      input$SNP_metab <- input$SNP_eqtl
      
      #
      result <- coloc.abf(
        dataset1 = list(beta = input$BETA_eqtl, varbeta = (input$SE_eqtl)^2, pvalues = input$P_eqtl, N = input$N_eqtl, 
                        MAF = input$FRQ, type = "cc", s = ss, snp = input$SNP_eqtl, position = input$BP),
        dataset2 = list(beta = input$BETA_metab, varbeta = (input$SE_metab)^2, pvalues = input$P_metab, N = input$N_metab, 
                        MAF = input$MAF, type = "quant", snp = input$SNP_metab, position = input$BP),
        p12 = 1e-5
      )  
      
      #
      pp4_value <- max(result$results$SNP.PP.H4)
      
      #
      colocdata <- result$results[which.max(result$results$SNP.PP.H4), ]
      colocdata$EqtlFile <- basename(eqtl_file)
      colocdata$MetabFile <- basename(metab_file)
      colocdata$noverlap <- nrow(result$results)
      sum1=result$summary
      #
      sum1 <- as.data.frame(t(sum1))
      
      #
      return(cbind(
        sum1,
        data.frame(EqtlFile = basename(eqtl_file), MetabFile = basename(metab_file), PP4 = pp4_value),
        colocdata
      ))
      
    } else {
      #
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
        SNP = character(0),  
        position = integer(0),  
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
    
    cat("Error in processing files:", eqtl_file, ";", metab_file, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    NULL
  })
}

stopImplicitCluster()

write.table(results, file = "mr_bbj220/coloc/result/summary_bbj_cc_group.txt",
            append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

