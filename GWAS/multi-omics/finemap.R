
#install.packages("susieR")

rm(list = ls())
library(susieR)
library(data.table)
######LHPP#########################################
rm(list=ls())
highlight_snp="rs11245086"


make_pip_df_with_cs <- function(susie_fit, snp_vec, random_matrix) {
  
  # --- 基础pip_df ---
  pip_df <- data.frame(
    SNP = snp_vec,
    PIP = susie_fit$pip,
    index = seq_along(snp_vec)
  )
  
  # --- 如果没有CS（保险处理） ---
  if (is.null(susie_fit$sets$cs)) {
    return(pip_df)
  }
  
  # --- 提取CS信息 ---
  cs_summary <- summary(susie_fit)$cs
  cs_list <- susie_fit$sets$cs
  coverage <- susie_fit$sets$coverage
  
  # --- 初始化列 ---
  pip_df$cs_id <- NA
  pip_df$cs_log10bf <- NA
  pip_df$cs_avg_r2 <- NA
  pip_df$cs_min_r2 <- NA
  pip_df$cs_coverage <- NA
  
  # --- 填充 ---
  for (i in seq_along(cs_list)) {
    
    snp_idx <- cs_list[[i]]
    
    pip_df$cs_id[snp_idx] <- cs_summary$cs[i]
    pip_df$cs_log10bf[snp_idx] <- cs_summary$cs_log10bf[i]
    pip_df$cs_avg_r2[snp_idx] <- cs_summary$cs_avg_r2[i]
    pip_df$cs_min_r2[snp_idx] <- cs_summary$cs_min_r2[i]
    pip_df$cs_coverage[snp_idx] <- coverage[i]
  }
  
  # --- 是否在CS ---
  pip_df$in_cs <- !is.na(pip_df$cs_id)
  
  #r2
  pip_df=merge(pip_df,random_matrix,by="SNP")
  
  return(pip_df)
}

#####读入原始数据
aa=fread("localplot/cg_lhpp.txt")
bb=fread("localplot/LHPP.txt")
cc=fread("localplot/Q9H008_rs11245086.txt")

#####计算z
aa$z=aa$BETA/aa$SE
bb$z=bb$BETA/bb$SE
cc$z=cc$BETA/cc$SE

#########################

aa=aa[,c("SNP_Chr","SNP_Bp","SNP","z","N_sample")]
colnames(aa)=c("CHR","BP","SNP","z","N")
bb=bb[,c("CHR","BP","SNP","z","N")]
cc=cc[,c("CHR","BP","SNP","z","N")]


#####读入ld r2
random_matrix=fread("localplot/lhppmatrix.txt")
random_matrix=as.data.frame(random_matrix)
random_matrix$SNP=colnames(random_matrix)
random_matrix=random_matrix[,c(highlight_snp,"SNP")]
colnames(random_matrix)[1]="r2"



######################for-cpg####

aar=merge(aa,random_matrix,by="SNP")

ld_matrix=fread("localplot/lhppmatrix.txt")
ld_matrix <- as.data.frame(ld_matrix)
rownames(ld_matrix)=colnames(ld_matrix)

ld_matrix2 <- ld_matrix[aar$SNP, aar$SNP]

all(aar$SNP == rownames(ld_matrix2))
ld_matrix2=as.matrix(ld_matrix2)

####susier
nnn=max(aar$N)

fitted_rss3 <- susie_rss(aar$z, ld_matrix2, n=nnn, L=10, refine=T, max_iter = 1000)
susie_plot(fitted_rss3, y="PIP")
summary(fitted_rss3)$cs

pip_df <- make_pip_df_with_cs(
  susie_fit = fitted_rss3,
  snp_vec = aar$SNP,
  random_matrix=random_matrix
)
pip_df$trait="LHPP_cg10700560"

######################for-gene####

ld_matrix=fread("eur_1kg_ld/lhppmatrix.txt")
ld_matrix <- as.data.frame(ld_matrix)
rownames(ld_matrix)=colnames(ld_matrix)
ld_matrix1=ld_matrix
ld_matrix1$SNP=colnames(ld_matrix1)
ld_matrix1=ld_matrix1[,c(highlight_snp,"SNP")]


bbr=merge(bb,ld_matrix1,by="SNP")

ld_matrix2 <- ld_matrix[bbr$SNP, bbr$SNP]

all(bbr$SNP == rownames(ld_matrix2))
ld_matrix2=as.matrix(ld_matrix2)

####susier
nnn=max(bbr$N)

fitted_rss3 <- susie_rss(bbr$z, ld_matrix2, n=nnn, L=10, refine=T, max_iter = 1000)
susie_plot(fitted_rss3, y="PIP")
summary(fitted_rss3)$cs

pip_df2 <- make_pip_df_with_cs(
  susie_fit = fitted_rss3,
  snp_vec = bbr$SNP,
  random_matrix=random_matrix
)
pip_df2$trait="LHPP_mRNA"



######################for-protein####

ccr=merge(cc,random_matrix,by="SNP")

ld_matrix=fread("localplot/lhppmatrix.txt")
ld_matrix <- as.data.frame(ld_matrix)
rownames(ld_matrix)=colnames(ld_matrix)

ld_matrix2 <- ld_matrix[ccr$SNP, ccr$SNP]

all(ccr$SNP == rownames(ld_matrix2))
ld_matrix2=as.matrix(ld_matrix2)

####susier
nnn=max(ccr$N)

fitted_rss3 <- susie_rss(ccr$z, ld_matrix2, n=nnn, L=10, refine=T, max_iter = 1000)
susie_plot(fitted_rss3, y="PIP")
summary(fitted_rss3)$cs

pip_df3 <- make_pip_df_with_cs(
  susie_fit = fitted_rss3,
  snp_vec = ccr$SNP,
  random_matrix=random_matrix
)
pip_df3$trait="LHPP_protein"



all_pip=rbind(pip_df,pip_df2,pip_df3)

all_pip=all_pip[!is.na(all_pip$cs_id),]

fwrite(all_pip,"finemap/result_lhpp.csv")









