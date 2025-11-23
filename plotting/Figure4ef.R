
rm(list=ls())

library(dplyr)
library(geni.plots)


########cpg-protein#################

rm(list=ls())
library(data.table)
aa=fread("cg_lhpp.txt")
bb=fread("LHPP.txt")
cc=fread("Q9H008_rs11245086.txt")

aa=aa[,c("SNP_Chr","SNP_Bp","P_value")]
bb=bb[,c("CHR","BP","P")]
cc=cc[,c("SNP","CHR","BP","P")]

d1=merge(aa,bb,by.x=c("SNP_Chr","SNP_Bp"),by.y=c("CHR","BP"))

d2=merge(d1,cc,by.x=c("SNP_Chr","SNP_Bp"),by.y=c("CHR","BP"))

merged_data=d2[,c(5,1,2,3,4,6)]
colnames(merged_data)=c("marker","chr","pos" ,"pvalue_1","pvalue_2","pvalue_3")

random_matrix=fread("lhppmatrix.txt")
random_matrix=as.data.frame(random_matrix)
rownames(random_matrix)=colnames(random_matrix)

# 提取 merged_data 中的 marker 名称
selected_markers <- merged_data$marker

# 确保 random_matrix 的行名和列名都包含这些 marker
selected_markers <- selected_markers[selected_markers %in% rownames(random_matrix) & 
                                       selected_markers %in% colnames(random_matrix)]

# 提取对应的行和列
filtered_matrix <- random_matrix[selected_markers, selected_markers, drop = FALSE]

library(ggplot2)
p <- fig_region_stack(
  data = merged_data,
  traits = c("cg10700560", "LHPP (mRNA)","LHPP (protein)"),
  corr = filtered_matrix,
  r2 = TRUE,
  point_size = 2, 
  build = 37,
  alpha=0.8,  
  top_marker = "rs11245086",
  title_center = T,
  thresh=c(5e-8),
  thresh_colour = "#e59027"
)

ggsave("fig4e.png", 
       plot = p, width = 6.2, height = 8.2, units = "in", dpi = 300)


##############cpg-metab#######################

rm(list=ls())
library(data.table)
aa=fread("cg_fads2.txt")
bb=fread("FADS2.txt")
cc=fread("m209_rs174559.txt")

aa=aa[,c("SNP_Chr","SNP_Bp","P_value")]
bb=bb[,c("CHR","BP","P")]
cc=cc[,c("SNP","CHR","BP","P")]

d1=merge(aa,bb,by.x=c("SNP_Chr","SNP_Bp"),by.y=c("CHR","BP"))

d2=merge(d1,cc,by.x=c("SNP_Chr","SNP_Bp"),by.y=c("CHR","BP"))

merged_data=d2[,c(5,1,2,3,4,6)]
colnames(merged_data)=c("marker","chr","pos" ,"pvalue_1","pvalue_2","pvalue_3")

random_matrix=fread("fads2matrix.txt")
random_matrix=as.data.frame(random_matrix)
rownames(random_matrix)=colnames(random_matrix)

# 提取 merged_data 中的 marker 名称
selected_markers <- merged_data$marker

# 确保 random_matrix 的行名和列名都包含这些 marker
selected_markers <- selected_markers[selected_markers %in% rownames(random_matrix) & 
                                       selected_markers %in% colnames(random_matrix)]

# 提取对应的行和列
filtered_matrix <- random_matrix[selected_markers, selected_markers, drop = FALSE]

library(ggplot2)
p <- fig_region_stack(
  data = merged_data,
  traits = c("cg21029357", "FADS2 (mRNA)","Docosapentaenoic acid"),
  corr = filtered_matrix,
  r2 = TRUE,
  point_size = 2, 
  build = 37,
  alpha=0.8,  
  top_marker = "rs174559",
  title_center = T,
  thresh=c(5e-8),
  thresh_colour = "#e59027"
)

ggsave("fig4f.png", 
       plot = p, width = 6.2, height = 8.2, units = "in", dpi = 300)

