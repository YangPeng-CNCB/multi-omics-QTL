
#####boxplot###########
rm(list=ls())


geno=fread("readable_geno.raw")
colnames(geno) <- gsub("_.*", "", colnames(geno))
geno=geno[,-3:-6]

pro=fread("pro_input.txt")

map=fread("readable_geno.map")
map=map[,c(1,2,4)]
colnames(map)=c("CHR","snp","BP")


geno=geno[,c("FID","IID","rs2210912")]
pro=pro[,c("FID","IID","Q96P31")]

aa=merge(geno,pro,by=c("FID","IID"))

library(ggplot2)

# 
aa$rs2210912 <- factor(aa$rs2210912, 
                       levels = c(0, 1, 2), 
                       labels = c("AA (372)", "CA (519)", "CC (163)"))


p=ggplot(aa, aes(x = rs2210912, y = Q96P31)) +
  geom_boxplot(fill = "white", color = "black", outlier.colour = "#2A9D8F") +  
  theme_classic() +  
  labs(#title = "", 
    x = "Genotype at rs2210912", 
    y = "Protein FCRL3 level") +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 14, color = "black"),  
    plot.title = element_text(size = 16, hjust = 0.5, color = "black")  
  )
p
ggsave("boxplot.png", p,
       width = 6, height = 6, dpi = 300, bg = "white")  





#######################


rm(list=ls())

library(dplyr)
library(geni.plots)

##############fcrl3_gd#######################

rm(list=ls())
library(data.table)
aa=fread("fcrl3_gd.txt")
bb=fread("Q96P31_rs2210912.txt")
bb <- bb %>% distinct(ID, .keep_all = TRUE)
aa <- aa %>% distinct(SNPID, .keep_all = TRUE)

aa=aa[,c("CHR","POS","p.value")]
bb=bb[,c("CHR","BP","P","SNP")]

d1=merge(aa,bb,by.x=c("CHR","POS"),by.y=c("CHR","BP"))
d1 <- d1 %>% distinct(SNP, .keep_all = TRUE)
merged_data=d1[,c(5,1,2,4,3)]
colnames(merged_data)=c("marker","chr","pos" ,"pvalue_1","pvalue_2")

random_matrix=fread("fcrl3_2matrix.txt")
random_matrix=as.data.frame(random_matrix)
rownames(random_matrix)=colnames(random_matrix)

# 
selected_markers <- merged_data$marker

# 
selected_markers <- selected_markers[selected_markers %in% rownames(random_matrix) & 
                                       selected_markers %in% colnames(random_matrix)]


filtered_matrix <- random_matrix[selected_markers, selected_markers, drop = FALSE]

merged_data=merged_data[merged_data$marker %in% colnames(filtered_matrix),]


library(ggplot2)
p <- fig_region_stack(
  data = merged_data,
  traits = c("FCRL3","Graves’ disease"),
  corr = filtered_matrix,
  r2 = TRUE,
  point_size = 2, 
  build = 37,
  alpha=0.8,  
  top_marker = "rs2210912",
  title_center = T,
  thresh=c(5e-8),
  thresh_colour = "#e59027"
)

ggsave("fcrl3_gd.png", 
       plot = p, width = 6.2, height = 6.2, units = "in", dpi = 300)


#####cc_ugt2b17_dca###########
rm(list=ls())
library(data.table)
cc=fread("cc_4.txt")
gene=fread("UGT2B17.txt")
metab=fread("m317.txt")

cc=cc[,c("CHR","POS","p.value")]
gene=gene[,c("CHR","BP","P")]
metab=metab[,c("SNP","CHR","BP","P")]

d1=merge(metab,cc,by.x=c("CHR","BP"),by.y=c("CHR","POS"))

d2=merge(d1,gene,by=c("CHR","BP"), all.x = TRUE)

merged_data=d2[,c(3,1,2,6,4,5)]
colnames(merged_data)=c("marker","chr","pos" ,"pvalue_1","pvalue_2","pvalue_3")


random_matrix=fread("ugtmatrix.txt")
random_matrix=as.data.frame(random_matrix)
rownames(random_matrix)=colnames(random_matrix)

# 
selected_markers <- merged_data$marker

# 
selected_markers <- selected_markers[selected_markers %in% rownames(random_matrix) & 
                                       selected_markers %in% colnames(random_matrix)]


filtered_matrix <- random_matrix[selected_markers, selected_markers, drop = FALSE]


p <- fig_region_stack(
  data = merged_data,
  traits = c("UGT2B17", "Deoxycholic acid 3-glucuronide","Colorectal cancer"),
  corr = filtered_matrix,
  r2 = TRUE,
  point_size = 2, 
  build = 37,
  alpha=0.8,  
  top_marker = "rs2603153",
  title_center = T,
  thresh=c(5e-8),
  thresh_colour = "#e59027"
)
p
library(ggplot2)
ggsave("crc-ugt2b17-dca.png", 
       plot = p, width = 6.2, height = 8.2, units = "in", dpi = 300)



