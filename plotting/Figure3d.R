install.packages('RIdeogram')
require(RIdeogram)
library(data.table)
library(dplyr)
rm(list=ls())

human_karyotype1=fread("每条染色体的长度CHR.txt")
data(human_karyotype, package="RIdeogram")
human_karyotype=human_karyotype[1:22,]
human_karyotype[,1:3]=human_karyotype1[,1:3]


label=fread("metab的novel的结果用于画染色体图.txt")
label$Chr=as.character(label$Chr)
#label=label[1:10,]
# 假设 label 是你的数据框
label_dedup <- label %>%
  group_by(Chr, End) %>%
  arrange(Type) %>%  # 将 Type 按字母顺序排列，"novel" 会在 "known" 之前
  filter(row_number() == 1) %>%  # 只保留每组的第一行
  ungroup() %>%  # 解除分组
  distinct()  # 去重
#####
label_dedup$color[label_dedup$color=="33a02c"]<-"2A9D8F"
label_dedup$color[label_dedup$color=="6a3d9a"]<-"FCA311"
label_dedup$Type[label_dedup$Type == "known"] <- "Known"
label_dedup$Type[label_dedup$Type == "novel"] <- "Novel"
label_dedup$Shape="circle"
ideogram(karyotype = human_karyotype, label = label_dedup, label_type = "marker")

convertSVG("chromosome.svg",file = "a1" ,device = "pdf",width=8,height=5, dpi = 600)


