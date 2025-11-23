rm(list=ls())
library(data.table)
library(CMplot)

metab=fread("metabQTL_mrQTL.txt")
metab=metab[,c(1:3,15)]
colnames(metab)=c("Chromosome","Position","SNP","metab")

pro=fread("pQTL_prQTL.txt")
pro=pro[,c("CHR","BP","SNP","P")]
colnames(pro)=c("Chromosome","Position","SNP","protein")

meth=fread("C:/Users/yp/Desktop/合并写文章/甲基化/meQTL_5e8_unique_for_cmplot.txt")
colnames(meth)[2:3]=c("Chromosome","Position")
meth=as.data.frame(meth)


cmdata1=merge(meth[,c(2:4)],pro[,c(1:2,4)],by=c("Chromosome","Position"),all=T)
cmdata=merge(cmdata1,metab[,c(1:2,4)],by=c("Chromosome","Position"),all=T)

cmdata=cmdata[,c(3,1,2,6,5,4)]

pr=fread("prqtl.txt")
mr=fread("mrqtl.txt")

SNPs <- list(
  mr$SNP,
  pr$SNP,
  NA
)


pdf("cmplot_output.pdf", width = 12, height = 12)  # 设置宽度和高度，依据需要调整


CMplot(
  cmdata,
  plot.type = "c",
  cex = 0.6,  # 调整点的大小，0.3表示比默认值小
  r = 0.4,
  outward = FALSE,
  cir.chr.h = 1.3,
  chr.den.col = "black",
  file.output = FALSE,
  verbose = T,
  highlight=SNPs,            # 高亮的 SNP 名称（list 对应多个 trait）
  highlight.type="p",
  highlight.pch = 18                   # 设置为菱形（pch = 5）
)


dev.off()
  











