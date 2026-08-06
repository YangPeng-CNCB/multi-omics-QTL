
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


rm(list = ls())
{
  
  library(data.table)
  library(locuszoomr)
  library(EnsDb.Hsapiens.v75)
  library(ggplot2)
  library(rtracklayer)
  library(patchwork)
  
  
  ################设置主体
  
  
  
  border_theme <- theme_bw() +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.6),
      panel.grid = element_blank(),
      plot.title = element_text(hjust=0.5, size=12, face="bold")
    )
  
  no_x <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
  
  sig_line <- geom_hline(
    yintercept = -log10(5e-8),
    colour = "#e59027",
    linewidth = 0.6,
    linetype = "dashed"
  )
  
  LD_colors <- c(
    "grey85",     # no LD
    "lightblue2",    # 0-0.2
    "lightgreen",    # 0.2-0.4
    "yellow2",    # 0.4-0.6
    "orange2",    # 0.6-0.8
    "red3",    # 0.8-1
    "purple"      # index SNP
  )
  
  
  recomb.hg19 <- import.bw("hapMapRelease24CombinedRecombMap.bw")
  
  
}

########cpg-protein#################


highlight_snp="rs2210912"
#####读入原始数据
aa=fread("FCRL3.txt")
bb=fread("Q96P31_rs2210912.txt")
cc=fread("fcrl3_gd.txt")


#####读入ld r2
random_matrix=fread("fcrl3_2matrix.txt")
random_matrix=as.data.frame(random_matrix)
random_matrix$SNP=colnames(random_matrix)
random_matrix=random_matrix[,c(highlight_snp,"SNP")]
colnames(random_matrix)[1]="r2"


#########################

aa=aa[,c("CHR","BP","SNP","P")]
bb=bb[,c("CHR","BP","SNP","P")]
cc=cc[,c("CHR","POS","p.value")]
cc1=merge(bb[,1:3],cc,by.x=c("CHR","BP"),by.y=c("CHR","POS"))

##########


aar=merge(aa,random_matrix,by="SNP")
bbr=merge(bb,random_matrix,by="SNP")
ccr=merge(cc1,random_matrix,by="SNP")


##########
loc1 <- locus(aar, gene="FCRL3", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc2 <- locus(bbr, gene="FCRL3", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc3 <- locus(ccr, gene="FCRL3", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "p.value",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc1 <- link_recomb(loc1, recomb = recomb.hg19)
loc2 <- link_recomb(loc2, recomb = recomb.hg19)
loc3 <- link_recomb(loc3, recomb = recomb.hg19)



p1 <- gg_scatter(loc1,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 40) +
  sig_line +
  labs(title="FCRL3(mRNA)") +
  border_theme +
  no_x

p2 <- gg_scatter(loc2,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 20) +
  sig_line +
  labs(title="FCRL3(protein)") +
  border_theme +
  no_x

p3 <- gg_scatter(loc3,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 2) +
  sig_line +
  labs(title="Graves’ disease") +
  border_theme +
  no_x

g <- gg_genetracks(loc3, highlight = "FCRL3", maxrows = 3,filter_gene_biotype = 'protein_coding',
                   cex.text=0.8) +
  border_theme


combined_plot <- (p1 / p2 / p3 / g) +
  plot_layout(heights = c(3,3,3,2))


ggsave("fig6a_fcrl3.png", 
       plot = combined_plot, width = 7, height = 7.5, units = "in", dpi = 300)


#####cc_ugt2b17_dca###########
########cpg-protein#################


highlight_snp="rs2603153"
#####读入原始数据
aa=fread("UGT2B17.txt")
bb=fread("m317.txt")
cc=fread("cc_4.txt")


#####读入ld r2
random_matrix=fread("ugtmatrix.txt")
random_matrix=as.data.frame(random_matrix)
random_matrix$SNP=colnames(random_matrix)
random_matrix=random_matrix[,c(highlight_snp,"SNP")]
colnames(random_matrix)[1]="r2"


#########################

aa=aa[,c("CHR","BP","SNP","P")]
bb=bb[,c("CHR","BP","SNP","P")]
cc=cc[,c("CHR","POS","p.value")]
cc1=merge(bb[,1:3],cc,by.x=c("CHR","BP"),by.y=c("CHR","POS"))

##########


aar=merge(aa,random_matrix,by="SNP")
bbr=merge(bb,random_matrix,by="SNP")
ccr=merge(cc1,random_matrix,by="SNP")


##########
loc1 <- locus(aar, gene="UGT2B17", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc2 <- locus(bbr, gene="UGT2B17", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc3 <- locus(ccr, gene="UGT2B17", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "p.value",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc1 <- link_recomb(loc1, recomb = recomb.hg19)
loc2 <- link_recomb(loc2, recomb = recomb.hg19)
loc3 <- link_recomb(loc3, recomb = recomb.hg19)



p1 <- gg_scatter(loc1, index_snp = "rs145450963", labels="rs145450963",  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 20) +
  sig_line +
  labs(title="UGT2B17") +
  border_theme +
  no_x

p2 <- gg_scatter(loc2,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 10) +
  sig_line +
  labs(title="Deoxycholic acid 3-glucuronide") +
  border_theme +
  no_x

p3 <- gg_scatter(loc3,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 1) +
  sig_line +
  labs(title="Colorectal cancer") +
  border_theme +
  no_x

g <- gg_genetracks(loc3, highlight = "UGT2B17", maxrows = 3,filter_gene_biotype = 'protein_coding',
                   cex.text=0.8) +
  border_theme


combined_plot <- (p1 / p2 / p3 / g) +
  plot_layout(heights = c(3,3,3,2))


ggsave("fig6c_ugt2b17.png", 
       plot = combined_plot, width = 7, height = 7.5, units = "in", dpi = 300)



