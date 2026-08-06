
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

#rm(list=ls())
highlight_snp="rs11245086"
#####读入原始数据
aa=fread("localplot/cg_lhpp.txt")
bb=fread("localplot/LHPP.txt")
cc=fread("localplot/Q9H008_rs11245086.txt")


#####读入ld r2
random_matrix=fread("localplot/lhppmatrix.txt")
random_matrix=as.data.frame(random_matrix)
random_matrix$SNP=colnames(random_matrix)
random_matrix=random_matrix[,c(highlight_snp,"SNP")]
colnames(random_matrix)[1]="r2"


#########################

aa=aa[,c("SNP_Chr","SNP_Bp","SNP","P_value")]
bb=bb[,c("CHR","BP","SNP","P")]
cc=cc[,c("CHR","BP","SNP","P")]

##########


aar=merge(aa,random_matrix,by="SNP")
bbr=merge(bb,random_matrix,by="SNP")
ccr=merge(cc,random_matrix,by="SNP")


##########
loc1 <- locus(aar, gene="LHPP", flank=5e5,  chrom = "SNP_Chr",pos = "SNP_Bp", labs = "SNP", p = "P_value",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc2 <- locus(bbr, gene="LHPP", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc3 <- locus(ccr, gene="LHPP", flank=5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc1 <- link_recomb(loc1, recomb = recomb.hg19)
loc2 <- link_recomb(loc2, recomb = recomb.hg19)
loc3 <- link_recomb(loc3, recomb = recomb.hg19)



p1 <- gg_scatter(loc1,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 1) +
  sig_line +
  labs(title="cg10700560") +
  border_theme +
  no_x

p2 <- gg_scatter(loc2,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 20) +
  sig_line +
  labs(title="LHPP (mRNA)") +
  border_theme +
  no_x

p3 <- gg_scatter(loc3,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 1) +
  sig_line +
  labs(title="LHPP (protein)") +
  border_theme +
  no_x

g <- gg_genetracks(loc3, highlight = "LHPP", maxrows = 3,filter_gene_biotype = 'protein_coding',
                   cex.text=0.8) +
  border_theme

library(patchwork)

combined_plot <- (p1 / p2 / p3 / g) +
  plot_layout(heights = c(3,3,3,2))


ggsave("locizoom/lhpp.png", 
       plot = combined_plot, width = 7, height = 7.5, units = "in", dpi = 300)



##########fads2#####################################

########cpg-protein#################

#rm(list=ls())
highlight_snp="rs174559"
#####读入原始数据
aa=fread("localplot/cg_fads2.txt")
bb=fread("localplot/FADS2.txt")
cc=fread("localplot/m209_rs174559.txt")


#####读入ld r2
random_matrix=fread("localplot/fads2matrix.txt")
random_matrix=as.data.frame(random_matrix)
random_matrix$SNP=colnames(random_matrix)
random_matrix=random_matrix[,c(highlight_snp,"SNP")]
colnames(random_matrix)[1]="r2"


#########################

aa=aa[,c("SNP_Chr","SNP_Bp","SNP","P_value")]
bb=bb[,c("CHR","BP","SNP","P")]
cc=cc[,c("CHR","BP","SNP","P")]

##########


aar=merge(aa,random_matrix,by="SNP")
bbr=merge(bb,random_matrix,by="SNP")
ccr=merge(cc,random_matrix,by="SNP")


##########


loc1 <- locus(aar, gene="FADS2", flank = 5e5,  chrom = "SNP_Chr",pos = "SNP_Bp", labs = "SNP", p = "P_value",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc2 <- locus(bbr, gene="FADS2", flank = 5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc3 <- locus(ccr, gene="FADS2", flank = 5e5,  chrom = "CHR",pos = "BP", labs = "SNP", p = "P",LD = "r2",
              ens_db="EnsDb.Hsapiens.v75")

loc1 <- link_recomb(loc1, recomb = recomb.hg19)
loc2 <- link_recomb(loc2, recomb = recomb.hg19)
loc3 <- link_recomb(loc3, recomb = recomb.hg19)



p1 <- gg_scatter(loc1,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 1 ) +
  sig_line +
  labs(title="cg21029357") +
  border_theme +
  no_x

p2 <- gg_scatter(loc2,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 40) +
  sig_line +
  labs(title="FADS2 (mRNA)") +
  border_theme +
  no_x

p3 <- gg_scatter(loc3,  index_snp = highlight_snp, labels=highlight_snp,  LD_scheme=LD_colors, recomb_col = "#6baed6", nudge_y = 1.5) +
  sig_line +
  labs(title="Docosapentaenoic acid") +
  border_theme +
  no_x

g <- gg_genetracks(loc3, highlight = "FADS2", maxrows = 5,filter_gene_biotype = 'protein_coding',
                   cex.text=0.8) +
  border_theme

library(patchwork)

combined_plot <- (p1 / p2 / p3 / g) +
  plot_layout(heights = c(3,3,3,2))


ggsave("locizoom/fads2.png", 
       plot = combined_plot, width = 7, height = 7.5, units = "in", dpi = 300)
