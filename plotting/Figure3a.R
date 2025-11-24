rm(list=ls())

#######

class =read.table("代谢物信息.txt", header=T,
                     sep = "\t",quote = "",fill = T)
superclass=class[c("m_id","SuperClass")]
colnames(superclass)=c("tag","superclass")

AB<-read.table(file = "m841_na_1e-5.txt",header=T)

AB$tag <- sub("^[^.]+\\.(.*?)\\..*$", "\\1", AB$tag)

colnames(AB)[1] <- "CHR"
AB$AB[AB$CHR == 1] <- AB$POS[AB$CHR == 1] 
AB$AB[AB$CHR == 2] <- AB$POS[AB$CHR == 2] +248956422
AB$AB[AB$CHR == 3] <- AB$POS[AB$CHR == 3] +491149951
AB$AB[AB$CHR == 4] <- AB$POS[AB$CHR == 4] +689445510
AB$AB[AB$CHR == 5] <- AB$POS[AB$CHR == 5] +879660065
AB$AB[AB$CHR == 6] <- AB$POS[AB$CHR == 6] +1061198324
AB$AB[AB$CHR == 7] <- AB$POS[AB$CHR == 7] +1232004303
AB$AB[AB$CHR == 8] <- AB$POS[AB$CHR == 8] +1391350276
AB$AB[AB$CHR == 9] <- AB$POS[AB$CHR == 9] +1536488912
AB$AB[AB$CHR == 10] <- AB$POS[AB$CHR == 10] +1674883629
AB$AB[AB$CHR == 11] <- AB$POS[AB$CHR == 11] +1808681051
AB$AB[AB$CHR == 12] <- AB$POS[AB$CHR == 12] +1943767673
AB$AB[AB$CHR == 13] <- AB$POS[AB$CHR == 13] +2077042982
AB$AB[AB$CHR == 14] <- AB$POS[AB$CHR == 14] +2191407310
AB$AB[AB$CHR == 15] <- AB$POS[AB$CHR == 15] +2298451028
AB$AB[AB$CHR == 16] <- AB$POS[AB$CHR == 16] +2400442217
AB$AB[AB$CHR == 17] <- AB$POS[AB$CHR == 17] +2490780562
AB$AB[AB$CHR == 18] <- AB$POS[AB$CHR == 18] +2574038003
AB$AB[AB$CHR == 19] <- AB$POS[AB$CHR == 19] +2654411288
AB$AB[AB$CHR == 20] <- AB$POS[AB$CHR == 20] +2713028904
AB$AB[AB$CHR == 21] <- AB$POS[AB$CHR == 21] +2777473071
AB$AB[AB$CHR == 22] <- AB$POS[AB$CHR == 22] +2824183054



AB=merge(AB,superclass,by="tag")

selected_columns <- AB[c("AB", "P", "superclass")]


# 创建自定义横轴刻度
custom_ticks <- c(0, 248956422, 491149951,689445510,879660065,1061198324,1232004303,1391350276,
                  1536488912,1674883629,1808681051,1943767673,2077042982,2191407310,2298451028,
                  2400442217,2490780562,2574038003,2654411288,2713028904,2777473071,2824183054,
                  2875001522)

# 创建自定义刻度标签
custom_labels <- c(1:23)

# 额外添加的刻度位置和标签
extra_ticks <- c(124478211, 370053187, 590297731, 784552788, 970429195, 1146601314, 1311677290, 1463919594, 1605686271,
                 1741782340, 1876224362, 2010405328, 2134225146, 2244929169, 2349446623, 2445611390, 2532409283, 2614224646,
                 2683720096, 2745250988, 2800828063, 2849592288)

extra_labels <- c(1:22)

# 定义偏移量
offset <- 10000000


library(ggplot2)



p <- ggplot(AB, aes(x = AB, y = -log10(P), color = superclass)) +
  geom_point(size = 1, shape = 19) +
  scale_color_manual(values = c("Alkaloids"="chartreuse4","Benzenoids"="#e78ac3",
                                "Lipids and lipid-like molecules"="#66c2a5","Nucleosides"="#ffd92f",
                                "Organic acids and derivatives"="#fc8d62","Organic nitrogen compounds"="#8da0cb",
                                "Organic oxygen compounds"="#a6d854","Organoheterocyclic compounds"="#e5c494",
                                "Phenylpropanoids and polyketides"="cadetblue4","Unknown"="#1f78b4",
                                "z1" = "lightgrey", "z2" = "darkgray")) +
  labs(x = "Chromosome", y = expression(-phantom() * log[10](P))) + 
  scale_x_continuous(breaks = custom_ticks, labels = custom_labels) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black", size = 1),  
    axis.text.y = element_text(margin = margin(r = 10), size = 20),  
    axis.ticks.y = element_line(size = 1),  
    axis.title = element_text(size = 22),  
    plot.title = element_text(hjust = 0.5, size = 16),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14),  
    legend.key.size = unit(1.5, "lines"),  
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank()   
    )+
  guides(color = guide_legend(title = "SuperClass")) +  
  geom_hline(yintercept = -log10(5e-8), linetype = "solid", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1)  +
  annotate("text", x = extra_ticks + offset, y = -0.5, label = extra_labels, size = 6, hjust = 0.5, vjust = 1.5)+  
  annotate("text", x = 1437500761, y = -6, label = "Chromosome", size = 9, hjust = 0.5, vjust = 1.5)  

#p
ggsave("曼哈顿图.png", p, width = 20, height = 9, limitsize = FALSE, bg = "white")







