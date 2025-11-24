

#########提琴图###########
rm(list=ls())

hsnp=read.table("hsnp.txt",header = T,fill = T)
info=read.table("代谢物信息.txt", header=T,
              sep = "\t",quote = "",fill = T)
hsnp$tag <- sapply(strsplit(hsnp$file, "\\."), function(x) x[1])

aa=merge(hsnp,info,by.x="tag",by.y="m_id")
aamid=median(aa$hsnp)

# 
median_values <- aggregate(hsnp ~ SuperClass, data = aa, FUN = median)


p <- ggplot(aa, aes(x = SuperClass, y = hsnp, fill = SuperClass)) +
  geom_violin(width = 1) +  
  theme_minimal() +
  scale_fill_manual(values = c("Alkaloids"="#a6d854","Benzenoids"="#e78ac3",
                               "Lipids and lipid-like molecules"="#66c2a5","Nucleosides"="#ffd92f",
                               "Organic acids and derivatives"="#fc8d62","Organic nitrogen compounds"="#8da0cb",
                               "Organic oxygen compounds"="chartreuse4","Organoheterocyclic compounds"="#e5c494",
                               "Phenylpropanoids and polyketides"="cadetblue4","Unknown"="#1f78b4")) +
  labs(#title = "Metabolite heritability explained by assayed genotypes", 
       x = "SuperClass", y = "Heritability explained by genetic variants") +
  guides(fill = "none") +  
  geom_crossbar(data = median_values, aes(ymin = hsnp, ymax = hsnp), 
                position = position_dodge(width = 0.75), width = 0.5)+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        panel.background = element_rect(fill = "white", color = NA), 
        panel.grid = element_blank(),  
        axis.line = element_line(color = "black"), 
        plot.title = element_text(hjust = 0.5),  
        legend.position = "none") + 
  geom_hline(yintercept = median(aa$hsnp), linetype = "dashed", color = "blue")  
p
ggsave("fig3_violin.png", p, width = 8, height = 5, dpi = 300, bg = "white")  

