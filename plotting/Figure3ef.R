
#############################

rm(list=ls())
library(data.table)
library(ggplot2)
maf=fread("metab_eur_maf.csv")

library(dplyr)



maf <- na.omit(maf)


maf=unique(maf)
sub_maf=maf[,c("A1_FREQ","eur_freq")]
colnames(sub_maf)=c("freq","eumaf")

sub_maf <- sub_maf %>%
  mutate(eumaf = ifelse(eumaf == 0, 0.00001, eumaf))

sub_maf <- sub_maf %>%
  mutate(more = ifelse(10*eumaf < freq, 1, 0))




# 创建 ggplot 图
p=ggplot(data = sub_maf, aes(x = eumaf, y = freq)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  
  geom_point(aes(color = as.factor(more))) +
  scale_color_manual(values = c("0" = "gray", "1" = "blue")) +
  
  scale_x_continuous(trans = "log10", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5),
                     labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), "0.1", "0.5"),
                     limits = c(1e-5, 0.5)) +
  
  scale_y_continuous(trans = "log10", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5),
                     labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), "0.1", "0.5"),
                     limits = c(1e-2, 0.5)) +
  
  labs(x = "MAF in EUR", y = "MAF in CAS")+
  theme_bw()+ 
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),  
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")  

p
ggsave("CASMAF vs EURAMF.png", p, width = 8, height = 5, dpi = 300, bg = "white")  # 调整图像的宽度和高度

#############################

rm(list=ls())
library(data.table)
library(ggplot2)
maf=fread("metab_eur_maf.csv")

library(dplyr)



maf <- na.omit(maf)


maf=unique(maf)
sub_maf=maf[,c("A1_FREQ","eur_freq","Reported/Novel")]
colnames(sub_maf)=c("freq","eumaf","more")

sub_maf <- sub_maf %>%
  mutate(eumaf = ifelse(eumaf == 0, 0.00001, eumaf))





# 创建 ggplot 图
p=ggplot(data = sub_maf, aes(x = eumaf, y = freq)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +

  geom_point(aes(color = as.factor(more))) +
  scale_color_manual(values = c("Reported" = "gray", "Novel" = "#FCA311")) +

  scale_x_continuous(trans = "log10", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5),
                     labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), "0.1", "0.5"),
                     limits = c(1e-5, 0.5)) +
 
  scale_y_continuous(trans = "log10", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5),
                     labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), "0.1", "0.5"),
                     limits = c(1e-2, 0.5)) +
 
  labs(x = "MAF in EUR", y = "MAF in CAS")+
  theme_bw()+  
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),  
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none")   

p
ggsave("CASMAF vs EURAMF_novel.png", p, width = 8, height = 5, dpi = 300, bg = "white")  # 调整图像的宽度和高度

