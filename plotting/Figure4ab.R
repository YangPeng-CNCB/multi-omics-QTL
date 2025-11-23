

#################画图c
rm(list=ls())
library(data.table)

aa=fread("coloc/用于画图的数据/画图a的数据.txt")

aa$disconsistence=1-aa$consistence

library(tidyr)

# 转换数据为长格式
data_long <- aa %>%
  pivot_longer(cols = c(consistence, disconsistence), 
               names_to = "type", 
               values_to = "value")

data_long$value=round(data_long$value*100, 1)


library(ggplot2)


p=ggplot(data_long, aes(x = Dataset, y = value, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(x = "", y = "Percentage (%)", fill = "Type") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(),      # 去除主网格线
        panel.grid.minor = element_blank(),      # 去除次网格线
        axis.line = element_line(color = "black")) +  # 添加x轴和y轴
  scale_x_discrete(labels = function(x) gsub("\\(", "\n(", x),
                   limits = unique(data_long$Dataset)) +
  scale_fill_manual(values = c("consistence" = "#fc8d62", "disconsistence" = "#a6d854")) +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5, size = 5)
p
ggsave("coloc/用于画图的数据/fig4a.png", p,
       width = 8, height = 5.5, dpi = 300, bg = "white")  # 调整图像的宽度和高度


#################画图d
rm(list=ls())
library(data.table)

aa=fread("coloc/用于画图的数据/画图b的数据.txt")
aa=aa[1:3,]
aa$ncorr=1-aa$corr

library(tidyr)

# 转换数据为长格式
data_long <- aa %>%
  pivot_longer(cols = c(corr, ncorr), 
               names_to = "type", 
               values_to = "value")

data_long$value=round(data_long$value*100, 1)


library(ggplot2)


p=ggplot(data_long, aes(x = Dataset, y = value, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(x = "", y = "Percentage (%)", fill = "Type") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(),      # 去除主网格线
        panel.grid.minor = element_blank(),      # 去除次网格线
        axis.line = element_line(color = "black")) +  # 添加x轴和y轴
  scale_x_discrete(labels = function(x) gsub("\\(", "\n(", x),
                   limits = unique(data_long$Dataset)) +
  scale_fill_manual(values = c("Positive" = "#ffd92f", "Negative" = "#1f78b4")) +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.5, size = 5)
p
ggsave("coloc/用于画图的数据/fig4b.png", p,
       width = 8, height = 5.5, dpi = 300, bg = "white")  # 调整图像的宽度和高度



