library(data.table)
rm(list=ls())

loci=fread("loci.txt")
loci=loci[,c(17:19)]
loci=unique(loci)
head(loci)

library(dplyr)
library(ggplot2)

# 统计 unique_tags 每个数字出现的次数
unique_tags_count <- loci %>%
  group_by(unique_tags) %>%
  summarise(count = n())

# 手动设置因子的级别顺序
unique_tags_count$unique_tags <- factor(unique_tags_count$unique_tags, 
                                        levels = c(as.character(sort(as.numeric(unique_tags_count$unique_tags[unique_tags_count$unique_tags != ">=10"]))), ">=10"))



# 查看统计结果
print(unique_tags_count)

# 创建堆叠柱状图，设置背景为白色，隐藏灰色框线
p=ggplot(unique_tags_count, aes(x = unique_tags, y = count, fill = "blue2")) +
  geom_bar(stat = "identity") +  # 堆叠柱状图
  geom_text(aes(label = count), vjust = -0.3) +
  labs(#title = "The number of metabolites associated with each locus",
       x = "Number of associated metabolites per locus", y = "Number of loci") +  # 设置标题和 y 轴标签，y 轴标签为空
  theme(axis.text.x = element_text(size = 12, hjust = 0.5),  # 调整 x 轴标签的字体大小
        axis.text.y = element_text(size = 12),  # 调整 y 轴标签的字体大小
        axis.title.x = element_text(size = 14),  # 调整 x 轴标题的字体大小
        axis.title.y = element_text(size = 14),  # 调整 y 轴标题的字体大小
        panel.background = element_rect(fill = "white", color = NA),  # 设置背景为白色，隐藏灰色框线
        panel.grid = element_blank(),  # 隐藏网格线
        axis.line = element_line(color = "black"),  # 设置坐标轴线为黑色
        plot.title = element_text(hjust = 0.5),  # 设置标题居中
        legend.position = "none")
p
ggsave("fig2b.png", p, width = 8, height = 5, dpi = 300, bg = "white")  # 调整图像的宽度和高度


