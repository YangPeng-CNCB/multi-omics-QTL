

rm(list=ls())


library(data.table)
library(ggplot2)

aa=fread("mediation_result.txt")
bb=fread("相关性和偏相关性.txt")

ab=merge(aa,bb,by="tag",all.y = T)

#######画p的分布直方图#######
# 使用 cut() 函数将 sobel_p 列切割成 100 分
aa$sobel_p_cut <- cut(aa$sobel_p, breaks = seq(0, 1, length.out = 101), include.lowest = TRUE)

# 计算每个区间的频数
sobel_p_counts <- table(aa$sobel_p_cut)
sobel_p_data <- data.frame(
  sobel_p_range = seq(0, 1, length.out = 100),
  count = as.numeric(sobel_p_counts)
)



# 绘制柱状分布图
p=ggplot(sobel_p_data, aes(x = sobel_p_range, y = count)) +
  geom_bar(
    stat = "identity", 
    fill = "#1f78b4FF", 
    color = "#1f78b4FF",  # 修改柱子轮廓颜色
    width = 0.01          # 调整柱子宽度
  ) +
  labs(
    x = expression("Sobel " ~ italic(P) ~ "values"),
    y = "Count",
    #title = "Sobel test p-values for Mediation Analysis"  # 添加居中标题
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),          # 去掉背景灰线
    panel.background = element_blank(),    # 去掉背景色
    axis.text = element_text(size = 18),   # 坐标轴字体大小
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),  
    axis.title = element_text(size = 20), # 坐标轴标题样式
    plot.title = element_text(size = 16, hjust = 0.5), # 标题样式
    axis.line = element_line(color = "black"),          # 显示坐标轴线
    axis.ticks = element_line(color = "black")          # 显示坐标轴刻度线
  ) +
  scale_x_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) + # 设置 x 轴刻度
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue", size = 0.8) + # 添加蓝色虚线
  geom_vline(xintercept = 0.05 / 3665, linetype = "dashed", color = "red", size = 0.8) + # 添加红色虚线
  coord_cartesian(ylim = c(0, max(sobel_p_data$count))) # 保证虚线从 y=0 开始
p
ggsave("fig5a.png", p,
       width = 8, height = 5.5, dpi = 300, bg = "white")  # 调整图像的宽度和高度
###############

bb=aa[aa$sobel_p<0.05&aa$Proportion_Mediated_Estimate>0,]
cc=aa[aa$sobel_p<0.05/3665&aa$Proportion_Mediated_Estimate>0,]

# 绘制bb数据的直方图
p=ggplot(bb, aes(x = Proportion_Mediated_Estimate)) +
  # 绘制bb数据的直方图，设置宽度为 0.025，颜色为浅黄色，边框为灰色
  geom_histogram(binwidth = 0.0125, color = "grey50", fill = "#ffd92fFF", alpha = 0.5) +
  # 绘制cc数据的直方图，设置不同的颜色
  geom_histogram(data = cc, aes(x = Proportion_Mediated_Estimate), 
                 binwidth = 0.0125, color = "grey50", fill = "#66c2a5FF", alpha = 0.5) +
  # 添加红色加粗虚线，表示中位数
  geom_vline(aes(xintercept = median(Proportion_Mediated_Estimate, na.rm = TRUE)), 
             color = "blue", linetype = "dashed", size = 0.8) +
  geom_vline(aes(xintercept = 0.2516612), 
             color = "red", linetype = "dashed", size = 0.8) +
  # 设置主题样式
  theme_bw(base_size = 22) +
  # 自动根据数据设置 x 轴范围
  scale_x_continuous(expand = c(0.005, 0.005)) +
  theme(
    panel.grid = element_blank(),          # 去掉背景灰线
    panel.background = element_blank(),    # 去掉背景色
    legend.position = "bottom",
    axis.text = element_text(size = 18),   # 坐标轴字体大小
    #plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 标题样式
    #axis.title = element_text(face = "bold"),            # 坐标轴标题加粗
    text = element_text(size = 20)                      # 字体大小
  ) +
  # 添加标题和坐标轴标签
  labs(
    #title = "Mediation Proportion",
    x = "Mediation proportion",
    y = "Count"
  )
p
ggsave("fig5b.png", p,
       width = 8, height = 5.5, dpi = 300, bg = "white")  # 调整图像的宽度和高度


###画散点图#####


rm(list=ls())
aa=fread("相关性和偏相关性.txt")

library(ggplot2)
library(dplyr)
# 找到 corr.p < 0.05 的点的 corr 最大值和最小值
corr_p_subset <- aa %>% filter(corr.p > 0.05)
max_corr <- max(corr_p_subset$corr, na.rm = TRUE)
min_corr <- min(corr_p_subset$corr, na.rm = TRUE)

# 添加颜色分类，并根据 corr.p 调整透明度
aa <- aa %>%
  mutate(color_group = case_when(
    partial_corr.p > 0.05 ~ "> 0.05",
    partial_corr.p < 0.05/3665 ~ "< 0.05/3665",
    TRUE ~ "0.05/3665 - 0.05"
  ),
  alpha_value = ifelse(corr.p < 0.05, 1, 0.3))  # 设置透明度

# 绘制散点图
p <- ggplot(aa, aes(x = corr, y = partial_corr, color = color_group, alpha = alpha_value)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +  # 添加 y = x 蓝色虚线
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +  # x=0 的黑色实线
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +  # y=0 的黑色实线
  geom_vline(xintercept = max_corr, linetype = "dashed", color = "red") +  # 在最大值位置添加竖虚线
  geom_vline(xintercept = min_corr, linetype = "dashed", color = "red") +  # 在最小值位置添加竖虚线
  scale_color_manual(values = c("> 0.05" = "#ffd92fFF", 
                                "0.05/3665 - 0.05" = "#66c2a5FF", 
                                "< 0.05/3665" = "#fc8d62FF"),
                     breaks = c("> 0.05", "0.05/3665 - 0.05", "< 0.05/3665"),
                     labels = c(
                       "> 0.05" = expression("> 0.05"),
                       "0.05/3665 - 0.05" = expression(1.36 %*% 10^-5 ~ "- 0.05"),
                       "< 0.05/3665" = expression("<" ~ 1.36 %*% 10^-5)
                     )) +  # 指定图例顺序
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +  # 设置透明度范围，并隐藏 alpha 图例
  theme_minimal() +
  labs(x = "Correlation", y = "Partial correlation", color = expression("Partial correlation " ~ italic(P) ~ "-value")) +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 20), # 设置y轴标题大小为14，且加粗
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20), # 设置y轴标题大小为14，且加粗
        axis.text.x = element_text(size = 18))# 设置y轴刻度大小为12)  # 移除背景灰线
p
ggsave("fig5c.png", p,
       width = 9.6, height = 6, dpi = 300, bg = "white")  # 调整图像的宽度和高度



#######画韦恩图##################

library(VennDiagram)

ab$sobel_p[is.na(ab$sobel_p)] <- 1

library(dplyr)

ab <- ab %>%
  group_by(tag) %>%
  slice_min(sobel_p, with_ties = FALSE) 

grid.newpage()
venn.plot <- draw.triple.venn(
  area1 = nrow(ab[ab$sobel_p < (0.05)& ab$Proportion_Mediated_Estimate > 0,]),
  area2 = nrow(ab[ab$partial_corr.p < 0.05,]),
  area3 = nrow(ab[ab$corr.p < 0.05,]),
  n12 = nrow(ab[ab$sobel_p < 0.05 & ab$Proportion_Mediated_Estimate > 0& ab$partial_corr.p < 0.05,]),
  n13 = nrow(ab[ab$sobel_p < 0.05 & ab$Proportion_Mediated_Estimate > 0& ab$corr.p < 0.05,]),
  n23 = nrow(ab[ab$partial_corr.p < 0.05 & ab$corr.p < 0.05,]),
  n123 = nrow(ab[ab$sobel_p < 0.05 & ab$Proportion_Mediated_Estimate > 0& ab$partial_corr.p < 0.05 & ab$corr.p < 0.05,]),
  category = c(expression(atop("Sobel", italic(P) ~ " < 0.05")), 
               expression(atop("Partial corr", italic(P) ~ " < 0.05")),
               expression(atop("Original corr", italic(P) ~ " < 0.05"))),
  fill = c("#E7D4E8", "#A8DDB5", "#FEE440"),
  alpha = 0.5,
  cat.cex = 2.5,  # 调大 category 字体
  cex = 2.5,
  scaled = TRUE,
  cat.pos = c(310, 78, 170),  # 调整 category 角度，避免与数字重叠
  cat.dist = c(-0.1,-0.134,-0.150) , # 适当增加注释间距
  fontfamily     = rep("Arial", 7),  # 圆心数字的字体
  cat.fontfamily = rep("Arial", 3)   # category 标签的字体
)

png("fig5d.png",
    width = 10.5, height = 8.05, units = "in", res = 300)

grid.newpage()

# 添加全设备的黑色边框
grid.rect(
  gp = gpar(lwd = 4, col = "black"),
  width = unit(1, "npc"), 
  height = unit(1, "npc")
)

# 在设备中心推入一个正方形视口（8英寸×8英寸）
pushViewport(viewport(x = 0.5, y = 0.5, width = unit(10, "in"), height = unit(10, "in")))

# 绘制韦恩图，确保 scaled = FALSE 保证圆形（且此处坐标系为正方形）
grid.draw(venn.plot)
popViewport()

dev.off()






