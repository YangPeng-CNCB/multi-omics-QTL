
#############画图##########################

rm(list=ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

gtex=fread("coloc/用于画图的数据/画category图.txt")


# 转换格式
result <- gtex %>%
  pivot_longer(cols = c(unique, shared),
               names_to = "gene_type", values_to = "count") %>%
  mutate(gene_type = ifelse(gene_type == "unique", "Unique", "Shared")) %>%
  rename(Trait = phenotype) 



color_map <- data.frame(
  Trait = rep(c("CpG", "Protein", "Metabolite"), each = 2),
  gene_type = rep(c("Shared", "Unique"), times = 3),
  fill_color = c("#ffd92fFF", adjustcolor("#ffd92fFF", alpha.f = 0.4),
                 "#fc8d62FF", adjustcolor("#fc8d62FF", alpha.f = 0.4),
                 "#66c2a5FF", adjustcolor("#66c2a5FF", alpha.f = 0.4))
)
# 
result <- merge(result, color_map, by = c("Trait", "gene_type"), all.x = TRUE)

# 
result <- result %>%
  mutate(Trait = factor(Trait, levels = c("CpG", "Protein", "Metabolite")))




########

p= ggplot(result, aes(x = Trait, y = count*100, fill = gene_type)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = fill_color)) +
  facet_wrap(~category, scales = "free_x", nrow = 1, strip.position = "bottom") +
  scale_fill_identity() + 
  theme_minimal() +
  labs(
    x = NULL,
    y = "Percentage of colocalized pairs (%)"
  ) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.line.x = element_blank(),  
    strip.text = element_text(face = "plain", size = 14), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    plot.subtitle = element_text(size = 14, face = "bold"), 
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    axis.line.y = element_line(color = "black"),
    axis.title.y = element_text(size = 14), 
    axis.text.y = element_text(size = 12) 
  )
p

ggsave("coloc/用于画图的数据/fig4d.png", p,
       width = 8, height = 5.5, dpi = 300, bg = "white")  








