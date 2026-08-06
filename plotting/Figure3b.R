library(data.table)
rm(list=ls())

loci=fread("loci.txt")
loci=loci[,c(17:19)]
loci=unique(loci)
head(loci)

library(dplyr)
library(ggplot2)

# 
unique_tags_count <- loci %>%
  group_by(unique_tags) %>%
  summarise(count = n())

# 
unique_tags_count$unique_tags <- factor(unique_tags_count$unique_tags, 
                                        levels = c(as.character(sort(as.numeric(unique_tags_count$unique_tags[unique_tags_count$unique_tags != ">=10"]))), ">=10"))



# 
print(unique_tags_count)

# 
p=ggplot(unique_tags_count, aes(x = unique_tags, y = count, fill = "blue2")) +
  geom_bar(stat = "identity") +  
  geom_text(aes(label = count), vjust = -0.3) +
  labs(#title = "The number of metabolites associated with each locus",
       x = "Number of associated metabolites per locus", y = "Number of loci") +  
  theme(axis.text.x = element_text(size = 12, hjust = 0.5),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        panel.background = element_rect(fill = "white", color = NA),  
        panel.grid = element_blank(),  
        axis.line = element_line(color = "black"),  
        plot.title = element_text(hjust = 0.5),  
        legend.position = "none")
p
ggsave("fig3b.png", p, width = 8, height = 5, dpi = 300, bg = "white")  


