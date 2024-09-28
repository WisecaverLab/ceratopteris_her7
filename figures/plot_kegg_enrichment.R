library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)

df = read.table("Desktop/functional_enrichment_results_KEGG_wminexp_lfc1.txt", sep = '\t', header = TRUE)
unique_kegg = read.table("Desktop/plotkegglist.txt", sep = '\t', header = TRUE)

df_filtered <- df %>%
  filter(kegg %in% unique_kegg$kegg)


p <- ggplot(df_filtered, aes(x = reorder(kegg, desc(order)), y = seqfreq, fill = adjpval < 0.05)) +
  geom_bar(stat = "identity", width = 0.5, color = "#1465AC") +
  scale_fill_manual(values = c("TRUE" = "#1465AC", "FALSE" = "white"), guide = FALSE) +  # Manual fill colors
  labs(x = "KEGG Pathway", y = "Seqfreq", fill = "Significance (adjpval < 0.05)") +  # Axis and legend labels
  ggtitle("Seqfreq for KEGG Pathways (Ordered by 'order' column)")  + # Plot title
  facet_grid(col = vars(set)) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 0.3, 0.6)) +
  theme_classic()




ggsave2('Desktop/test.pdf', device = "pdf", width = 20, height = 12, units = "cm")

 p
 
