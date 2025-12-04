
library(data.table)
library(ggplot2)

A = fread(file = 'msea_ora_result.csv')

ggplot(data = head(A, 16), 
       mapping = aes(x = Erichment_Ratio, 
                     y = reorder(KEGG, -P),
                     colour = P, 
                     size = Hits)) + 
  geom_point() + 
  xlab(label = 'Enrichment Ratio') +
  ylab("KEGG") +
  theme(axis.text.y = element_text(size = 14))

ggsave(filename = 'meta_enrichment.tiff',
       width = 14, height = 8, units = 'in', dpi = 300)
