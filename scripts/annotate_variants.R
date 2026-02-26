# Load required packages
library(dplyr)
library(ggplot2)
library(readr)

# Read CSV of parsed variants
df <- read_csv('results/tables/pathogenic_variants.csv', show_col_types = FALSE)

# Summarize: number of unique pathogenic variants per gene
gene_counts <- df %>%
  filter(GENE != "" & !is.na(GENE)) %>%
  count(GENE, sort = TRUE)

# Save summary table
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(gene_counts, "results/tables/gene_pathogenic_counts.csv")

# Plot: Top 20 genes by pathogenic variant count
top_n <- 20
top_gene_counts <- gene_counts %>%
  top_n(top_n, n) %>%
  arrange(desc(n))

fig_dir <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
png(file.path(fig_dir, "top_genes_pathogenic_variant_counts.png"), width=900, height=500)
ggplot(top_gene_counts, aes(x = reorder(GENE, n), y = n)) +
  geom_bar(stat = 'identity', fill = "steelblue") +
  coord_flip() +
  labs(title = paste("Top", top_n, "Genes by Pathogenic Variant Count"),
       x = "Gene", y = "Pathogenic Variant Count") +
  theme_minimal(base_size = 14)
dev.off()

# Optional: Disease-gene association plot
gene_disease_counts <- df %>%
  filter(GENE != "", DISEASE != "") %>%
  count(GENE, DISEASE, sort = TRUE)

write_csv(gene_disease_counts, "results/tables/gene_disease_pathogenic_counts.csv")

# Plot: Top 20 disease-gene associations
top_gene_disease <- gene_disease_counts %>% top_n(top_n, n) %>% arrange(desc(n))

png(file.path(fig_dir, "top_gene_disease_associations.png"), width=900, height=500)
ggplot(top_gene_disease, aes(x = reorder(paste(GENE, DISEASE, sep=" - "), n), y = n)) +
  geom_bar(stat = 'identity', fill = "darkred") +
  coord_flip() +
  labs(title = paste("Top", top_n, "Gene-Disease Associations (Pathogenic Variants)"),
       x = "Gene - Disease", y = "Variant Count") +
  theme_minimal(base_size = 12)
dev.off()