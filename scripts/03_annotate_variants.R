# 03_analyze_variants.R — Comprehensive tidyverse analysis of ClinVar pathogenic variants.
#
# Reads:  results/tables/pathogenic_variants.csv  (produced by 02_parse_variants.py)
# Writes: results/tables/*.csv  and  results/figures/*.png
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
})
# ── Paths ─────────────────────────────────────────────────────────────────────
TABLE_DIR  <- "results/tables"
FIG_DIR    <- "results/figures"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)
# ── Load & clean data ─────────────────────────────────────────────────────────
df_raw <- read_csv(
  file.path(TABLE_DIR, "pathogenic_variants.csv"),
  col_types = cols(
    CHROM        = col_character(),
    POS          = col_double(),
    ID           = col_character(),
    REF          = col_character(),
    ALT          = col_character(),
    GENE         = col_character(),
    DISEASE      = col_character(),
    CLNSIG       = col_character(),
    VARIANT_TYPE = col_character()
  )
)
df <- df_raw %>%
  mutate(
    # Replace underscores with spaces in disease names
    DISEASE = str_replace_all(DISEASE, "_", " "),
    # Convert placeholder strings to NA
    DISEASE = na_if(DISEASE, "not provided"),
    DISEASE = na_if(DISEASE, "not specified"),
    GENE    = na_if(GENE, "")
  )
cat("\n════════════════════════════════════════════\n")
cat("  Dataset overview\n")
cat("════════════════════════════════════════════\n")
cat(sprintf("  Rows     : %d\n", nrow(df)))
cat(sprintf("  Columns  : %d\n", ncol(df)))
cat(sprintf("  Missing GENE    : %d (%.1f%%)\n",
            sum(is.na(df$GENE)),    mean(is.na(df$GENE))    * 100))
cat(sprintf("  Missing DISEASE : %d (%.1f%%)\n",
            sum(is.na(df$DISEASE)), mean(is.na(df$DISEASE)) * 100))
cat("════════════════════════════════════════════\n\n")
# ── Analysis 1: Gene-level pathogenic variant counts ─────────────────────────
cat("Running Analysis 1: Gene-level counts...\n")
gene_counts <- df %>%
  filter(!is.na(GENE)) %>%
  group_by(GENE) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
write_csv(gene_counts, file.path(TABLE_DIR, "gene_pathogenic_counts.csv"))
top20_genes <- gene_counts %>% slice_head(n = 20)
p1 <- ggplot(top20_genes,
             aes(x = fct_reorder(GENE, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  labs(
    title = "Top 20 Genes by Pathogenic Variant Count",
    x = "Gene",
    y = "Number of Pathogenic Variants"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(FIG_DIR, "top_genes_pathogenic_variant_counts.png"),
       p1, width = 10, height = 6, dpi = 300)
# ── Analysis 2: Gene–Disease association ─────────────────────────────────────
cat("Running Analysis 2: Gene–Disease associations...\n")
gene_disease_counts <- df %>%
  filter(!is.na(GENE), !is.na(DISEASE)) %>%
  group_by(GENE, DISEASE) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
write_csv(gene_disease_counts,
          file.path(TABLE_DIR, "gene_disease_pathogenic_counts.csv"))
top20_gd <- gene_disease_counts %>%
  slice_head(n = 20) %>%
  mutate(label = paste(GENE, DISEASE, sep = " — "))
p2 <- ggplot(top20_gd,
             aes(x = fct_reorder(label, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  labs(
    title = "Top 20 Gene–Disease Associations (Pathogenic Variants)",
    x = "Gene — Disease",
    y = "Variant Count"
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(FIG_DIR, "top_gene_disease_associations.png"),
       p2, width = 12, height = 7, dpi = 300)
# ── Analysis 3: Chromosome distribution ──────────────────────────────────────
cat("Running Analysis 3: Chromosome distribution...\n")
chrom_levels <- c(as.character(1:22), "X", "Y", "MT")
chrom_counts <- df %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
  group_by(CHROM) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(CHROM)
write_csv(chrom_counts, file.path(TABLE_DIR, "chromosome_distribution.csv"))
p3 <- ggplot(chrom_counts %>% filter(!is.na(CHROM)),
             aes(x = CHROM, y = n, fill = n)) +
  geom_col() +
  scale_fill_viridis_c(option = "viridis", guide = "none") +
  labs(
    title = "Pathogenic Variant Count by Chromosome",
    x = "Chromosome",
    y = "Number of Variants"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(FIG_DIR, "chromosome_distribution.png"),
       p3, width = 10, height = 5, dpi = 300)
# ── Analysis 4: Variant type breakdown ───────────────────────────────────────
cat("Running Analysis 4: Variant type breakdown...\n")
variant_type_summary <- df %>%
  group_by(VARIANT_TYPE) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(desc(count))
write_csv(variant_type_summary,
          file.path(TABLE_DIR, "variant_type_summary.csv"))
p4 <- ggplot(variant_type_summary,
             aes(x = fct_reorder(VARIANT_TYPE, count), y = count, fill = VARIANT_TYPE)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),
            hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_viridis_d(option = "cividis") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Variant Type Distribution",
    x = "Variant Type",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(FIG_DIR, "variant_type_distribution.png"),
       p4, width = 9, height = 5, dpi = 300)
# ── Analysis 5: Gene pleiotropy ───────────────────────────────────────────────
cat("Running Analysis 5: Gene pleiotropy...\n")
pleiotropy <- df %>%
  filter(!is.na(GENE), !is.na(DISEASE)) %>%
  group_by(GENE) %>%
  summarise(
    n_variants  = n(),
    n_diseases  = n_distinct(DISEASE),
    .groups = "drop"
  ) %>%
  filter(n_variants >= 5) %>%
  arrange(desc(n_diseases)) %>%
  slice_head(n = 30)
write_csv(pleiotropy, file.path(TABLE_DIR, "gene_pleiotropy.csv"))
top20_pleio <- pleiotropy %>% slice_head(n = 20)
p5 <- ggplot(top20_pleio,
             aes(x = fct_reorder(GENE, n_diseases), y = n_diseases, fill = n_diseases)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "inferno", guide = "none") +
  labs(
    title = "Top 20 Most Pleiotropic Genes",
    subtitle = "Genes associated with the most distinct diseases (≥5 pathogenic variants)",
    x = "Gene",
    y = "Number of Distinct Diseases"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(FIG_DIR, "gene_pleiotropy.png"),
       p5, width = 10, height = 6, dpi = 300)
# ── Analysis 6: Pipeline summary statistics ───────────────────────────────────
cat("Running Analysis 6: Summary statistics table...\n")
most_common_gene    <- gene_counts$GENE[1]
most_common_disease <- gene_disease_counts %>%
  group_by(DISEASE) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1) %>%
  pull(DISEASE)
variants_per_gene <- gene_counts$n
summary_stats <- tibble(
  statistic = c(
    "Total variants",
    "Unique genes",
    "Unique diseases",
    "Unique chromosomes",
    "Most common gene",
    "Most common disease",
    "Median variants per gene",
    "Mean variants per gene"
  ),
  value = c(
    as.character(nrow(df)),
    as.character(n_distinct(df$GENE, na.rm = TRUE)),
    as.character(n_distinct(df$DISEASE, na.rm = TRUE)),
    as.character(n_distinct(df$CHROM)),
    most_common_gene,
    most_common_disease,
    as.character(median(variants_per_gene)),
    sprintf("%.2f", mean(variants_per_gene))
  )
)
write_csv(summary_stats, file.path(TABLE_DIR, "pipeline_summary_stats.csv"))
cat("\n════════════════════════════════════════════\n")
cat("  Pipeline Summary Statistics\n")
cat("════════════════════════════════════════════\n")
for (i in seq_len(nrow(summary_stats))) {
  cat(sprintf("  %-30s %s\n",
              summary_stats$statistic[i],
              summary_stats$value[i]))
}
cat("════════════════════════════════════════════\n")
cat("\nAll analyses complete.\n")
cat(sprintf("Tables  → %s/\n", TABLE_DIR))
cat(sprintf("Figures → %s/\n", FIG_DIR))
