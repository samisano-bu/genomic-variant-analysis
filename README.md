# Genomic Variant Analysis Pipeline

A reproducible, three-stage bioinformatics pipeline that downloads the NCBI ClinVar GRCh38 VCF, filters for pathogenic variants, parses them into tidy tabular data, and performs six distinct analyses with publication-quality visualizations. The project is designed to demonstrate end-to-end data-engineering and analytical skills — Bash preprocessing, Python data wrangling with pandas, and R tidyverse statistical analysis — in the context of real clinical genomics data.

---

## Directory Structure

```
genomic-variant-analysis/
│
├── data/
│   └── raw/                    # ClinVar VCF download target (gitignored)
│
├── results/
│   ├── tables/                 # CSV output from Python & R steps
│   └── figures/                # PNG plots produced by R
│
├── scripts/
│   ├── 01_preprocess.sh        # Bash: download, decompress, filter VCF
│   ├── 02_parse_variants.py    # Python: parse VCF → tidy CSV
│   └── 03_analyze_variants.R   # R: six analyses + visualizations
│
├── install_r_packages.R        # Helper to install R dependencies
├── requirements.txt            # Python dependencies
├── .gitignore
└── README.md
```

---

## Prerequisites

| Tool | Version |
|------|---------|
| Bash | any modern shell |
| curl + gunzip | standard Unix utilities |
| Python | 3.8+ |
| R | 4.0+ |

## Installation

```bash
# Clone the repository
git clone https://github.com/samisano-bu/genomic-variant-analysis.git
cd genomic-variant-analysis

# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
Rscript install_r_packages.R
```

---

## Data Acquisition

The raw ClinVar GRCh38 VCF is downloaded automatically by `01_preprocess.sh`. If you prefer to download it manually first:

```bash
curl -L -o data/raw/clinvar.vcf.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip data/raw/clinvar.vcf.gz
```

---

## Pipeline Execution

Run the three scripts in order from the project root:

```bash
# Step 1 — Download & filter VCF (Bash)
bash scripts/01_preprocess.sh

# Step 2 — Parse VCF into tidy CSV (Python)
python scripts/02_parse_variants.py

# Step 3 — Analyze & visualize (R)
Rscript scripts/03_analyze_variants.R
```

---

## Output Files

### Tables (`results/tables/`)

| File | Description |
|------|-------------|
| `pathogenic_variants.csv` | All pathogenic variants with CHROM, POS, REF, ALT, GENE, DISEASE, CLNSIG, VARIANT_TYPE |
| `gene_pathogenic_counts.csv` | Pathogenic variant count per gene, sorted descending |
| `gene_disease_pathogenic_counts.csv` | Variant counts for every gene–disease pair |
| `chromosome_distribution.csv` | Variant counts per chromosome (ordered 1–22, X, Y, MT) |
| `variant_type_summary.csv` | Count and proportion of SNV / Insertion / Deletion / Complex |
| `gene_pleiotropy.csv` | Top 30 genes ranked by number of distinct associated diseases |
| `pipeline_summary_stats.csv` | High-level key statistics for the full dataset |

### Figures (`results/figures/`)

| File | Description |
|------|-------------|
| `top_genes_pathogenic_variant_counts.png` | Horizontal bar chart — top 20 genes by variant count |
| `top_gene_disease_associations.png` | Horizontal bar chart — top 20 gene–disease pairs |
| `chromosome_distribution.png` | Bar chart of variant counts by chromosome |
| `variant_type_distribution.png` | Bar chart of variant type proportions |
| `gene_pleiotropy.png` | Top 20 most pleiotropic genes |

---

## Analysis Highlights

1. **Gene-level counts** — identifies the most mutation-burdened genes in the ClinVar pathogenic set.
2. **Gene–Disease associations** — reveals which gene–disease pairs accumulate the highest number of curated pathogenic variants.
3. **Chromosome distribution** — shows genome-wide spread of pathogenic variants with correct cytogenetic ordering.
4. **Variant type breakdown** — classifies every variant as SNV, Insertion, Deletion, or Complex based on REF/ALT lengths.
5. **Gene pleiotropy** — ranks genes by the number of distinct diseases they are associated with, highlighting highly pleiotropic disease genes.
6. **Summary statistics** — a concise tidy table of pipeline-wide metrics suitable for reporting.

---

## Skills Demonstrated

**Bash**
- `set -euo pipefail` for safe scripting
- Conditional downloads with `curl` / `gunzip`
- VCF header-preserving filtering (`grep '^#'` + `grep 'CLNSIG=Pathogenic'`)
- Unix text-processing pipeline (`wc -l`, `cut`, `sort`, `uniq -c`, `awk`)

**Python**
- `pandas` DataFrame construction and export
- `re`-based INFO field parsing with a reusable `parse_info_field()` helper
- Variant classification logic (`classify_variant_type()`)
- QC summary reporting to stdout
- `if __name__ == '__main__':` guard and module docstring

**R / tidyverse**
- `readr` for type-safe CSV import
- `dplyr` for grouped summaries and filtering
- `stringr` for disease-name cleaning
- `forcats::fct_reorder()` for ordered categorical axes
- `ggplot2` + `scale_fill_viridis_c/d()` for publication-quality plots
- `ggsave()` at 300 DPI

---

## Contact

- Sean Amisano — sean.matteo.amisano@gmail.com
