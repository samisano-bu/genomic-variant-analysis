# Genomic Variant Analysis Pipeline

This mini-project demonstrates how to integrate Unix shell scripting, Python, and R to analyze real genomic variant data from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/). The workflow is robust, modular, and fully reproducible—ideal for bioinformatics review, interviews, and academic applications.

---

## Project Overview

*Objective:*  
To extract, parse, summarize, and visualize pathogenic variants from the ClinVar VCF dataset, identifying gene–disease associations and demonstrating practical skills at each stage.

*Skills Demonstrated:*
- Shell scripting for rapid, high-throughput data filtering
- Python/pandas for data wrangling and conversion to tidy tables
- R for statistical summarization and professional plotting
- Clean project organization and reproducibility
- Bioinformatics literacy (variant format, gene/disease annotation, summarization)

---

## Directory Structure

Your project is organized like a professional bioinformatics pipeline:

```
genomic-variant-analysis/
│
├── data/
│   └── raw/              # Original data: ClinVar VCF (not included in repo)
│
├── results/
│   ├── tables/           # Key CSV output files (variant tables, summary counts)
│   └── figures/          # Visualization PNGs made in R
│
├── scripts/              # Bash, Python, and R scripts for each workflow step
│
├── .gitignore            # Prevents committed data and system files
├── README.md             # You're reading it!
```

---

## Data Acquisition

**Raw Data:**  
The primary input is the [ClinVar GRCh38 VCF file](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz).  
Due to file size, raw data is NOT included in the repository.

**To obtain:**
```bash
# Download (update the URL if a newer release is available)
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz
mkdir -p data/raw
mv clinvar.vcf data/raw/clinvar.vcf
```

---

## Workflow Steps

### 1. Bash/Unix Preprocessing

Efficiently filter the VCF for only pathogenic variants:

```bash
grep '^#' data/raw/clinvar.vcf > data/raw/pathogenic_variants.vcf
grep 'CLNSIG=Pathogenic' data/raw/clinvar.vcf >> data/raw/pathogenic_variants.vcf
```

Result:  
- `data/raw/pathogenic_variants.vcf` — filtered variants for downstream analysis

---

### 2. Python Parsing & Structuring

Extract key information, convert VCF to a tidy CSV using pandas:

```bash
python scripts/parse_variants.py
```

Features:
- Parses columns: chromosome, position, ID, ref/alt alleles
- Extracts annotation: gene name (`GENEINFO`), disease (`CLNDN`), clinical significance (`CLNSIG`)
- Creates: `results/tables/pathogenic_variants.csv`

---

### 3. R Summarization & Visualization

Summarize, annotate, visualize using R and ggplot2:

```bash
Rscript scripts/annotate_variants.R
```

- Summarizes pathogenic variant counts per gene
- Highlights top gene–disease associations
- Outputs summary tables (`gene_pathogenic_counts.csv`, `gene_disease_pathogenic_counts.csv`)
- Saves publication-quality plots (PNG files) to `results/figures/`

---

## Quickstart / How to Reproduce

1. **Download raw ClinVar VCF:**  
   See above (bash commands).
2. **Filter pathogenic variants:**  
   ```bash
   grep '^#' data/raw/clinvar.vcf > data/raw/pathogenic_variants.vcf
   grep 'CLNSIG=Pathogenic' data/raw/clinvar.vcf >> data/raw/pathogenic_variants.vcf
   ```
3. **Parse with Python:**  
   ```bash
   python scripts/parse_variants.py
   ```
4. **Analyze and visualize with R:**  
   ```bash
   Rscript scripts/annotate_variants.R
   ```
5. **Review outputs:**  
   - Tables: `results/tables/`  
   - Figures: `results/figures/`  

*All steps are modular, file-driven, and can be debugged or extended as needed.*

---

## Expected Output Files

**`results/tables/`:**
- `pathogenic_variants.csv` — all pathogenic variants extracted
- `gene_pathogenic_counts.csv` — gene-level summary count
- `gene_disease_pathogenic_counts.csv` — gene–disease association summary

**`results/figures/`:**
- `top_genes_pathogenic_variant_counts.png`
- `top_gene_disease_associations.png`

---

## Example Project Walkthrough for Interview

> "This pipeline demonstrates my fluency bridging Unix, Python, and R to perform real-world genetics analysis. I start with raw ClinVar VCF data, filter for pathogenic variants using shell scripting, parse and structure the filtered output with pandas, and then conduct annotation, summarization, and visualization in R. Output tables and figures are organized for quick review, and the modular workflow is ready for scale, automation, or deeper statistical analysis. This is the kind of reproducible, robust pipeline used in modern genomics and bioinformatics."

---

## Best Practices

- No large or sensitive data in repo; all users must download raw input files.
- Scripts are modular and documented for each step.
- Directory structure ensures clarity and reproducibility.
- Results are versioned and easy to showcase.

---

## Contact

- Sean Amisano
- sean.matteo.amisano@gmail.com

---

## Notes

- For updates, new data releases, or expanded analyses, swap in newer ClinVar VCF files and re-run the pipeline.
- For custom annotation or expanded output, extend scripts in the `scripts/` directory.
- The filtering strategy can be adapted for variant significance, population frequency, or custom targets.

---

**Thank you for reviewing this project!**
*Designed to be fully reproducible, modular, and practical for both academic and industry settings.*
