# Genomic Variant Analysis Project: Interview Summary

In this mini-project, I demonstrate my ability to integrate Unix shell scripting, Python, and R to analyze **real genomic variant data** from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/). The pipeline is robust, reproducible, and modular—suitable for larger-scale genomics projects and for clear review in interviews.

## Project Organization

I use a standard data science directory structure:
- `data/raw` – input VCF files
- `results/tables` – summary tables (CSV)
- `results/figures` – output plots
- `scripts` – all analysis scripts (Bash, Python, R)

## 1. Input Data

The primary data source is the **public ClinVar GRCh38 VCF file**, downloaded from NCBI. This real-world dataset contains both pathogenic and benign variants.

**How to get the data:**
```bash
# Download the latest ClinVar VCF (update URL if needed, see NCBI FTP)
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz
mkdir -p data/raw
mv clinvar.vcf data/raw/clinvar.vcf
```

## 2. Bash/Unix Scripting for Preprocessing

I use a Bash script to efficiently filter for pathogenic variants based on the `CLNSIG=Pathogenic` flag in the VCF INFO field:
```bash
grep '^#' data/raw/clinvar.vcf > data/raw/pathogenic_variants.vcf
grep 'CLNSIG=Pathogenic' data/raw/clinvar.vcf >> data/raw/pathogenic_variants.vcf
```
This produces a manageable VCF of only pathogenic variants, demonstrating high-throughput data triage.

## 3. Python for Data Parsing and Structuring

Using Python and pandas, I've parsed the filtered VCF to extract:
- Key columns: chromosome, position, ID, ref/alt alleles
- Fields from INFO: gene name (`GENEINFO`), disease (`CLNDN`), clinical significance (`CLNSIG`)
Output is a tidy CSV, ready for statistical and graphical analysis.
```bash
python scripts/parse_variants.py
```

## 4. R for Data Analysis and Visualization

In R:
- Summarize counts of pathogenic variants for each gene.
- Summarize top gene–disease associations.
- Create publication-quality plots with ggplot2.
- All output is saved to `results/tables/` and `results/figures/`.
```bash
Rscript scripts/annotate_variants.R
```

## 5. General Workflow & Best Practices

- **Modular:** Each step (Bash, Python, R) is decoupled and file-driven—easy to debug or extend.
- **Reproducible:** Full workflow can be rerun anytime datasets are updated.
- **Well-organized:** Directories are clean and ready for version control/code review.
- **Clearly documented:** This README serves as a step-by-step guide.

## **Skills & Concepts Demonstrated**

- Robust shell scripting for high-speed data filtering.
- Python for parsing semi-structured data (VCF) into analysis-ready tables.
- R for summarization and professional-quality visualization.
- File organization, pipeline modularity, automation, and bioinformatics literacy.
- Understanding and using public human disease variant data.

## **Example Project Walkthrough Script**

_"In this project, I designed and executed a reproducible, multi-language pipeline to analyze real genomic variant data from the ClinVar database. I organized the project using standard data science folder principles. Using Unix shell scripting, I filtered the ClinVar VCF file to isolate pathogenic variants. With Python, I parsed and exported key gene and disease information into tidy CSV tables. Finally, I used R to summarize and visualize gene–disease associations, producing figures of pathogenic variant counts per gene and top disease associations. This workflow demonstrates how I connect the strengths of Bash, Python, and R—going from raw data to publication-quality insight. The pipeline is extensible for larger datasets, more complex statistics, or future automation.”_

---

## **Reproducibility**

To reproduce this analysis:
1. Download the latest ClinVar VCF and place in `data/raw/`.
2. Run the filtering, parsing, and R scripts as above.
3. Review summary tables and output figures in `results/`.

---

## **Note**
- The large VCF file is **not included** in version control—please follow instructions above to obtain it for analysis.

---