"""
02_parse_variants.py — Parse a pathogenic ClinVar VCF into a tidy CSV.
Reads:  data/raw/pathogenic_variants.vcf
Writes: results/tables/pathogenic_variants.csv
Columns emitted:
    CHROM, POS, ID, REF, ALT, GENE, DISEASE, CLNSIG, VARIANT_TYPE
"""
import os
import re
import pandas as pd
# ── Helper functions ──────────────────────────────────────────────────────────
def parse_info_field(info_str: str, key: str) -> str:
    """Return the value for *key* from a VCF INFO string, or '' if absent."""
    match = re.search(rf'(?:^|;){re.escape(key)}=([^;]+)', info_str)
    return match.group(1) if match else ''
def classify_variant_type(ref: str, alt: str) -> str:
    """Classify a variant as SNV, Insertion, Deletion, or Complex."""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    if len(ref) < len(alt):
        return 'Insertion'
    if len(ref) > len(alt):
        return 'Deletion'
    return 'Complex'
def parse_vcf(vcf_path: str) -> pd.DataFrame:
    """Parse a VCF file and return a tidy pandas DataFrame."""
    rows = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info = fields[:8]
            gene_raw = parse_info_field(info, 'GENEINFO')
            # GENEINFO format: SYMBOL:EntrezID  — keep only the symbol
            gene = gene_raw.split(':')[0] if gene_raw else ''
            disease = parse_info_field(info, 'CLNDN')
            clnsig = parse_info_field(info, 'CLNSIG')
            variant_type = classify_variant_type(ref, alt)
            rows.append({
                'CHROM': chrom,
                'POS': int(pos),
                'ID': vid,
                'REF': ref,
                'ALT': alt,
                'GENE': gene,
                'DISEASE': disease,
                'CLNSIG': clnsig,
                'VARIANT_TYPE': variant_type,
            })
    return pd.DataFrame(rows)
# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    VCF_PATH = 'data/raw/pathogenic_variants.vcf'
    OUT_PATH = 'results/tables/pathogenic_variants.csv'
    os.makedirs('results/tables', exist_ok=True)
    print(f"Parsing {VCF_PATH} …")
    df = parse_vcf(VCF_PATH)
    # ── QC summary ────────────────────────────────────────────────────────────
    total = len(df)
    unique_genes = df['GENE'].replace('', pd.NA).dropna().nunique()
    missing_gene = df['GENE'].replace('', pd.NA).isna().sum()
    missing_pct = missing_gene / total * 100 if total > 0 else 0.0
    chroms = sorted(df['CHROM'].unique(), key=lambda c: (c.isdigit() is False, c))
    print("\n════════════════════════════════════════════")
    print("  QC Summary")
    print("════════════════════════════════════════════")
    print(f"  Total variants         : {total:,}")
    print(f"  Unique genes           : {unique_genes:,}")
    print(f"  Missing gene annotation: {missing_gene:,} ({missing_pct:.1f}%)")
    print(f"  Chromosomes represented: {', '.join(chroms)}")
    print("\n  Top 10 genes by variant count:")
    top_genes = (
        df['GENE'].replace('', pd.NA)
        .dropna()
        .value_counts()
        .head(10)
    )
    for gene, cnt in top_genes.items():
        print(f"    {gene:<12} {cnt:>6,}")
    print("\n  Variant type breakdown:")
    for vtype, cnt in df['VARIANT_TYPE'].value_counts().items():
        print(f"    {vtype:<12} {cnt:>6,} ({cnt/total*100:.1f}%)")
    print("════════════════════════════════════════════\n")
    df.to_csv(OUT_PATH, index=False)
    print(f"Saved {total:,} variants to {OUT_PATH}")
