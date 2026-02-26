import os
import csv

infile = 'data/raw/pathogenic_variants.vcf'
os.makedirs('results/tables', exist_ok=True)
outfile = 'results/tables/pathogenic_variants.csv'

with open(infile) as vcf, open(outfile, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    # Write CSV header
    writer.writerow(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENE', 'DISEASE', 'CLNSIG'])

    for line in vcf:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
        
        # Parse INFO: look for GENEINFO, CLNDN, CLNSIG
        gene, disease, clnsig = '', '', ''
        for entry in info.split(';'):
            if entry.startswith('GENEINFO='):
                # GENEINFO format: BRAF:673
                gene = entry.split('=')[1].split(':')[0] if ':' in entry else entry.split('=')[1]
            if entry.startswith('CLNDN='):
                disease = entry.split('=')[1]
            if entry.startswith('CLNSIG='):
                clnsig = entry.split('=')[1]
        writer.writerow([chrom, pos, vid, ref, alt, gene, disease, clnsig])