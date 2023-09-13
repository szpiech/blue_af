# blue_af
usage: calc_blue_freq.py [-h] -v VCF_FILE -k KIN_MAT -i KING_IDS

Calculate BLUE of allele frequency. Can handle missing data. Outputs as a GARLIC .freq file.

options:
  -h, --help            show this help message and exit
  -v VCF_FILE, --vcf VCF_FILE
                        Gzipped VCF file
  -k KIN_MAT, --king KIN_MAT
                        Square kinship matrix. e.g. from plink --make-king square
  -i KING_IDS, --ids KING_IDS
                        List of individual IDs in order of kinship matrix. Should also be same order of samples in VCF. Assumes first row is a header.