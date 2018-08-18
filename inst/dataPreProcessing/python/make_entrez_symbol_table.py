import pandas as pd

# download this:
# ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# make sure that this file is in same path
entrez = pd.read_csv("Homo_sapiens.gene_info.gz", sep="\t", usecols=[1,2], index_col=0)

entrez.to_csv("Entrez_symbol.tsv", sep="\t")

