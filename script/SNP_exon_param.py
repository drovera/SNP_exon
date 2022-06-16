# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# common parameters: dir, files, sources, intermediate, results ...

#########################################

# Structure of directories to update #
root = 'C:/Users/danie/Documents/SNP_exon/'
ref_dir = root + 'ref/'
dta_dir = root + 'dta/'
result_dir = root + 'result/'

#################################################################

# GWAS sources and their features, choice the source by GWAS_src #

# Name of analysed data: data
# Positions of used fields, File with header, first column = 0
# chromosome: SNP_chr, name of SNP: SNP_name , base position: SNP_bp, p-value: SNP_pv
# inf_pv, sup_pv are limits of p-value to eliminate eventually outlier p-values
# start values are inf_pv, sup_pv = 0.0, 1.0
# Correlation coefficient between SNP p-value: SNP_R
### update inf_pv, sup_pv and SNP_R by SNP_analyse ###
# Parameters shape and scale, gamma law min distance SNP - exons computed on SNPs inside genes
### update shape and scale by SNP_exon_analyse.py ###
datas = ['ebi_006719']      # 0
         # http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006719

data_src = 0
if data_src == 0:
    data = datas[data_src]
    SNP_file = dta_dir + 'ebi_006719_download.txt'
    title_line = True
    SNP_chr, SNP_name, SNP_bp, SNP_pv = 0, 1, 2, 6
    inf_pv, sup_pv = 0, 1
    SNP_R = 0.9994
    shape = 0.4538
    scale = 34094.3
else:
    data = 'unknown'

##############################################################################

# Choice of option for result and result files #

# SNP_on_genes_top.py and SNP_exon_dist_Z.py:
# Top number to display result in SNP_on_genes_top.py for gene and in SNP_exon_dist_Z.py for SNP
# 0 if computing from a list of genes, limit the number of genes used for computing detailed results
top = 0
# if top equals to 0, gene list file in input for detailed contribution of SNP
# as array by SNP_on_genes_top.py
# as graph by SNP_exon_dist_Z.py
gene_list = 'ebi_006719_top12'
gene_list_file = result_dir + gene_list + '.txt'

# Result for gene as table
geneZ_file = result_dir + data + '_Zgn.txt'
# Chr: chromosome  # 0
# Gene: gene, can be twice if internal and splicing effect   # 1
# BPbeg: begin of gene as base position  # 2
# BPend: end of gene as base position  # 3
# Z: computed Z-score of gene  # 4
# PV: p-value computed from Z-score   # 5
# InOut: 0 if internal effect, 1 if splicing effect  # 6

# SNP_Z_on_gene.py: the SNP Z threshold, choice by SNP_exon_dist_Z
Z_SNP_thr = 4.0

######################################################

# Intermediate files

# Intermediate file for computing Z score
# 0: SNP_name
# 1: SNP p-value
# 2: SNP Z score
pvZ_file = dta_dir + data +'_pvZ.txt'

# Intermediate result for distance between SNP and exons, SNP exon distance
sed_file = dta_dir + data + '_sed.txt'
# Every SNP is twice, even the records are identical
# fields for twice SNP
# 0: chromosome
# 1: relative position case
# 2: SNP
# 3: gene
# 4: exon order
# 5: distance SNP to exon

# Labels of relative position cases
cases = ['SNP inside exon',  # 0
         'SNP inside the same gene',  # 1
         'SNP inside the the 2 nearest genes',  # 2
         'SNP inside this gene but outside the other nearest gene',  # 3
         'SNP outside this gene but inside the other nearest gene',  # 4
         'SNP between exons and outside the 2 nearest genes ',  # 5
         'SNP at the extremities and outside all genes']  # 6

########################################################################################

# Position of exons and genes from NIH #

# Extraction by chromosome from NIH: parameter and output
### launch in first exons_from_NIH_nh after extracting sequence files ###
# refSeq from NIH genome - https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens+%5Borgn%5D
extract_from_NIH = ref_dir + 'NIHsequence/sequence_chr'

# Description of exon_pos.txt
# 0: chromosome
# 1: begin base position
# 2: end base position
# 3: NIH name of gene
# 4: exon order
exon_file = ref_dir + 'exon_pos.txt'

# Description of gene_pos.txt
# 0: chromosome
# 1: begin base position
# 2: end base position
# 3: NIH name
# 4: exon number
gene_file = ref_dir + 'gene_pos.txt'

# Abnormality of extraction are listed in
NIH_abn = ref_dir + 'exon_abn.txt'

# Chromosome size
chr_size = [
248956422,
242193529,
198295559,
190214555,
181538259,
170805979,
159345973,
145138636,
138394717,
133797422,
135086622,
133275309,
114364328,
107043718,
101991189,
90338345,
83257441,
80373285,
58617616,
64444167,
46709983,
50818468,
156040895,
57227415]