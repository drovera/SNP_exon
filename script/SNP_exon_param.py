# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# common parameters: dir, files, values ...

### Parameters ###

# Common part of all directories to update
root = 'C:/Users/danie/Documents/SNP_exon/'
ref_dir = root + 'ref/'

# Extract by chromosome from NIH
# refSeq from NIH genome - https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens+%5Borgn%5D
extract_from_NIH = root + 'ref/NIHsequence/sequence_chr'

# Source of biological data about SNP
SNP_file = root + 'dta/icogs_bcac_public_results_euro.genesis.assoc.txt'
# Positions of used fields, File with header, first column = 0
title_line = True
# chromosome, name of SNP, base position, p-value
SNP_chr, SNP_name, SNP_bp, SNP_pv = 0, 1, 2, 8
# limit of p-value to eliminate outlier p-values in SNP_analyse:
inf_pv, sup_pv = 1.3805e-15, 0.99
# Correlation coefficient between SNP p-value
SNP_R = 0.8478

# Name of analysed data
data = 'BCAC'

# Parameters of gamma law min distance SNP - exons computed on SNP inside genes
shape, scale = 4.54746696e-01, 2.97027070e+04

# Choice of option for result
if_top = True
# Top number to display result
top = 10
# If not top, gene list file in input, use for detailed contribution of SNP
gene_list = 'VEGA50'
gene_list_file = root + 'result/GL_' + gene_list + '.txt'

### Other parameters and information ###

# Final result for gene as table
geneZ_file = root + 'result/' + data + '_Zgn.txt'
# Chr: chromosome
# Gene: gene, can be twice if internal and splicing effect
# BPbeg: begin of gene as base position
# BPend: end of gene as base position
# Z: computed Z-score of gene
# PV: p-value computed from Z-score
# InOut: 0 if internal effect, 1 if splicing effect

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

# Intermediate file for computing Z score
# 0: SNP_name
# 1: SNP Z score
pvZ_file = root + 'dta/' + data +'_pvZ.txt'

# Intermediate result for distance between SNP and exons, SNP exon distance
sed_file = root + 'dta/' + data + '_sed.txt'
# Every SNP is twice, even the records are identical
# fields for twice SNP
# 0: chromosome
# 1: relative position case
# 2: SNP
# 3: gene
# 4: exon order
# 5: distance SNP to exon

cases = ['SNP inside exon',
         'SNP inside the same gene',
         'SNP inside the the 2 nearest genes',
         'SNP inside this gene but outside the other nearest gene',
         'SNP outside this gene but inside the other nearest gene',
         'SNP between exons and outside the 2 nearest genes ',
         'SNP at the extremities and outside all genes']

# Chromosome size from NIH
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
156040895]




