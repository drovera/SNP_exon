# Complete process for genes
# Generate the p-value of genes after fixing parameters of the data in SNP_exon_param

import SNP_exon_param as P
import SNP_compute_Z
import SNP_exon_distance
import SNP_on_genes_file

if __name__ == "__main__":
    SNP_compute_Z.SNP_compute_Z()
    SNP_exon_distance.SNP_exon_distance().proceed()
    SNP_on_genes_file.SNP_on_gene_file()
    print('process completed, 3 files created:')
    print(P.data + '_pvZ.txt, ' + P.data + '_sed.txt and ' + P.data + '_Zgn.txt')