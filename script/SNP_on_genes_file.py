# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Write result of computing Z-scores of genes from p-values of SNPs as file containing:
# Chr: chromosome
# Gene: gene name
# BPbeg: begin of gene in base positions
# BPend: end of gene in base positions
# Z: resulting Z-score
# PV:associated p-value
# Int_Spl: flag for the two cases, structural: 0, splicing: 1

from scipy import stats
import SNP_exon_param as P
import SNP_on_genes_utils as SG

def SNP_on_gene_file():
    sg = SG.SNP_on_gene_utils()
    SNP_Z = sg.read_SNP_Z()
    gene_SNP_in_w = sg.read_dist_in(SNP_Z)
    gene_Z_in, _ = sg.gene_Z(SNP_Z, gene_SNP_in_w, sg.denom_in(gene_SNP_in_w))
    gene_SNP_out_w = sg.read_dist_out(SNP_Z)
    gene_Z_out, _ = sg.gene_Z(SNP_Z, gene_SNP_out_w, sg.denom_out(gene_SNP_out_w))
    gene_Z_in = sg.tuple_list_to_dict(gene_Z_in)
    gene_Z_out = sg.tuple_list_to_dict(gene_Z_out)
    file = P.geneZ_file
    out = open(file, mode='w')
    out.write('Chr\tGene\tBPbeg\tBPend\tZ\tPV\tInt_Spl\n')
    for line in open(P.gene_file):
        sl = line.split()
        gene = sl[3]
        if gene in gene_Z_in:
            z = gene_Z_in[gene]
            ln = sl[0] + '\t' + gene + '\t' + sl[1] + '\t' + sl[2] + '\t'
            ln += str(z) + '\t' + str(1 - stats.norm.cdf(z)) + '\t0\n'
            out.write(ln)
        if gene in gene_Z_out:
            z = gene_Z_out[gene]
            ln = sl[0] + '\t' + gene + '\t' + sl[1] + '\t' + sl[2] + '\t'
            ln += str(z) + '\t' + str(1 - stats.norm.cdf(z)) + '\t1\n'
            out.write(ln)
    out.close()
    print('Result for genes written in', P.geneZ_file)

if __name__ == "__main__":
    SNP_on_gene_file()
