# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Write result of computing Z-scores of genes from p-values of SNPs
# For top Z-scores of genes
# of
# For a list of genes

# SNP_gene_Z dict of dict, SNP_gene_Z[snp][gene]=z
# gene_Z: list(z, gene)

import sys
from scipy import stats
import SNP_exon_param as P
import SNP_on_genes_utils as SG


class SNP_on_gene_analyse:

    def __init__(self):
        self.sg = SG.SNP_on_gene_utils()

    def gene_SNP_top(self, top):
        gene_Z_in, SNP_gene_Z_in = self.sg.gene_Z(self.sg.gene_SNP_in_w,  self.sg.denom_in(self.sg.gene_SNP_in_w))
        gene_Z_out, SNP_gene_Z_out = self.sg.gene_Z(self.sg.gene_SNP_out_w, self.sg.denom_out(self.sg.gene_SNP_out_w))
        result_file = P.root + 'result/' + P.data + '_T' + str(top) + '.txt'
        out = open(result_file, mode='w')
        out.write(str(top) + 'top genes by internal effect\n')
        top_gene = [gene_Z_in[i][1] for i in range(top)]
        top_gene_Z = [gene_Z_in[i][0] for i in range(top)]
        self.write_detail_array(out, top_gene, top_gene_Z, SNP_gene_Z_in)
        out.write(str(top) + 'top genes by splicing effect\n')
        top_gene = [gene_Z_out[i][1] for i in range(top)]
        top_gene_Z = [gene_Z_out[i][0] for i in range(top)]
        self.write_detail_array(out, top_gene, top_gene_Z, SNP_gene_Z_out)
        out.write(str(top) + 'top genes by added Z-scores, internal effect in first and splicing effect in second\n')
        gene_Z_in_d = self.sg.tuple_list_to_dict(gene_Z_in)
        gene_Z_out_d = self.sg.tuple_list_to_dict(gene_Z_out)
        gene_Z_io = list()
        for sgi in gene_Z_in:
            gene = sgi[1]
            if gene in gene_Z_out_d:
                gene_Z_io.append((sgi[0] + gene_Z_out_d[gene], gene))
        gene_Z_io.sort(reverse=True)
        top_gene = [gene_Z_io[i][1] for i in range(top)]
        out.write('1- only internal effect\n')
        top_gene_Z = [gene_Z_in_d[gene] for gene in top_gene]
        self.write_detail_array(out, top_gene, top_gene_Z, SNP_gene_Z_in)
        out.write('2- only splicing effect\n')
        top_gene_Z = [gene_Z_out_d[gene] for gene in top_gene]
        self.write_detail_array(out, top_gene, top_gene_Z, SNP_gene_Z_out)
        out.close()
        print('Result of top', str(top), 'in', result_file)

    def gene_SNP_list(self):
        gene_Z_in, SNP_gene_Z_in = self.sg.gene_Z(self.sg.gene_SNP_in_w,  self.sg.denom_in(self.sg.gene_SNP_in_w))
        gene_Z_out, SNP_gene_Z_out = self.sg.gene_Z(self.sg.gene_SNP_out_w, self.sg.denom_out(self.sg.gene_SNP_out_w))
        gene_Z_in = self.sg.tuple_list_to_dict(gene_Z_in)
        gene_Z_out = self.sg.tuple_list_to_dict(gene_Z_out)
        gene_list_in, gene_list_Z_in = list(), list()
        gene_list_out, gene_list_Z_out = list(), list()
        list_file = P.root + 'result/GL_' + P.gene_list + '.txt'
        result_file = P.root + 'result/' + P.data + '_' + P.gene_list + '.txt'
        out = open(result_file, mode='w')
        print('List of genes red from', list_file)
        for gene in open(list_file):
            gene = gene.replace('\r', '').replace('\n', '')
            if gene not in gene_Z_in and gene not in gene_Z_out:
                out.write(gene + '\tnot found effect\n')
            else:
                if gene in gene_Z_in:
                    gene_list_in.append(gene)
                    gene_list_Z_in.append(gene_Z_in[gene])
                if gene in gene_Z_out:
                    gene_list_out.append(gene)
                    gene_list_Z_out.append(gene_Z_out[gene])
        out.write('1- only internal effect\n')
        self.write_detail_array(out, gene_list_in , gene_list_Z_in, SNP_gene_Z_in)
        out.write('2- only splicing effect\n')
        self.write_detail_array(out, gene_list_out, gene_list_Z_out, SNP_gene_Z_out)
        out.close()
        print('Result of', P.gene_list, 'in', result_file)

    def write_detail_array(self, out, gene_list, top_gene_Z, SNP_gene_Z):
        top_SNP = set()
        for snp in SNP_gene_Z:
            for gene in SNP_gene_Z[snp]:
                if gene in gene_list:
                    top_SNP.add(snp)
        top_SNP = list(top_SNP)
        top_SNP.sort(key=lambda x: self.sg.SNP_Z[x], reverse=True)
        form1 = '{:.4f}'
        out.write('R=' + str(P.SNP_R))
        for gene in gene_list:
            out.write('\t' + gene)
        out.write('\ngene Z')
        for Z in top_gene_Z:
            out.write('\t' + form1.format(Z))
        out.write('\tSNP_Z\n')
        for snp in top_SNP:
            out.write(snp)
            for gene in gene_list:
                if gene in SNP_gene_Z[snp]:
                    out.write('\t' + form1.format(SNP_gene_Z[snp][gene]))
                else:
                    out.write('\t0')
            out.write('\t' + form1.format(self.sg.SNP_Z[snp]) + '\n')

if __name__ == "__main__":
    sog = SNP_on_gene_analyse()
    if P.if_top:
        sog.gene_SNP_top(P.top)
    else:
        sog.gene_SNP_list()



