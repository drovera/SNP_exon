# Compute the contribution of gene_Z of highest Z SNP
# The choice of SNP Z threshold and max of min distance SNP exon is done by SNP_exon_dist_Z
# The input value is the number of highest Z SNP see SNP_exon_param

# SNP_w_gene: list of tuple (SNP, weight, SNP Z, gene)
# result : list of tuple (SNP, SNP Z, effect, gene, gene Z)

import SNP_exon_param as P
import SNP_on_genes_utils as SG

log_dir = 'C:/Users/danie/Documents/data/algo/'

class SNP_Z_on_gene:

    def write_log(self, gene_SNP_out_w):
        out = open(log_dir + P.datas[P.data_src] + '_gsw.log', mode='w')
        for gene in gene_SNP_out_w:
            for sw in gene_SNP_out_w[gene]:
                txt = gene + '\t' + sw[0] + '\t' + str(sw[1])
                out.write(txt + '\n')
        out.close()
        print('gene_SNP_out_w written')

    def read_Zgn(self):
        file = P.geneZ_file
        in_ = open(file)
        gene_Z = dict()
        eof = False
        in_.readline()  # title line
        ln = in_.readline()
        while ln:
            sln = ln.split()
            if sln[6] == '1':
                gene_Z[sln[1]] = sln[4]
            ln = in_.readline()
        print(file, 'is red')
        return gene_Z

    def proceed(self):
        sg = SG.SNP_on_gene_utils()
        SNP_Z = sg.read_SNP_Z()
        gene_SNP_out_w = sg.read_dist_out(SNP_Z)
        gene_Z = self.read_Zgn()
        gene_denom = sg.denom_out(gene_SNP_out_w)
        SNP_w_gene = list()
        for gene in gene_SNP_out_w:
            for sw in gene_SNP_out_w[gene]:
                if SNP_Z[sw[0]] > P.Z_SNP_thr:
                    SNP_w_gene.append((sw[0], SNP_Z[sw[0]], sw[1], gene))
        result = list()
        it = iter(SNP_w_gene)
        prev = next(it)
        sum_ = prev[1] * prev[2]
        while True:
            try:
                item = next(it)
                if prev[0] == item[0] and prev[3] == item[3]:
                    sum_ += item[1] * item[2]
                else:
                    result.append((prev[0], SNP_Z[prev[0]], sum_ / gene_denom[prev[3]], prev[3], gene_Z[prev[3]]))
                    #print(prev[0], sum_, prev[3], sep='\t')
                    prev = item
                    sum_ = prev[1] * prev[2]
            except StopIteration:
                result.append((prev[0], SNP_Z[prev[0]], sum_ / gene_denom[prev[3]], prev[3], gene_Z[prev[3]]))
                #print(prev[0], sum_, prev[3], sep='\t')
                break
        print('SNP\tSNP_Z\tZ_to_gene\tgene\tgene_Z')
        for r in result:
            print(r[0], r[1], r[2], r[3], r[4], sep='\t')


if __name__ == "__main__":
    sog = SNP_Z_on_gene()
    sog.proceed()