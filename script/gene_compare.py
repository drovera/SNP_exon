# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compare several lists of gene by Z-score, by rank and by Z

import math
import SNP_exon_param as P
import numpy as np

data_srcs = [3, 4, 6]  # number of GWAS data corresponding to P.data_scr

class gene_compare:

    def __init__(self):
        self.nd = len(data_srcs)
        self.empty = np.empty(self.nd + 2, dtype=float)
        self.empty.fill(None)


    def add_gene(self, data_scr, gene, z, gene_z):
        if gene not in gene_z:
            gene_z[gene] = np.array(self.empty)
        gene_z[gene][data_scr] = z

    def complete(self, gene_z):
        for gn in gene_z:
            count = 0
            z = 0.0
            for ds in range(self.nd):
                if not math.isnan(gene_z[gn][ds]):
                    count += 1
                    z += gene_z[gn][ds]
            gene_z[gn][self.nd] = count
            gene_z[gn][self.nd + 1] = z / math.sqrt(count)

    def read_all_Zgn(self):
        gene_0_z = dict()
        gene_1_z = dict()
        for ds in range(self.nd):
            file = P.result_dir + P.datas[data_srcs[ds]] + '_Zgn.txt'
            in_ = open(file)
            in_.readline()
            ln = in_.readline()
            while ln:
                sln = ln.split()
                in_out = sln[6]
                if in_out == '0':
                    self.add_gene(ds, sln[1], float(sln[4]), gene_0_z)
                if in_out == '1':
                    self.add_gene(ds, sln[1], float(sln[4]), gene_1_z)
                ln = in_.readline()
            print('genes red from', file)
        self.complete(gene_0_z)
        self.complete(gene_1_z)
        return gene_0_z, gene_1_z


    def save_gene_z(self, gene_z, marker):
        file = P.result_dir + 'meta_' + marker + '.txt'
        out = open(file, mode='w')
        ln = 'gene\t'
        for ds in range(self.nd):
            ln += P.datas[data_srcs[ds]] + '\t'
        ln += 'nb\tresultZ\n'
        out.write(ln)
        for gn in gene_z:
            ln = gn + '\t'
            for ds in range(self.nd):
                ln += str(gene_z[gn][ds]) + '\t'
            ln += str(int(gene_z[gn][self.nd])) + '\t' + str(gene_z[gn][self.nd + 1]) + '\n'
            out.write(ln)
        out.close()
        print('Result of comparison in', file)

if __name__ == "__main__":
    gc = gene_compare()
    gene_0_z, gene_1_z = gc.read_all_Zgn()
    marker = str(gc.nd) + 'int'  # default value for in exon, direct effect
    gc.save_gene_z(gene_0_z, marker)
    marker = str(gc.nd) + 'spl'  # default value for in gene and out exon, splicing effect
    gc.save_gene_z(gene_1_z, marker)
