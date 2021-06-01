# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compute Z-scores of genes from Z-scores of SNP and weights

# SNP_Z: dict of SNP giving Z

# Gene_SNP_in or Gene_SNP_out, list of list
# 0: gene
# 1: list of
#   1.0: SNP
#   1.1: weight

import math
import copy as c
import SNP_exon_param as P
import SNP_exon_utils as U

class SNP_on_gene_utils:

    def __init__(self):
        self.SNP_Z = dict()
        self.gene_SNP_in_w = dict()
        self.gene_SNP_out_w = dict()
        self.read_SNP_Z()
        self.read_dist_in()
        self.read_dist_out()


    def read_SNP_Z(self):
        for ln in open(P.pvZ_file):
            sln = ln.split()
            self.SNP_Z[sln[0]] = float(sln[2])

    def read_dist_in(self):
        file = P.sed_file
        in_ = open(file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            in_.readline()
            if sln1[2] in self.SNP_Z:
                if int(sln1[1]) == 0:
                    if sln1[3] not in self.gene_SNP_in_w:
                        self.gene_SNP_in_w[sln1[3]] = list()
                    self.gene_SNP_in_w[sln1[3]].append((sln1[2], 1.0))
        print('Distances SNP to exons, SNP inside exons, from', file)

    def read_dist_out(self):
        file = P.sed_file
        in_ = open(file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            ln = in_.readline()
            sln2 = ln.split()
            if int(sln1[1]) in {1, 2}:
                if sln1[3] not in self.gene_SNP_out_w:
                    self.gene_SNP_out_w[sln1[3]] = list()
                if sln2[3] not in self.gene_SNP_out_w:
                    self.gene_SNP_out_w[sln2[3]] = list()
                if sln1[2] in self.SNP_Z:
                    self.gene_SNP_out_w[sln1[3]].append((sln1[2], U.weight_f(int(sln1[5]))))
                if sln2[2] in self.SNP_Z:
                    self.gene_SNP_out_w[sln2[3]].append((sln2[2], U.weight_f(int(sln2[5]))))
            if int(sln1[1]) == 3:
                if sln1[3] not in self.gene_SNP_out_w:
                    self.gene_SNP_out_w[sln1[3]] = list()
                if sln1[2] in self.SNP_Z:
                    self.gene_SNP_out_w[sln1[3]].append((sln1[2], U.weight_f(int(sln1[5]))))
            if int(sln2[1]) == 3:
                if sln2[3] not in self.gene_SNP_out_w:
                    self.gene_SNP_out_w[sln2[3]] = list()
                if sln2[2] in self.SNP_Z:
                    self.gene_SNP_out_w[sln2[3]].append((sln2[2], U.weight_f(int(sln2[5]))))
        print('Distances SNP to exons, SNP outside exons and inside genes, from', file)

    def denom_in(self, gene_SNP_w):
        gene_SNP_d = dict()
        for gene in gene_SNP_w:
            n = len(gene_SNP_w[gene])
            gene_SNP_d[gene] = math.sqrt(n + P.SNP_R * n * (n - 1))
        return gene_SNP_d

    def denom_out(self, gene_SNP_w):
        gene_SNP_d = dict()
        for gene in gene_SNP_w:
            denom = 0.0
            it1 = iter(gene_SNP_w[gene])
            while True:
                try:
                    snp1 = next(it1)
                    denom = denom + snp1[1] * snp1[1]
                    it2 = c.copy(it1)
                    while True:
                        try:
                            snp2 = next(it2)
                            denom = denom + P.SNP_R * snp1[1] * snp2[1]
                        except StopIteration:
                            break
                except StopIteration:
                    break
            gene_SNP_d[gene] = math.sqrt(denom)
        return gene_SNP_d

    def gene_Z(self, gene_SNP_w, gene_SNP_d):
        gene_Z = list()
        SNP_gene_Z = dict()
        for snp in self.SNP_Z:
            SNP_gene_Z[snp] = dict()
        for gene in gene_SNP_w:
            for sw in gene_SNP_w[gene]:
                SNP_gene_Z[sw[0]][gene] = 0.0
        for gene in gene_SNP_w:
            gz = 0.0
            it = iter(gene_SNP_w[gene])
            while True:
                try:
                    snp = next(it)
                    sz = snp[1] * self.SNP_Z[snp[0]] / gene_SNP_d[gene]
                    SNP_gene_Z[snp[0]][gene] = sz + SNP_gene_Z[snp[0]][gene]
                    gz = gz + sz
                except StopIteration:
                    break
            gene_Z.append((gz, gene))
        gene_Z.sort(reverse=True)
        return gene_Z, SNP_gene_Z

    def tuple_list_to_dict(self, list2):
        dict_ = dict()
        key, val = 1, 0
        for item in list2:
            dict_[item[key]] = item[val]
        return dict_

