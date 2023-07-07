# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Draw graph Z-score of SNP in function of min distance SNP exons
# Graw graph to help choosing a threshold of Z-score to limit analysis of SNP effect

# Z_dist: list of (Z, min distance to exon, SNP)

import SNP_exon_param as P
import SNP_exon_utils as U
import numpy as np
import matplotlib.pyplot as plt

class SNP_exon_dist_Z:

    def read_distance(self):
        ''' read only SNP with effect case 0, 1, 2 3'''
        min_dist_btw = dict()
        file = P.sed_file
        in_ = open(file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            ln = in_.readline()
            sln2 = ln.split()
            if int(sln1[1]) == 0:
                min_dist_btw[sln1[2]] = 0.0
            elif int(sln1[1]) in {1, 2}:
                min_dist_btw[sln1[2]] = min(int(sln1[5]), int(sln2[5]))
            elif int(sln1[1]) == 3:
                min_dist_btw[sln1[2]] = int(sln1[5])
            elif int(sln2[1]) == 3:
                min_dist_btw[sln2[2]] = int(sln2[5])
        print('distance SNP exons red from', file)
        print('SNP inside genes:', len(min_dist_btw))
        return min_dist_btw

    def read_distance_by_gene(self, gene_list):
        ''' distance SNPs inside genes to exons '''
        gene_set = set(gene_list)
        gene_snp, gene_d, gene_z = dict(), dict(), dict()
        for gene in gene_set:
            gene_snp[gene], gene_d[gene], gene_z[gene] = list(), list(), list()
        snp_z = U.read_SNP_Z()
        file = P.sed_file
        in_ = open(file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            c1, snp, gene1, d1 = int(sln1[1]), sln1[2], sln1[3], int(sln1[5])
            ln = in_.readline()
            sln2 = ln.split()
            c2, gene2, d2 = int(sln2[1]), sln2[3], int(sln2[5])
            if c1 in {0, 1}:
                if gene1 in gene_set:
                    gene_snp[gene1].append(snp) # gene1 == gene 2
                    gene_d[gene1].append(min(d1, d2))
                    gene_z[gene1].append(snp_z[snp])
            elif c1 == 2: # between 2 exons of different genes
                if gene1 in gene_set:
                    gene_snp[gene1].append(snp)
                    gene_d[gene1].append(d1)
                    gene_z[gene1].append(snp_z[snp])
                if gene2 in gene_set:
                    gene_snp[gene2].append(snp)
                    gene_d[gene2].append(d2)
                    gene_z[gene2].append(snp_z[snp])
            else: # inside one gene only
                if c1 == 3:
                    if gene1 in gene_set:
                        gene_snp[gene1].append(snp)
                        gene_d[gene1].append(d1)
                        gene_z[gene1].append(snp_z[snp])
                if c2 == 3:
                    if gene2 in gene_set:
                        gene_snp[gene2].append(snp)
                        gene_d[gene2].append(d2)
                        gene_z[gene2].append(snp_z[snp])
        print('distance SNP exons red from', file)
        return gene_snp, gene_d, gene_z

    def read_SNP_Z_dist(self):
        ''' combine Z and distance '''
        min_dist_btw = self.read_distance()
        file = P.pvZ_file
        Z_dist = list()
        for ln in open(file):
            sln = ln.split()
            if sln[0] in min_dist_btw:
                Z_dist.append((float(sln[2]), min_dist_btw[sln[0]], sln[0]))
        print('Z-score of SNP red from', file)
        return Z_dist

    def draw_d_z(self, Z_dist):
        size = len(Z_dist)
        d, z, i = np.zeros(size, dtype=int), np.zeros(size, dtype=float), 0
        for zd in Z_dist:
            d[i], z[i] = zd[1], zd[0]
            i += 1
        plt.scatter(d, z, s=1)
        plt.title('Z-score, min distance SNP exons for ' + P.source)
        plt.xlabel('min distance')
        plt.ylabel('Z-score')
        plt.axhline(c='black')

    def highest_d_z(self, Z_dist, part):
        ''' draw extreme values histogram'''
        md, mz, nb, mx = list(), list(), list(), 0
        n, mn = 0, int(len(Z_dist) * part)
        for z_d in Z_dist:
            mz.append(z_d[0])
            if z_d[1] > mx:
                mx = z_d[1]
            md.append(mx)
            n += 1
            nb.append(n)
            if n > mn:
                break
        fig, ax = plt.subplots(1, 2)
        fig.suptitle('Part ' + str(part) + ' of all SNP for ' + P.source)
        ax[0].set_xlabel('max of min distance')
        ax[0].set_ylabel('Z-score')
        ax[0].plot(md, mz)
        ax[1].set_xlabel('number')
        ax[1].plot(nb, mz)


    def draw_all_and_top(self, part_of_SNP):
        Z_dist = self.read_SNP_Z_dist()
        self.draw_d_z(Z_dist)
        Z_dist.sort(reverse=True)
        for part in part_of_SNP:
            self.highest_d_z(Z_dist, part)


if __name__ == "__main__":
    sed = SNP_exon_dist_Z()
    part_of_SNP = [0.01, 0.001, 0.0001]
    sed.draw_all_and_top(part_of_SNP)
    plt.show()