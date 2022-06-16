# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Draw graph Z-score of SNP in function of min distance SNP exons
# Graw graph to help choosing a threshold of Z-score to limit analysis of SNP effect
#

# Z_dist: list of (Z, min distance to exon, SNP)

import SNP_exon_param as P
import SNP_on_genes_utils as SG
import numpy as np
import matplotlib.pyplot as plt

class SNP_exon_dist_Z:

    def read_distance(self):
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
        gene_set = set(gene_list)
        gene_snp, gene_d, gene_z = dict(), dict(), dict()
        for gene in gene_set:
            gene_snp[gene], gene_d[gene], gene_z[gene] = list(), list(), list()
        snp_z = SG.SNP_on_gene_utils().read_SNP_Z()
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
        plt.title('Z-score, min distance SNP exons for ' + P.data)
        plt.xlabel('min distance')
        plt.ylabel('Z-score')
        plt.axhline(c='black')

    def highest_d_z(self, Z_dist, part):
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
        fig.suptitle('Part ' + str(part) + ' of all SNP for ' + P.data)
        ax[0].set_xlabel('max of min distance')
        ax[0].set_ylabel('Z-score')
        ax[0].plot(md, mz)
        ax[1].set_xlabel('number')
        ax[1].plot(nb, mz)

    def print_top(self, top, Z_dist):
        it = iter(Z_dist)
        print('SNP\tZ\tmin_dist')
        for _ in range(top):
            z_d = next(it)
            print(z_d[2], z_d[0], z_d[1], sep='\t')

    def draw_all_and_top(self, part_of_SNP):
        Z_dist = sed.read_SNP_Z_dist()
        sed.draw_d_z(Z_dist)
        Z_dist.sort(reverse=True)
        sed.print_top(P.top, Z_dist)
        for part in part_of_SNP:
            sed.highest_d_z(Z_dist, part)

    def draw_by_gene_from_list(self):
        print('List of genes red from', P.gene_list_file)
        gene_list = list()
        for gene in open(P.gene_list_file):
            gene = gene.replace('\r', '').replace('\n', '')
            gene_list.append(gene)
        gene_snp, gene_d, gene_z = self.read_distance_by_gene(gene_list)
        print('gene\tSNP\tmin_dist\tSNP_Z')
        for gene in gene_list:
            snp_i, d_i, z_i = iter(gene_snp[gene]), iter(gene_d[gene]), iter(gene_z[gene])
            while True:
                try:
                    print(gene, next(snp_i), next(d_i), next(z_i), sep='\t')
                except StopIteration:
                    break
        maxd, minz, maxz = 0, 0.0, 0.0
        for gene in gene_list:
            if gene_d[gene]:
                maxd = max(maxd, max(gene_d[gene]))
            if gene_z[gene]:
                minz = min(minz, min(gene_z[gene]))
                maxz = max(maxz, max(gene_z[gene]))
        lims = ', dist SNP-exon in [0, {0}], Z in [{1:.3f}, {2:.3f}]'.format(maxd, minz, maxz)
        ig, n = iter(gene_list), 0
        while True:
            try:
                if not n % 12:
                    fig, ax = plt.subplots(3, 4)
                    fig.suptitle(P.data + lims)
                    for i in range(3):
                        for j in range(4):
                            ax[i, j].tick_params(left=False, right=False, labelleft=False,
                                                 labelbottom=False, bottom=False)
                    for i in range(3):
                        for j in range(4):
                            gene = next(ig)
                            ax[i, j].scatter(gene_d[gene], gene_z[gene], s=1)
                            ax[i, j].set_xlim([0, maxd])
                            ax[i, j].set_ylim([minz, maxz])
                            ax[i, j].axhline(c='blue', lw=1)
                            ax[i, j].set_title(gene)
            except StopIteration:
                break

if __name__ == "__main__":
    sed = SNP_exon_dist_Z()
    if P.top > 0:
        part_of_SNP = [0.01, 0.001, 0.0001]
        sed.draw_all_and_top(part_of_SNP)
    else:
        try:
            sed.draw_by_gene_from_list()
        except FileNotFoundError:
            print('File Error, check name in parameters or directory')
    plt.show()