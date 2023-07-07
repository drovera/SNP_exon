# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compute gene Z from SNP Z, weights from gamma distribution and correlation between SNP
# Correlation from LD correlation extracted by LDlinkR as triangular matrix

# format gene file
# 0: snp
# 1: Z
# 2: weight

# format correlation file as superior triangular matrix
# 0,0: snp; 1,0: corr; 2,0: corr;  ...
# ...
# 0,n-1: snp; 1, n-1: corr
# 0,n: snp

import os
import numpy as np
from scipy import stats
import SNP_exon_param as P


def get_Z_W(gene_file):
    snps, Z, W = list(), list(), list()
    for line in open(gene_file):
        line = line[0:-2]
        sln = line.split(';')
        snps.append(sln[0])
        Z.append(float(sln[1]))
        W.append(float(sln[2]))
    return snps, Z, W


def comp_gene_Z_corr(corr_file, gene_file):
    corr_dict = P.make_corr_dict(corr_file)
    snps, Z, W = get_Z_W(gene_file)
    n = len(snps)
    R = np.zeros([n, n])
    for r in range(n):
        for c in range(n):
            try:
                R[r, c] = corr_dict[snps[r]][snps[c]]
            except KeyError as e:
                R[r, c] = P.def_LDcorr
    N = np.dot(W, Z)
    D2 = np.matmul(np.matmul(W, R), np.transpose(W))
    gZ = N / np.sqrt(D2)
    return gZ


def comp_gene_Z_default(gene_file):
    snps, Z, W = get_Z_W(gene_file)
    R = np.empty((len(Z), len(Z)))
    for r in range(len(snps)):
        for c in range(len(snps)):
            if snps[r] == snps[c]:
                R[r, c] = 1.0
            else:
                R[r, c] = P.def_LDcorr
    return np.dot(W, Z) / np.sqrt(np.matmul(np.matmul(W, R), np.transpose(W)))


def compute_gene_Z(default=False):
    gene_Z = dict()
    for entry in os.scandir(P.ldcor_dir):
        if entry.path.endswith('_g.txt') and not entry.name.startswith('_'):
            gene_file = entry.path
            gene_pos = entry.name[0:-6]
            corr_file = entry.path[0:-6] + '_c.txt'
            if default:
                gZ = comp_gene_Z_default(gene_file)
                print(gene_pos, 'd', sep='\t')
            else:
                if os.path.exists(corr_file):
                    gZ = comp_gene_Z_corr(corr_file, gene_file)
                    print(gene_pos, 'c', sep='\t')
                else:
                    gZ = comp_gene_Z_default(gene_file)
                    print(gene_pos, 'd', sep='\t')
            gene_Z[gene_pos] = gZ
    return gene_Z


def save_Z_genes(default=False):
    gene_Z = compute_gene_Z(default=default)
    g_in = open(P.gene_file)
    g_in.readline()  # title line
    line = g_in.readline()
    if default:
        file = P.gene_Z_file_default
    else:
        file = P.gene_Z_file_ldcorr
    gz_out = open(file, mode='w')
    print('Compute PV abd save in', file)
    gz_out.write('Chr\tGene\tBPbeg\tBPend\tZ\tPV\tInt_Spl\n')
    while line:
        spl = line.split()
        Chr, BPbeg, BPend, Gene = spl[0], spl[1], spl[2], spl[3]
        gene_pos = Gene + '_0'  # case interior
        if gene_pos in gene_Z:
            gZ = gene_Z[gene_pos]
            PV = str(1 - stats.norm.cdf(gZ))
            ol = Chr + '\t' + Gene + '\t' + BPbeg + '\t' + BPend + '\t' + str(gZ) + '\t' + PV + '\t0'
            gz_out.write(ol + '\n')
        gene_pos = Gene + '_1'  # case splicing
        if gene_pos in gene_Z:
            gZ = gene_Z[gene_pos]
            PV = str(1 - stats.norm.cdf(gZ))
            ol = Chr + '\t' + Gene + '\t' + BPbeg + '\t' + BPend + '\t' + str(gZ) + '\t' + PV + '\t1'
            gz_out.write(ol + '\n')
        line = g_in.readline()
    gz_out.close()
    print("Z and pv genes written in", file)


if __name__ == "__main__":
    save_Z_genes(default=P.default_option)
