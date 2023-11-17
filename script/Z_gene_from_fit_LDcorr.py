# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compute gene Z from SNP Z, weights from gamma distribution and correlation between SNP (LD)
# LD correlation computed from fitting in function of between distances in utils

import numpy as np
import scipy.stats as stats
import SNP_exon_param as P
import SNP_exon_utils as U


def read_sed():
    file = P.sed_file
    gene0_SNP = dict()
    gene1_SNP = dict()
    gene1_seW = dict()
    for line in open(file):
        sln = line.split()
        if sln[1] == '0':
            if sln[3] not in gene0_SNP:
                gene0_SNP[sln[3]] = set()
            gene0_SNP[sln[3]].add(sln[2])
        elif sln[1] in {'1', '2', '3'}:
            if sln[3] not in gene1_SNP:
                gene1_SNP[sln[3]] = list()
                gene1_seW[sln[3]] = list()
            gene1_SNP[sln[3]].append(sln[2])
            gene1_seW[sln[3]].append(U.weight_f(float(sln[5])))
    print('genes red from', P.sed_file)
    return gene0_SNP, gene1_SNP, gene1_seW


def read_SNP():
    SNP_Z = dict()
    SNP_pos = dict()
    file = P.pvZ_file
    for line in open(file):
        sln = line.split()
        SNP_Z[sln[0]] = float(sln[2])
        SNP_pos[sln[0]] = int(sln[4])
    print('SNPs red from', P.pvZ_file)
    return SNP_Z, SNP_pos


def print_mx(SNPs, R):
    sz = len(SNPs)
    print('.\t', end='')
    for i in range(sz):
        print(SNPs[i], end='\t')
    print()
    for i in range(sz):
        print(SNPs[i], end='\t')
        for j in range(sz):
            print(R[i, j], end='\t')
        print()


def LD_matrix(SNPs, SNP_pos):
    sz = len(SNPs)
    R = np.zeros((sz, sz))
    for i in range(1, sz):
        for j in range(i + 1, sz):
            R[i, j] = U.LD_corr(np.abs(SNP_pos[SNPs[i]] - SNP_pos[SNPs[j]]))
    return R


def compute_Z(W, Z, R):
    H = 0.0
    D = 0.0
    sz = len(Z)
    for i in range(0, sz):
        H = H + W[i] * Z[i]
        D = D + W[i] * W[i]
    for i in range(0, sz):
        for j in range(i + 1, sz):
            D = D + 2 * R[i, j] * W[i] * W[j]
    return H / np.sqrt(D)


def compute0_Z(gene, SNP_pos, SNP_Z, gene0_SNP):
    SNPs = list(gene0_SNP[gene])
    R = LD_matrix(SNPs, SNP_pos)
    W = np.full(len(SNPs), 1)
    Z = np.zeros(len(SNPs))
    for i in range(len(SNPs)):
        Z[i] = SNP_Z[SNPs[i]]
    return compute_Z(W, Z, R)


def compute1_Z(gene, SNP_pos, SNP_Z, gene1_SNP, gene1_seW):
    SNPs = gene1_SNP[gene]
    R = LD_matrix(SNPs, SNP_pos)
    W = np.array(gene1_seW[gene])
    Z = np.zeros(len(SNPs))
    for i in range(len(SNPs)):
        Z[i] = SNP_Z[SNPs[i]]
    return compute_Z(W, Z, R)


def save_Z_genes():
    SNP_Z, SNP_pos = read_SNP()
    gene0_SNP, gene1_SNP, gene1_seW = read_sed()
    g_in = open(P.gene_file)
    g_in.readline()  # title line
    line = g_in.readline()
    file = P.gene_Z_file_fitLD
    gz_out = open(file, mode='w')
    gz_out.write('Chr\tGene\tBPbeg\tBPend\tZ\tPV\tInt_Spl\n')
    while line:
        spl = line.split()
        Chr, BPbeg, BPend, Gene = spl[0], spl[1], spl[2], spl[3]
        print(Chr, Gene)
        if Gene in gene0_SNP:
            gZ = compute0_Z(Gene, SNP_pos, SNP_Z, gene0_SNP)
            PV = str(1 - stats.norm.cdf(gZ))
            ol = Chr + '\t' + Gene + '\t' + BPbeg + '\t' + BPend + '\t' + str(gZ) + '\t' + PV + '\t0'
            gz_out.write(ol + '\n')
        if Gene in gene1_SNP:
            gZ = compute1_Z(Gene, SNP_pos, SNP_Z, gene1_SNP, gene1_seW)
            PV = str(1 - stats.norm.cdf(gZ))
            ol = Chr + '\t' + Gene + '\t' + BPbeg + '\t' + BPend + '\t' + str(gZ) + '\t' + PV + '\t1'
            gz_out.write(ol + '\n')
        line = g_in.readline()
    gz_out.close()
    print("Z and pv genes written in", file)


if __name__ == "__main__":
    save_Z_genes()
