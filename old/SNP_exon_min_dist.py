# Draw 3 graphs
# -log10(p-value) in function of min distance SNP-exon
# idem by chromosome
# cumulative function of -log10(p-value) in function of min distance SNP-exon


# Description of exon_pos_file extracted from from NIH genome
# (https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens+%5Borgn%5D)
# 0: chromosome
# 1: base position
# 2: 1 for exon start, 2 for exon end
# 3: HUGO ID for Gene

import matplotlib.pyplot as plt
import math
from operator import itemgetter


# Name of data for graph title and Path of SNP file to type here
name_of_data = 'BCAC'
SNP_file = 'C:/Users/danie/Documents/Bio_Data/BCAC/icogs_bcac_public_results_euro.genesis.assoc.txt'
# File with header, input column number above, first column = 0
SNPchr = 0  # chromosome
SNPname = 1  # name of SNP
SNPbp = 2  # base position
SNPpv = 8  # p-value

exon_pos_file = 'exon_pos.txt'


def read_seq_pos(file):
    seq_pos = list()
    for line in open(file):
        sq = line.split()
        seq_pos.append((int(sq[0]), int(sq[1]), int(sq[2]), sq[3]))
    return seq_pos


def read_SNP(SNP_file):
    seq_pos = list()
    SNP_pv = dict()
    SNPi = open(SNP_file)
    SNPi.readline()  # title line
    SNPl = SNPi.readline()
    while SNPl:
        SNP = SNPl.split()
        pv = float(SNP[SNPpv])
        SNP_pv[SNP[SNPname]] = pv
        seq_pos.append((int(SNP[SNPchr]), int(SNP[SNPbp]), 0, SNP[SNPname]))
        SNPl = SNPi.readline()
    SNPi.close()
    return seq_pos, SNP_pv


def search_inf(seq_pos, bp):
    l, r = 0, len(seq_pos) - 1
    while l <= r:
        m = int(l + (r - l) / 2)
        if seq_pos[m][1] < bp:
            l = m + 1
        else:
            r = m - 1
    return r


def exon_SNP_nearest(exon_seq_pos, SNP_seq_pos):
    SNP_min_dist_exon = list()
    for chr in range(1, 24):
        exon_pos = list(filter(lambda x: x[0] == chr, exon_seq_pos))
        SNP_pos = list(filter(lambda x: x[0] == chr, SNP_seq_pos))
        for snpi in range(len(SNP_pos)):
            if SNP_pos[snpi][1] < exon_pos[0][1]:
                d_inf = exon_pos[0][1] - SNP_pos[snpi][1]
                d_min = d_inf
                gene = exon_pos[0][3]
            else:
                exon_i = search_inf(exon_pos, SNP_pos[snpi][1])
                d_inf = SNP_pos[snpi][1] - exon_pos[exon_i][1]
                if exon_i + 1 == len(exon_pos):
                    d_min = d_inf
                    gene = exon_pos[exon_i][3]
                else:
                    d_sup = exon_pos[exon_i + 1][1] - SNP_pos[snpi][1]
                    if d_inf < d_sup:
                        d_min = d_inf
                        gene = exon_pos[exon_i][3]
                    else:
                        d_min = d_sup
                        gene = exon_pos[exon_i + 1][3]
            SNP_min_dist_exon.append((chr, SNP_pos[snpi][1], d_min, SNP_pos[snpi][3], gene))
    print('exon_SNP_nearest done')
    SNP_min_dist_exon.sort(key=itemgetter(2))
    return SNP_min_dist_exon


def draw_dist_pv(SNP_min_dist_exon, SNP_pv):
    min_dist, log_pv, log_pv_cumul = list(), list(), list()
    chr_min_dist, chr_log_pv = list(), list()
    for _ in range(24):
        chr_min_dist.append(list())
        chr_log_pv.append(list())
    sum_log_pv = 0
    for smde in SNP_min_dist_exon:
        min_dist.append(smde[2])
        chr_min_dist[smde[0] - 1].append(smde[2])
        lpv = -math.log10(SNP_pv[smde[3]])
        log_pv.append(lpv)
        chr_log_pv[smde[0] - 1].append(lpv)
        sum_log_pv = sum_log_pv + lpv
        log_pv_cumul.append(sum_log_pv)
    print('\tmin', 'max', sep='\t')
    print('min distance', min(min_dist), max(min_dist), sep='\t')
    print('-log10(p-value)', min(log_pv), max(log_pv), sep='\t')
    plt.figure(1)
    g1 = plt.subplot()
    g1.scatter(min_dist, log_pv, marker='.', color='blue')
    g1.set_xlabel('min dist SNP-exon')
    g1.set_ylabel('-log10(p-value)')
    g1.set_title('SNP-exon on ' + name_of_data)
    plt.figure(2)
    g2 = plt.subplot()
    g2.set_xlabel('min dist SNP-exon')
    g2.set_ylabel('-log10(p-value) cumul')
    g2.set_title('SNP-exon on ' + name_of_data + ', cumulative -log10(p-value)')
    g2.plot(min_dist, log_pv_cumul)
    f, axs = plt.subplots(4, 6)
    chr = 0
    for i in range(4):
        for j in range(6):
            axs[i, j].scatter(chr_min_dist[chr], chr_log_pv[chr], marker='.', color='blue')
            axs[i, j].annotate('chr ' + str(chr + 1), xy=(0.9, 0.9), xycoords='axes fraction',
                               horizontalalignment='right', verticalalignment='top')
            axs[i, j].get_xaxis().set_visible(False)
            axs[i, j].get_yaxis().set_visible(False)
            chr += 1
    plt.show()


if __name__ == "__main__":
    SNP_seq_pos, SNP_pv = read_SNP(SNP_file)
    exon_seq_pos = read_seq_pos(exon_pos_file)
    SNP_min_dist_exon = exon_SNP_nearest(exon_seq_pos, SNP_seq_pos)
    draw_dist_pv(SNP_min_dist_exon, SNP_pv)
