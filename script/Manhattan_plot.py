# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Draw manhattan plot of genes from results of SNP_exons

# Fields of Z-score gene file
# 0: Chr: chromosome
# 1: Gene: gene influenced by SNPs
# 2: BPbeg: begin base position
# 3: BPend: end base position
# 4: Z: Z-score of gene
# 5: PV: p-value of gene
# 6: Int\_Spl: 0 or 1, 0: SNP inside exon, 1: SNP inside gene influencing splicing

import os
import numpy as np
import SNP_exon_param as P
import matplotlib.pyplot as plt


class Manhattan_plot:

    def __init__(self):
        self.org = np.zeros(24)
        self.pos = np.zeros(24)
        self.chr = list()
        self.lim = P.chr_size[0]
        for c in range(22):
            self.org[c + 1] = self.org[c] + P.chr_size[c]
            self.pos[c] = self.org[c] + P.chr_size[c] // 2
            self.chr.append(str(c + 1))
            self.lim = self.lim + P.chr_size[c + 1]
        self.chr.append('23')
        self.pos[22] = self.org[22] + P.chr_size[22] // 2

    def axes(self, ax, title):
        ax.set_xlim([0, self.lim])
        ax.set_ylim(bottom=0)
        ax.set_xticks(self.pos[0:-1])
        ax.set_xticklabels(self.chr)
        ax.set_ylabel('-log10(pvalue)')
        ax.set_title(title)


    def plot_gz(self, gz, ax, cols, title):
        size = len(gz)
        chr = np.zeros(size, dtype=int)
        for i in range(size):
            chr[i] = int(gz[i][0]) - 1
        bp = np.zeros(size, dtype=np.uint)
        lpv = np.zeros((24, size))
        for i in range(size):
            bp[i] = self.org[chr[i]] + (int(gz[i][2]) + int(gz[i][3])) // 2
            lpv[chr[i]][i] = -np.log10(float(gz[i][5]))
        for c in range(23):
            ax.scatter(bp, lpv[c], c=cols[c % len(cols)], s=2)
        self.axes(ax, title)
        self.threshold(ax)


    def plot(self, file):
        int_gz, spl_gz = list(), list()
        g_in = open(file)
        g_in.readline()  # title line
        line = g_in.readline()
        while line:
            spl = line.split()
            if spl[6] == '0':
                int_gz.append((spl[0], spl[1], spl[2], spl[3], spl[4], spl[5]))
            else:  #  spl[6] == '1'
                spl_gz.append((spl[0], spl[1], spl[2], spl[3], spl[4], spl[5]))
            line = g_in.readline()
        fig, axs = plt.subplots(2, 1)
        title = 'P-value of genes resulting from effect of ' + P.source + '\nSNP inside exons'
        self.plot_gz(int_gz, axs[0], ['red', 'orange', 'yellow'], title)
        title = 'SNP inside genes and outside exons'
        self.plot_gz(spl_gz, axs[1], ['green', 'blue', 'indigo'], title)

    def threshold(self, ax):
        sig_thresh1 = 0.05
        sig_thresh2 = 0.01
        print("Significance threshold = ", sig_thresh1, "and", sig_thresh2)
        ax.plot([0, self.lim], [-np.log10(sig_thresh1), -np.log10(sig_thresh1)], lw=1, color='black')
        ax.plot([0, self.lim], [-np.log10(sig_thresh2), -np.log10(sig_thresh2)], lw=1, color='black')

if __name__ == "__main__":
    fig1 = False
    if os.path.isfile(P.gene_Z_file_fitLD):
        fig1 = True
        print('fig1: ', end='')
        print('Manhattan plot for LD correlation coefficient computed by fitting')
        Manhattan_plot().plot(P.gene_Z_file_fitLD)
    if os.path.isfile(P.gene_Z_file_LDcorr):
        if fig1:
            print('fig2: ', end='')
        print('Manhattan plot for LD correlation coefficient extracted from repository')
        Manhattan_plot().plot(P.gene_Z_file_LDcorr)
    plt.show()