# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Draw manhattan plot of genes from results

import pandas as pd
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
        ax.set_xticks(self.pos)
        ax.set_xticklabels(self.chr)
        ax.set_ylabel('-log10(pvalue)')
        ax.set_title(title)

    def plot_df(self, ax, df, cols):
        for c in range(23):
            chr_df = df[(df['Chr'] == c + 1)]
            ax.scatter(chr_df['BP'], chr_df['-ln10(pv)'], c=cols[c % len(cols)], s=2)

    def plot(self, file):
        df = pd.read_table(file)
        df['BP'] = self.org[df['Chr'] - 1] + (df.BPbeg + df.BPend) // 2
        df['-ln10(pv)'] = -np.log10(df.PV)
        col0 = ['red', 'orange', 'yellow']
        col1 = ['green', 'blue', 'indigo']
        fig, axs = plt.subplots(2, 1)
        self.plot_df(axs[0], df[df['Int_Spl'] == 0], col0)
        self.axes(axs[0], 'P-value of genes resulting from effect of ' + P.data + '\nSNP inside exons')
        self.plot_df(axs[1], df[df['Int_Spl'] == 1], col1)
        self.axes(axs[1], 'SNP inside genes and outside exons')

    def threshold(self, ax):
        sig_thresh1 = 0.05
        sig_thresh2 = 0.01
        print("Significance threshold = ", sig_thresh1, "and", sig_thresh2)
        ax.plot([0, self.lim], [-np.log10(sig_thresh1), -np.log10(sig_thresh1)], lw=1, color='black')
        ax.plot([0, self.lim], [-np.log10(sig_thresh2), -np.log10(sig_thresh2)], lw=1, color='black')


if __name__ == "__main__":
    Manhattan_plot().plot(P.geneZ_file)
    plt.show()
