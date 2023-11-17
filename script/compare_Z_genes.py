# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# compare Z-score of genes from two files in result directory
#


import numpy as np
import matplotlib.pyplot as plt
import SNP_exon_param as P

################ sources for comparison in result dir ##############
# if format different from SNP exon with title line, input fields
# source: file, SNP exon format or not, gene_field, Z_field
source_1 = [P.gene_Z_file_LDcorr, True, -1, -1]
# source 2:
source_2 = [P.result_dir + 'ebi_006719.genes.out$.txt', False, 0, 7]
#####################################################################

def snp_exon_read_gene_Z(source):
    int_gz, spl_gz, gz = dict(), dict(), dict()
    g_in = open(source[0])
    g_in.readline()  # title line
    line = g_in.readline()
    while line:
        sl = line.split()
        if sl[6] == '0':
            int_gz[sl[1]] = float(sl[4])
        else:
            spl_gz[sl[1]] = float(sl[4])
        line = g_in.readline()
    print(source[0], 'red')
    genes = set(int_gz.keys()).union(spl_gz.keys())
    for g in genes:
        if g in int_gz and g in spl_gz:
            gz[g] = (int_gz[g] + spl_gz[g]) / np.sqrt(2)
        elif g in spl_gz:
            gz[g] = spl_gz[g]
        else:
            gz[g] = int_gz[g]
    return gz


def read_gene_Z(source):
    if source[1]:
        return snp_exon_read_gene_Z(source)
    else:
        gz = dict()
        g_in = open(source[0])
        g_in.readline()  # title line
        line = g_in.readline()
        while line:
            sl = line.split()
            gz[sl[source[2]]] = float(sl[source[3]])
            line = g_in.readline()
        print(source[0], 'red')
        return gz


def draw_comp(gz_1, gz_2):
    bins = 500
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95)
    axs[0, 0].set_title('Histogram Z-score from 1')
    axs[0, 0].hist(gz_1.values(), bins=bins, density=True)
    axs[0, 1].set_title('Histogram Z-score from 2')
    axs[0, 1].hist(gz_2.values(), bins=bins, density=True)
    axs[1, 0].set_title('Histogram of gap 2 - 1')
    axs[1, 0].hist([gz_2[g] - gz_1[g] for g in gz_1], bins=bins)
    zg_1 = [(gz_1[g], g) for g in gz_1]
    zg_1.sort()
    z_1 = [zg[0] for zg in zg_1]
    z_2 = [gz_2[zg[1]] for zg in zg_1]
    axs[1, 1].set_title('Sorted Z from 1 and Z from 2')
    axs[1, 1].plot(z_2)
    axs[1, 1].plot(z_1)


def common_genes_only(gz_1, gz_2):
    inter = set(gz_1.keys()).intersection(gz_2.keys())
    print('1 and 2 intersection / 1 size and 2 size:', len(inter), '/', len(gz_1), 'and', len(gz_2))
    cgz_1 = {key: value for key, value in gz_1.items() if key in inter}
    cgz_2 = {key: value for key, value in gz_2.items() if key in inter}
    return cgz_1, cgz_2


if __name__ == "__main__":
    gz_1 = read_gene_Z(source_1)
    gz_2 = read_gene_Z(source_2)
    cgz_1, cgz_2 = common_genes_only(gz_1, gz_2)
    draw_comp(cgz_1, cgz_2)
    plt.show()
