# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Analyse LD corr to find defaut LD corr

import os
import numpy as np

app_dir = "C:/Users/gisel/Documents/SNP_exon_1/"
prj = "ebi_006719"
prj_dir = app_dir + prj
ld_dir = prj_dir + "/LD_corr/"
gz_file = prj_dir + "/result/" + prj + '_Zgn.txt'

snp_all_nb = 2378717  # must be entered to avoid counting
trgl_list = True


def stat_LD_corr_all():
    corrs = list()
    snps = set()
    for entry in os.scandir(ld_dir):
        if entry.path.endswith('_c.txt'):
            for line in open(entry.path):
                line = line[0:-2]
                sln = line.split(';')
                snps.add(sln[0])
                for i in range(1, len(sln)):
                    try:
                        corrs.append(np.sqrt(float(sln[i])))
                    except ValueError:
                        pass
    print('SNPs from LDlinkR / all SNPs', len(snps), snp_all_nb, len(snps) / snp_all_nb)
    print('Average of LD correlation:', np.average(corrs))
    hist, bins = np.histogram(corrs, bins=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    for i in range(hist.size):
        print('[', bins[i], bins[i + 1], ']', round(hist[i] / len(corrs), 3))


if __name__ == "__main__":
    stat_LD_corr_all()
