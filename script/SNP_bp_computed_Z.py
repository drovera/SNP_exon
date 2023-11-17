# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compute Z-score of SNPs and save them with position of SNP

import math
from scipy import stats
import SNP_exon_param as P
import SNP_exon_utils as U

def SNP_compute_Z():
    seq_pos, SNP_pv = U.read_SNP(log=False)
    out = open(P.pvZ_file, mode='w')
    print('Infinity result of computing Z score')
    inf_pos, inf_neg = 0, 0
    for chr in range(24):
        for bp_snp in seq_pos[chr]:
            snp = bp_snp[2]
            bp = bp_snp[0]
            z = stats.norm.ppf(1 - SNP_pv[snp])
            if math.isinf(z):
                if z == math.inf:
                    inf_pos += 1
                else:
                    inf_neg += 1
                print(snp, z, SNP_pv[snp], sep='\t')
            else:
                out.write(snp + '\t' + str(SNP_pv[snp]) + '\t' + str(z) + '\t' + str(chr + 1) + '\t' + str(bp) + '\n')
    print('Z computed from p-values of SNP:', inf_pos, '+inf, ', inf_neg, '-inf /', len(SNP_pv))
    print('Z and position of SNPs saved in', P.pvZ_file)
    out.close()

if __name__ == "__main__":
    SNP_compute_Z()