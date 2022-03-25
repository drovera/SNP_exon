# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Compute Z-score of SNPs

import math
from scipy import stats
import SNP_exon_param as P
import SNP_exon_utils as U

def SNP_compute_Z():
    # log = True to get details of excluded p-values
    log = False
    _, SNP_pv = U.read_SNP(log)
    out = open(P.pvZ_file, mode='w')
    print('Infinity result of computing Z score')
    inf_pos, inf_neg = 0, 0
    for snp in SNP_pv:
        z = stats.norm.ppf(1 - SNP_pv[snp])
        if math.isinf(z):
            if z == math.inf:
                inf_pos += 1
            else:
                inf_neg += 1
            print(snp, z, SNP_pv[snp], sep='\t')
        else:
            out.write(snp + '\t' + str(SNP_pv[snp]) + '\t' + str(z) + '\n')
    print('Z computed from p-values of SNP:', inf_pos, '+inf, ', inf_neg, '-inf /', len(SNP_pv))
    out.close()

if __name__ == "__main__":
    SNP_compute_Z()