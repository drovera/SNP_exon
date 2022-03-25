# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Information to eliminate outlier p-values: histogram and quantiles
# Compute correlation coefficient between p-values

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import SNP_exon_param as P
import SNP_exon_utils as U


class SNP_analyse:


    def __init__(self):
        _, SNP_pv_dict = U.read_SNP(log=log)
        self.SNP_pv = list(SNP_pv_dict.values())

    def outliers(self):
        SNP_pv = self.SNP_pv
        print('Step 1- Choice of fork of p-values to eliminate outlier p-values')
        print('Quantiles are computed for p-values in', '[', P.inf_pv, ', ', P.sup_pv, ']')
        print('Once choosen, update the parameter to compute correlation coefficient at step 2')
        print('Warning: Python cannot compute the Z-score for p-value < 5.552e-17')
        quant = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999]
        print('quantile\tvalue')
        for q in quant:
            print(q, np.quantile(SNP_pv, q), sep='\t')
        U.plot_histo(SNP_pv, 1000, title=P.data + ' Histogram of p-values')

    def fit_lin_log(self, inf_pv, sup_pv):
        print('\nStep 2- Compute correlation coefficient between log of SNP p-values')
        print('Fixed limits of p-values', P.inf_pv, P.sup_pv)
        self.SNP_pv.sort()
        cx, cy = U.cumul_number(self.SNP_pv)
        cx, cy = cx[1:], cy[1:]
        i_min = np.where(cx > inf_pv)[0][0]
        i_max = np.where(cx < sup_pv)[0][-1]
        lx = np.log10(cx[i_min:i_max])
        ly = np.log10(cy[i_min:i_max])
        plt.figure()
        plt.title(P.data + ' linear adjustment of log p-value and log cumulative number')
        plt.xticks(ticks=[], labels=[])
        plt.yticks(ticks=[], labels=[])
        plt.plot(lx, ly)
        slope, intercept, r_value, p_value, std_err = stats.linregress(lx, ly)
        plt.plot(lx, lx * slope + intercept, label='R=' + '{:.4f}'.format(r_value))
        plt.legend()
        print('number of p-values', len(lx), 'out of a total of', len(self.SNP_pv))
        print('Real limits of p_values: [', cx[i_min], ', ', cx[i_max], ']', sep='')
        print('r_value:', r_value)
        print('p_value:', p_value)
        print('std_err:', std_err)

if __name__ == "__main__":
    # log = True to get details of excluded p-values
    log = False
    sa = SNP_analyse()
    sa.outliers()
    sa.fit_lin_log(P.inf_pv, P.sup_pv)
    plt.show()
