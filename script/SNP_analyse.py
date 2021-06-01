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
        self.SNP_pv = self.read_all_pv()

    def read_all_pv(self):
        SNP_pv = list()
        in_ = open(P.SNP_file)
        in_.readline()
        ln = in_.readline()
        while ln:
            sln = ln.split()
            SNP_pv.append(float(sln[P.SNP_pv]))
            ln = in_.readline()
        SNP_pv.sort()
        return SNP_pv

    def outliers(self):
        print('Choice of fork of p-values to eliminate outlier p-values')
        print('Once choosen, update the parameter to compute correlation coefficient')
        quant = [0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999]
        print('quantile\tvalue')
        for q in quant:
            print(q, np.quantile(self.SNP_pv, q), sep='\t')
        U.plot_histo(self.SNP_pv, 1000, title='Histogram of p-values')

    def fit_lin_log(self, inf_pv, sup_pv):
        print('\nSNP p-values')
        cx, cy = U.cumul_number(self.SNP_pv)
        cx, cy = cx[1:], cy[1:]
        i_min = np.where(cx > inf_pv)[0][0]
        i_max = np.where(cx < sup_pv)[0][-1]
        lx = np.log10(cx[i_min:i_max])
        ly = np.log10(cy[i_min:i_max])
        plt.figure()
        plt.title('linear adjustment of log p-value and log cumulative number')
        plt.plot(lx, ly)
        slope, intercept, r_value, p_value, std_err = stats.linregress(lx, ly)
        plt.plot(lx, lx * slope + intercept, label='R=' + '{:.4f}'.format(r_value))
        plt.legend()
        print('number of p-values', len(lx), 'out of a total of', len(self.SNP_pv))
        print('fork of p_values: [', cx[i_min], ', ', cx[i_max], ']', sep='')
        print('r_value:', r_value)
        print('p_value:', p_value)
        print('std_err:', std_err)

if __name__ == "__main__":
    sa = SNP_analyse()
    sa.outliers()
    sa.fit_lin_log(P.inf_pv, P.sup_pv)
    plt.show()
