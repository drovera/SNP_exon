# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Analyse of relative positions of SNP to exons
# Fit gamma distribution and give parameters shape and scale
# Visualize this analyse with graphs


from scipy import stats
import matplotlib.pyplot as plt
import SNP_exon_utils as U
import SNP_exon_param as P


class SNP_exon_analyse:

    def read_distance(self):
        min_dist_all = list()
        min_dist_btw = list()
        min_dist_out = list()
        in_ = open(P.sed_file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            ln = in_.readline()
            sln2 = ln.split()
            if int(sln1[1]) in {1, 2, 3, 4, 5, 6}:
                min_dist_all.append(min(int(sln1[5]), int(sln2[5])))
            if int(sln1[1]) in {1, 2}:
                min_dist_btw.append(min(int(sln1[5]), int(sln2[5])))
            if int(sln1[1]) == 3:
                min_dist_btw.append(int(sln1[5]))
            if int(sln2[1]) == 3:
                min_dist_btw.append(int(sln2[5]))
            if int(sln1[1]) == 4:
                min_dist_out.append(int(sln1[5]))
            if int(sln2[1]) == 4:
                min_dist_out.append(int(sln2[5]))
            if int(sln1[1]) == 5:
                min_dist_out.append(min(int(sln1[5]), int(sln2[5])))
            if int(sln2[1]) == 6:
                min_dist_out.append(int(sln2[5]))
        print('Minimal distance only concerns SNP outside exons')
        print('all:', len(min_dist_all), 'inside genes:', len(min_dist_btw), 'outside genes:', len(min_dist_out), 'SNP')
        min_dist_all.sort()
        min_dist_btw.sort()
        min_dist_out.sort()
        return min_dist_all, min_dist_btw, min_dist_out

    def fit_draw(self):
        min_dist_all, min_dist_btw, min_dist_out = self.read_distance()
        print()
        cxa, cya = U.cumul_number(min_dist_all)
        cxb, cyb = U.cumul_number(min_dist_btw)
        cxo, cyo = U.cumul_number(min_dist_out)
        plt.title(P.data + ' minimal distance between SNP and exon en x, normalized cumulative number en y')
        plt.xscale('log')
        plt.plot(cxa, cya, label='all')
        plt.plot(cxb, cyb, label='in genes')
        plt.plot(cxo, cyo, label='out genes')
        plt.legend()
        print('Fitting for all SNPs')
        param = U.curve_fit(cxa, cya, U.gamma_cdf)
        title = P.data + ' min distance from all SNPs to exons / normalized cumulative number\n'
        title = title + 'Gamma law, shape=' + '{:.5f}'.format(param[0]) + ', scale=' + '{:.0f}'.format(param[1])
        std = U.gamma_draw_fit(cxa, cya, U.gamma_cdf, param, title=title)
        print('shape\tscale\tstd')
        print(param[0], param[1], std, sep='\t')
        print('\nFitting for SNPs inside genes')
        param = U.curve_fit(cxb, cyb, U.gamma_cdf)
        title = P.data + ' min distance from SNPs inside genes to exons / normalized cumulative number\n'
        title = title + 'Gamma law, shape=' + '{:.5f}'.format(param[0]) + ', scale=' + '{:.0f}'.format(param[1])
        std = U.gamma_draw_fit(cxb, cyb, U.gamma_cdf, param, title=title)
        print('shape\tscale\tstd')
        print(param[0], param[1], std, sep='\t')
        bin_nb = 1000
        U.plot_histo(cxb, bin_nb, P.data + ' minimal distance from SNPs inside gene to exons', noyticks=True)
        ax = plt.figure().add_subplot()
        stats.probplot(cxb, dist=stats.gamma, sparams=param, plot=ax)
        ax.set_title(P.data + ' minimal distance from SNPs inside gene to exons\ngamma with parameters ' + str(param))
        ax.get_lines()[0].set_markersize(3.0)
        print()
        print('Result as parameters of gamma distribution, to copy in parameters:')
        print("shape =", param[0])
        print("scale =", param[1])


if __name__ == "__main__":
    SNP_exon_analyse().fit_draw()
    plt.show()
