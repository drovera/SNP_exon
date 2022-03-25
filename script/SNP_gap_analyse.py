# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Complement: fitting distance between SNP, not used to compute p-value of genes

import matplotlib.pyplot as plt
import SNP_exon_utils as U
import SNP_exon_param as P

class SNP_distance:
    def read_distance(self):
        distances = list()
        in_ = open(P.SNP_file)
        in_.readline()
        ln = in_.readline()
        p_chr = 1
        p_d = 0
        while ln:
            sln = ln.split()
            chr = sln[P.SNP_chr]
            if chr == p_chr:
                if p_d == 0:
                    p_d = float(sln[P.SNP_bp])
                else:
                    d = float(sln[P.SNP_bp])
                    dist = d - p_d
                    distances.append(dist)
                    p_d = d
            else:
                p_d = 0
                p_chr = chr
            ln = in_.readline()
        distances.sort()
        return distances

    def fit_distance(self):
        print('Distance Between SNP, all SNP')
        dist = self.read_distance()
        cx, cy = U.cumul_number(dist)
        param = U.curve_fit(cx, cy, U.gamma_cdf)
        U.gamma_draw_fit(cx, cy, U.gamma_cdf, param, title='Fit distance between SNP, ' + P.datas[P.data_src])

if __name__ == "__main__":
    SNP_distance().fit_distance()
    plt.show()