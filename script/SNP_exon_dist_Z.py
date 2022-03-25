# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Draw graph Z-score of SNP in function of min distance SNP exons
# only in the case of splicing effect
# Graw graph to help choosing a threshold of Z-score to limit analysis of SNP effect

# Z_dist: list of (Z, min distance to exon, SNP)

import SNP_exon_param as P
import numpy as np
import matplotlib.pyplot as plt

class SNP_exon_dist_Z:

    def read_distance(self):
        min_dist_btw = dict()
        file = P.sed_file
        in_ = open(file)
        while True:
            ln = in_.readline()
            if ln == '':
                break
            sln1 = ln.split()
            ln = in_.readline()
            sln2 = ln.split()
            if int(sln1[1]) in {1, 2}:
                min_dist_btw[sln1[2]] = min(int(sln1[5]), int(sln2[5]))
            if int(sln1[1]) == 3:
                min_dist_btw[sln1[2]] = int(sln1[5])
            if int(sln2[1]) == 3:
                min_dist_btw[sln2[2]] = int(sln2[5])
        print('distance SNP exons red from', file)
        print('Minimal distance only concerns SNP outside exons')
        print('SNP inside genes:', len(min_dist_btw))
        return min_dist_btw

    def read_SNP_Z_dist(self):
        min_dist_btw = self.read_distance()
        file = P.pvZ_file
        Z_dist = list()
        for ln in open(file):
            sln = ln.split()
            if sln[0] in min_dist_btw:
                Z_dist.append((float(sln[2]), min_dist_btw[sln[0]], sln[0]))
        print('Z-score of SNP red from', file)
        return Z_dist

    def d_z_array(self, Z_dist):
        size = len(Z_dist)
        d, z = np.zeros(size, dtype=int), np.zeros(size, dtype=float)
        i = 0
        for z_d in Z_dist:
            z[i], d[i] = z_d[0], z_d[1]
            i += 1
        return d, z

    def draw_d_z(self, Z_dist):
        d, z = self.d_z_array(Z_dist)
        plt.scatter(d, z, s=1)
        plt.title('Z-score, min distance SNP exons for ' + P.data)
        plt.xlabel('min distance')
        plt.ylabel('Z-score')
        plt.axhline(c='black')

    def highest_d_z(self, Z_dist, part):
        md, mz, nb, mx = list(), list(), list(), 0
        n, mn = 0, int(len(Z_dist) * part)
        for z_d in Z_dist:
            mz.append(z_d[0])
            if z_d[1] > mx:
                mx = z_d[1]
            md.append(mx)
            n += 1
            nb.append(n)
            if n > mn:
                break
        fig, ax = plt.subplots(1, 2)
        fig.suptitle('Part ' + str(part) + ' of all SNP for ' + P.data)
        ax[0].set_xlabel('max of min distance')
        ax[0].set_ylabel('Z-score')
        ax[0].plot(md, mz)
        ax[1].set_xlabel('number')
        ax[1].plot(nb, mz)

    def print_top(self, top, Z_dist):
        it = iter(Z_dist)
        print('SNP\tZ\tmin_dist')
        for _ in range(top):
            z_d = next(it)
            print(z_d[2], z_d[0], z_d[1], sep='\t')

if __name__ == "__main__":
    sed = SNP_exon_dist_Z()
    Z_dist = sed.read_SNP_Z_dist()
    sed.draw_d_z(Z_dist)
    Z_dist.sort(reverse=True)
    sed.print_top(100, Z_dist)
    sed.highest_d_z(Z_dist, 0.01)
    sed.highest_d_z(Z_dist, 0.001)
    sed.highest_d_z(Z_dist, 0.0001)
    plt.show()