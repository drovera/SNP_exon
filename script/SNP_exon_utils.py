# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# utils functions about file and graph

import numpy as np
from scipy import optimize, stats
import matplotlib.pyplot as plt
import SNP_exon_param as P


def read_gene():
    file = P.gene_file
    seq_pos = [list() for i in range(24)]
    for line in open(file):
        sq = line.split()
        seq_pos[int(sq[0]) - 1].append((int(sq[1]), 1, sq[3]))
        seq_pos[int(sq[0]) - 1].append((int(sq[2]), 2, sq[3]))
    return seq_pos

def read_exon():
    file = P.exon_file
    seq_pos = [list() for i in range(24)]
    for line in open(file):
        sq = line.split()
        seq_pos[int(sq[0]) - 1].append((int(sq[1]), 1, sq[3], int(sq[4])))
        seq_pos[int(sq[0]) - 1].append((int(sq[2]), 2, sq[3], int(sq[4])))
    return seq_pos


def read_SNP(log=False):
    file = P.SNP_file
    seq_pos = [list() for i in range(24)]
    snp_nb = 0
    SNP_pv = dict()
    in_ = open(file)
    if P.title_line:
        in_.readline()
        record = 2
    else:
        record = 1
    SNP_l = in_.readline()
    excluded, duplicated = np.zeros(24, dtype=int), np.zeros(24, dtype=int)
    ne, nd = 0, 0
    while SNP_l:
        try:
            SNP = SNP_l.split()
            c = int(SNP[P.SNP_chr]) - 1
            pv = float(SNP[P.SNP_pv])
            name = SNP[P.SNP_name]
            if name in SNP_pv:
                duplicated[c], nd = duplicated[c] + 1, nd + 1
                old_name, name = name, name + 'b'
                if log:
                    print('Duplicated SNP:', old_name, ', renamed:', name)
            if P.inf_pv < pv < P.sup_pv:
                SNP_pv[name] = pv
                seq_pos[c].append((int(SNP[P.SNP_bp]), 0, name))
                snp_nb += 1
            else:
                excluded[c], ne = excluded[c] + 1, ne + 1
                if log:
                    print('Excluded SNP:', name, 'pv=', pv)
        except Exception as excpt:
            print('Abnormality at record', record, excpt.args)
        SNP_l = in_.readline()
        record += 1
    print(snp_nb, 'p_values in [', P.inf_pv, ', ', P.sup_pv, ']')
    print('red from', record, 'lines from', file)
    if nd > 0 or ne > 0:
        print('Warning about reading SNP:', nd, 'duplicated,', ne, 'out limits')
        print('To get details set log = True ')
    if log:
        print('chr\tSNP\tNumber of excluded SNP, p-value out limits', P.inf_pv, P.sup_pv)
        for c in range(24):
            print(c + 1, len(seq_pos[c]), excluded[c], duplicated[c], sep='\t')
        print('sum', ne, nd, sep='\t')
    return seq_pos, SNP_pv

def plot_histo(x, bin_nb, title=None, min_x=None, max_x=None, funct=None, param=None, noyticks=False):
    if not min_x:
        min_x = min(x)
    if not max_x:
        max_x = max(x)
    bins = np.linspace(min_x, max_x, num=bin_nb)
    hist, _ = np.histogram(x, bins=bins, density=True)
    width = np.diff(bins)
    center = (bins[:-1] + bins[1:]) / 2
    _, ax = plt.subplots()
    if title:
        ax.set_title('Histogram of ' + title)
    if noyticks:
        ax.set_yticks([])
    ax.bar(center, hist, align='center', width=width)
    if funct:
        ax.plot(center, funct(center, *param), color='r')

def curve_fit(xdata, ydata, func, bounds=None):
    if bounds:
        popt, pcov = optimize.curve_fit(func, xdata, ydata, bounds=bounds)
    else:
        popt, pcov = optimize.curve_fit(func, xdata, ydata)
    print('Parameters:', popt, 'with bounds', bounds)
    perr = np.sqrt(np.diag(pcov))
    print('Estimated covariance, standard deviation errors on the parameter:', perr)
    return popt

def std_adjust_2param(xdata, ydata, func, start, delta, shrink, max_iter):
    delta = np.array(delta)
    min_std = 1
    go_on, iter = True, 0
    while go_on and iter < max_iter:
        go_on, iter = False, iter + 1
        for i1 in range(-1, 2):
            for i2 in range(-1, 2):
                p0, p1 = start[0] + delta[0] * i1, start[1] + delta[1] * i2
                std = np.sum(np.square(func(xdata, p0, p1) - ydata)) / len(xdata)
                if std < min_std:
                    min_std, start = std, (p0, p1)
                    go_on = True
        delta = shrink * delta
    if iter == max_iter:
        print('Result by max iteration:', end='\t')
    else:
        print('Result by shrinkage:', end='\t')
    print('param:', start, 'std:', min_std, sep='\t')
    return start

def cumul_number(xd):
    cx = np.zeros(len(xd) + 1)
    cy = np.zeros(len(xd) + 1)
    for i in range(1, len(xd) + 1):
        cx[i] = xd[i - 1]
        cy[i] = i / len(xd)
    return cx, cy

def gamma_draw_fit(cx, cy, funct, param, title=None):
    plt.figure(figsize=(7.0, 7.0))
    if title:
        plt.title(title)
    ey = funct(cx, *param)
    gy = ey - cy
    std = np.sum(np.square(gy)) / cy.size
    print('Standard deviation errors on the data:', std)
    plt.xscale('log')
    plt.plot(cx, cy, label='data')
    plt.plot(cx, ey, label='fit')
    plt.plot(cx, gy, label='gap')
    plt.legend()
    return std

def gamma_cdf(x, shape, scale):
    return stats.gamma.cdf(x, shape, scale=scale)

def gamma_pdf(x, shape, scale):
    return stats.gamma.pdf(x, shape, scale=scale)

# Weight function, above shape and scale parameter
def weight_f(d):
    return gamma_pdf(d, P.shape, P.scale)

