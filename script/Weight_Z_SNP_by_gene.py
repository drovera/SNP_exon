# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Create files by genes containing SNP, Z-score and weight
# produce a list usable to search the LD correlation coefficients

import SNP_exon_param as P
import SNP_exon_utils as U

def save_by_gene():
    SNP_Z = dict()
    file = P.pvZ_file
    for ln in open(file):
        sln = ln.split()
        SNP_Z[sln[0]] = sln[2]
    print(file, 'red')
    file_in = P.sed_file
    prev_SNP = ''
    snp_zws = list()
    for ln in open(file_in):
        sln = ln.split()
        if sln[1] == '0':
            if sln[2] != prev_SNP:  # keep only one SNP when internal
                snp_zws.append((sln[3] + '_0', sln[2], SNP_Z[sln[2]], '1.0'))
            prev_SNP = sln[2]
        if sln[1] in {'1', '2', '3'}:
            wgt = str(U.weight_f(float(sln[5])))
            snp_zws.append((sln[3] + '_1', sln[2], SNP_Z[sln[2]], str(wgt)))
    snp_zws.sort()  # sort by gene, twice SNPs side by side as in sef_file
    print(file_in, 'red')
    it = iter(snp_zws)
    snp_zw = next(it)
    while True:
        gene = snp_zw[0]
        g_file = P.ldcor_dir + gene + '_g.txt'
        gn_out = open(g_file, mode='w')  # one file by gene containing SNP, Z, W
        try:
            while gene == snp_zw[0]:
                gn_out.write(snp_zw[1] + ';' + snp_zw[2] + ';' + snp_zw[3] + '\n')
                snp_zw = next(it)
            gn_out.close()
            print(g_file, 'written')
        except StopIteration:
            gn_out.close()
            print(g_file, 'written')
            break
    gl_out = open(P.gene_list, mode='w')  # list of gene_position
    gl_out.write("gene_int_spl;snp_nb;flag\n")
    it = iter(snp_zws)
    snpzw = next(it)
    gene, ante_SNP, snp_nb = '', '', 0
    while True:
        try:
            if gene != snpzw[0]:
                if gene != '':
                    if snp_nb == 1 or (snp_nb == 2 and ante_SNP == prev_SNP):
                        gl_out.write(gene + ';' + str(snp_nb) + ';' + '1' + '\n')
                    else:
                        gl_out.write(gene + ';' + str(snp_nb) + ';' + '0' + '\n') # only one SNP
                gene = snpzw[0]
                snp_nb = 0
            snp_nb += 1
            ante_SNP = prev_SNP
            prev_SNP = snpzw[1]
            snpzw = next(it)
        except StopIteration:
            break
    if snp_nb == 1 or (snp_nb == 2 and ante_SNP == prev_SNP):
        gl_out.write(gene + ';' + str(snp_nb) + ';' + '1' + '\n')
    else:
        gl_out.write(gene + ';' + str(snp_nb) + ';' + '0' + '\n')
    gl_out.close()
    print('gene files and list written')

if __name__ == "__main__":
    save_by_gene()