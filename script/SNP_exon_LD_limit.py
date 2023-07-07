# Create a modified files due to limitations to search LD correlation
# The resulting file must be copied in LD directory

import os

root = "C:/Users/gisel/Documents/SNP_exon/LD/"
abn_data = "ebi_006719/abnor/"  # path of abnormal files
too_snp = ['CNTNAP2_1', 'CSMD1_1', 'EYS_1', 'LRP1B_1', 'RBFOX1_1', 'SGCZ_1']
threshold = 0.95  # threshold of part of sum Z-score * weight
create_file = True


def read_snp_list(sfile):
    ''' Read SNP and compute Z * weight '''
    snp_list = list()
    zws = 0.0
    for line in open(sfile):
        line = line[0:-1]
        sln = line.split(';')
        zw = abs(float(sln[1]) * float(sln[2]))
        zws = zws + zw
        sln.insert(0, zw)
        snp_list.append(sln)
    snp_list.sort(reverse=True)
    return zws, snp_list


def search_in_lim(too_spn):
    ''' Compute limits from threshold of sum Z * w'''
    in_limits = list()
    for gene in too_snp:
        sfile = root + abn_data + gene + "_g.txt"
        zws, snp_list = read_snp_list(sfile)
        cum, cums, = 0.0, list()
        snps = set()
        for s in snp_list:
            cum = cum + s[0] / zws
            snps.add(s[1])
            if cum > threshold:
                in_limits.append(snps)
                break
    print('gene_int_spl\tlimited_snp')
    for i in range(len(too_spn)):
        print(too_spn[i], len(in_limits[i]), sep='\t')
    return in_limits


def limit_snp_file(in_limits):
    ''' Create snp file from theshold, initial snp file is kept '''
    for i in range(len(too_snp)):
        sfile = root + abn_data + too_snp[i] + "_g.txt"
        ifile = root + abn_data + too_snp[i] + "_full_g.txt"
        if os.path.isfile(ifile):
            print('the initial file has been kept as', ifile)
            print('rename it to create the under threshold file')
            return
        os.rename(sfile, ifile)
        out = open(sfile, mode='w')
        for line in open(ifile):
            sln = line.split(';')
            if sln[0] in in_limits[i]:
                out.write(line)
        out.close()
        print('under threshold file in', sfile)
        print('initial snp file in', ifile)


if __name__ == "__main__":
    in_limits = search_in_lim(too_snp)
    if create_file:
        limit_snp_file(in_limits)