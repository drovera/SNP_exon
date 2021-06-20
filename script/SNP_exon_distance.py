# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Distance SNP exon saved in a _sed files, described in SNP_exon_param :

# Fields of sequence position (seq_pos)
# chromosome in first index of list
# 0: base position
# 1: 1 for exon start, 2 for exon end and 0 for SNP
# 2: identifier for SNP, index of exon list

# Fields of direct reverse sequence position (dr_seq_pos) by chromosome
# 0: base position of SNP
# 1: name of SNP
# 2: base position of last extremity of exon when running chromosome
# 3: name of gene
# 4: exon order of gene

import numpy as np
import math
from operator import itemgetter
import SNP_exon_utils as U
import SNP_exon_param as P


class SNP_exon_distance:

    def create_seq(self):
        self.seq_pos = U.read_exon()
        self.SNP_pos, _ = U.read_SNP()
        for c in range(24):
            self.seq_pos[c].extend(self.SNP_pos[c])
        [self.seq_pos[i].sort() for i in range(24)]
        self.SNP_io = self.SNP_in_gene()

    def SNP_in_gene(self):
        seq_pos = U.read_gene()
        SNP_io = [dict() for i in range(24)]
        for c in range(24):
            seq_pos[c].extend(self.SNP_pos[c])
        [seq_pos[i].sort() for i in range(24)]
        for c in range(24):
            in_gene = set()
            for sq in seq_pos[c]:
                if sq[1] == 1:
                    in_gene.add(sq[2])
                elif sq[1] == 2:
                    in_gene.remove(sq[2])
                else:
                    SNP_io[c][sq[2]] = set(in_gene)
        return SNP_io

    def direct_reverse(self, chr_seq_pos):
        exon = None
        dir_seq = list()
        rev_seq = list()
        for item in chr_seq_pos:
            if item[1] > 0:
                exon = item
            else:
                if exon:
                    dir_seq.append((item[0], item[2], exon[0], exon[2], exon[3]))
        exon = None
        for item in reversed(chr_seq_pos):
            if item[1] > 0:
                exon = item
            else:
                if exon:
                    rev_seq.append((item[0], item[2], exon[0], exon[2], exon[3]))
        dr_seq_chr = dir_seq
        dr_seq_chr.extend(rev_seq)
        dr_seq_chr.sort(key=itemgetter(1, 2))
        return dr_seq_chr

    def to_str(self, item, case):
        if case == 0:
            dist = 0
        else:
            dist = int(math.fabs(item[2] - item[0]))
            if dist == 0:
                case = 0
        s = str(self.c + 1) + '\t' + str(case) + '\t' + item[1] + '\t' + item[3] + '\t' + str(item[4])
        s = s + '\t' + str(dist) + '\n'
        return dist, s

    def sed_write(self, item1, case1, item2, case2):
        d1, s1 = self.to_str(item1, case1)
        d2, s2 = self.to_str(item2, case2)
        if d1 > 0 and d2 > 0:
            self.o_sed.write(s1)
            self.o_sed.write(s2)
        elif d1 == 0:
            case1, case2 = 0, 0
            self.o_sed.write(s1)
            self.o_sed.write(s1)
        else:
            case1, case2 = 0, 0
            self.o_sed.write(s2)
            self.o_sed.write(s2)
        self.count[self.c][case1] += 1
        self.count[self.c][case2] += 1

    def proceed_chr(self, chr_dr_seq, chr_SNP_io):
        it = iter(chr_dr_seq)
        prev = next(it)
        while True:
            try:
                item = next(it)
                if prev[1] == item[1]:
                    if prev[3] == item[3] and prev[4] == item[4]:
                        self.sed_write(prev, 0, prev, 0)
                    else:
                        if prev[3] == item[3]:
                            self.sed_write(prev, 1, item, 1)
                        else:
                            if prev[3] in chr_SNP_io[prev[1]] and item[3] in chr_SNP_io[prev[1]]:
                                self.sed_write(prev, 2, item, 2)
                            elif prev[3] in chr_SNP_io[prev[1]]:
                                self.sed_write(prev, 3, item, 4)
                            elif item[3] in chr_SNP_io[prev[1]]:
                                self.sed_write(prev, 4, item, 3)
                            else:
                                self.sed_write(prev, 5, item, 5)
                    prev = next(it)
                else:
                    self.sed_write(prev, 6, prev, 6)
                    prev = item
            except StopIteration:
                break

    def proceed(self):
        self.o_sed = open(P.sed_file, mode='w')
        self.count = np.zeros((24, len(P.cases)), dtype=int)
        self.create_seq()
        for self.c in range(24):
            if self.SNP_pos[self.c]:
                self.proceed_chr(self.direct_reverse(self.seq_pos[self.c]), self.SNP_io[self.c])
        self.o_sed.close()
        print('distance between SNP and exons written in', P.sed_file)
        self.print_count()

    def print_count(self):
        print('chr\tSNP\toutside gene\tin one gene\tin more one')
        for c in range(24):
            n0, n1, n2 = 0, 0, 0
            for snp in self.SNP_io[c]:
                if len(self.SNP_io[c][snp]) == 0:
                    n0 += 1
                elif len(self.SNP_io[c][snp]) == 1:
                    n1 += 1
                else:
                    n2 += 1
            print(str(c + 1), len(self.SNP_io[c]), n0, n1, n2, sep='\t')
        for case in range(len(P.cases)):
            print(case, ':', P.cases[case])
        print('chr\tcase0\tcase1\tcase2\tcase3\tcase4\tcase5\tcase6')
        for c in range(24):
            print(c + 1, end='')
            for case in range(len(P.cases)):
                print('\t', self.count[c][case], sep='', end='')
            print()


if __name__ == "__main__":
    SNP_exon_distance().proceed()
