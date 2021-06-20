# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# complement: Analyse position of genes and exons: overlapping

# Fields of sequence position (seq_pos)
# chromosome in first index of list
# 0: base position
# 1: 1 for exon start, 2 for exon end and 0 for SNP
# 2: name of gene
# 3: index of exon list if exon

import matplotlib.pyplot as plt
import SNP_exon_utils as U
import SNP_exon_param as P


class Exon_gene_analyse:

    def overlap_item(self, seq_pos):
        groups = [list() for i in range(24)]
        for c in range(24):
            group = set()
            seq_pos[c].sort()
            bp0, bp1 = 0, 0
            for item in seq_pos[c]:
                if item[1] == 1:
                    group.add(item[2])
                    bp0 = item[0]
                if item[1] == 2:
                    if item[0] != bp1:
                        bp1 = item[0]
                        groups[c].append((bp0, bp1, set(group)))
                    group.remove(item[2])
        return groups

    def list_items(self, groups, type):
        print('Overlapping', type)
        counts = [dict() for i in range(24)]
        for c in range(24):
            for g in groups[c]:
                length = len(g[2])
                if length in counts[c]:
                    counts[c][length] += 1
                else:
                    counts[c][length] = 1
                if len(g[2]) > 1:
                    print(c + 1, g[2])
        print('chr\toverlapping_size\t' + type + '_number')
        for c in range(24):
            for l in counts[c]:
                print(c + 1, l, counts[c][l], sep='\t')

    def overlap(self):
        groups = self.overlap_item(U.read_gene())
        self.list_items(groups, 'genes')
        exon_pos = U.read_exon()
        seq_pos = [list() for i in range(24)]
        for c in range(24):
            for ep in exon_pos[c]:
                seq_pos[c].append((ep[0], ep[1], ep[2] + '-' + str(ep[3])))
        groups = self.overlap_item(seq_pos)
        self.list_items(groups, 'exons')

    def exon_gap_between(self):
        exon_pos = [list() for i in range(24)]
        for ln in open(P.exon_file):
            sln = ln.split()
            exon_pos[int(sln[0]) - 1].append((float(sln[1]) + float(sln[2])) / 2)
        distances = list()
        for c in range(24):
            exon_pos[c].sort()
            it = iter(exon_pos[c])
            while True:
                try:
                    d1, d2 = next(it), next(it)
                    distances.append(d2 - d1)
                except StopIteration:
                    break
        distances.sort()
        egx, egy = U.cumul_number(distances)
        plt.title('Gap between exons en x log scale, normalized cumulative number en y')
        plt.xscale('log')
        plt.plot(egx, egy)
        threshold = 100000
        U.plot_histo(distances, 1000, title='gap between exons, max dist=' + str(threshold), max_x=threshold)
        plt.show()


if __name__ == "__main__":
    ega = Exon_gene_analyse()
    ega.overlap()
    ega.exon_gap_between()
