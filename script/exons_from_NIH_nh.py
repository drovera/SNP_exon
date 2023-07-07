# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# Create a file containing exon position and gene area extracted from NIH by chromosome
# use name in NIH
# See abnormalities in both cases

import SNP_exon_param as P

### launch in first exons_from_NIH_nh after extracting sequence files ###
# refSeq from NIH genome - https://www.ncbi.nlm.nih.gov/genome/?term=homo+sapiens+%5Borgn%5D
extract_from_NIH = P.ref_dir + 'NIHsequence/sequence_chr'

# Abnormality of extraction are listed in
NIH_abn = P.ref_dir + 'exon_abn.txt'

class Exons_from_NIH:

    def get_gene_name(self, line):
        gb = line.find('gene=')
        ge = line.find(']', gb)
        return line[gb + len('gene='):ge]

    def get_loc(self, chr, line):
        ''' localisation of exons '''
        locs = list()
        lb = line.find('location=')
        cp = line.find('complement', lb)
        jp = line.find('join', lb)
        if jp > -1:
            lb = line.find('(', jp)
            le = line.find(')', lb)
            location = line[lb + 1: le]
        else:
            if cp == -1:
                le = line.find(']', lb)
                location = line[lb + len('location='):le]
            else:
                lb = line.find('(', cp)
                le = line.find(')', lb)
                location = line[lb + 1: le]
        if '<' in location:
            location = location.replace('<', '')
        if '>' in location:
            location = location.replace('>', '')
        ls = location.split(',')
        for loc in ls:
            beg_end = loc.split('..')
            if len(beg_end) == 2:
                locs.append((beg_end[0], beg_end[1]))
            else:
                self.array_abn += 1
                self.abn_o.write(str(chr) + '\t' + self.get_gene_name(line) + '\t' + location + '\n')
        return locs

    def seq_to_exon(self):
        ''' write position of exons by gene and poisiton of genes '''
        repeated_gene = set()
        count_gene = list()
        count_exon = list()
        self.array_abn = 0
        self.abn_o.write('Abnormality in array of exons: \nchr\tgene\tpositions\n')
        for _ in range(24):
            count_gene.append(0)
            count_exon.append(0)
        for chr in range(1, 25):
            in_exon_file = extract_from_NIH + '{:02}'.format(chr) + '.txt'
            genes = set()
            inf = open(in_exon_file)
            line = inf.readline()
            while line:
                if not line.startswith('>'):
                    line = inf.readline()
                    continue
                gene = self.get_gene_name(line)
                if not gene:
                    line = inf.readline()
                    continue
                if gene in genes:
                    repeated_gene.add(gene)
                    line = inf.readline()
                    continue
                count_gene[chr - 1] += 1
                genes.add(gene)
                locs = self.get_loc(chr, line)
                for i in range(len(locs)):
                    count_exon[chr - 1] += len(locs)
                    self.pos_o.write(str(chr) + '\t' + str(locs[i][0]) + '\t' + str(locs[i][1])
                                     + '\t' + gene + '\t' + str(i) + '\n')
                self.gene_o.write(str(chr) + '\t' + str(locs[0][0]) + '\t' + str(locs[-1][1])
                                  + '\t' + gene + '\t' + str(len(locs)) + '\n')
                line = inf.readline()
        print(str(self.array_abn), 'Abnormality in array of exons\n')
        mess = str(len(repeated_gene)) + ' Repeated genes in sequences\n'
        print(mess)
        self.abn_o.write(mess)
        for gene in repeated_gene:
            self.abn_o.write(str(gene) + '\n')
        for i in range(24):
            print(i + 1, count_gene[i], 'genes', count_exon[i], 'exons', sep='\t')

    def write_seq(self):
        self.pos_o = open(P.exon_file, mode='w')
        self.gene_o = open(P.gene_file, mode='w')
        self.abn_o = open(NIH_abn, mode='w')
        self.seq_to_exon()
        self.pos_o.close()
        self.gene_o.close()
        self.abn_o.close()
        print('\nExon sequence save in', P.exon_file)
        print('Gene area save in', P.gene_file)
        print('Abnormalities save in', NIH_abn)


if __name__ == "__main__":
    Exons_from_NIH().write_seq()
