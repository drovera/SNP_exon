# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# common parameters: directories, files, data sources, parametres of scripts, results ...


############### directories #################
your_root = 'C:/Users/gisel/Documents/SNP_exon/'
source = 'ebi_006719'
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006719
ref_dir = your_root + 'ref/'
dta_dir = your_root + source + '/dta/'
ldcor_dir = your_root + source + '/LD_corr/'
result_dir = your_root + source + '/result/'
abnor_dir = your_root + source + '/abnor/'
#############################################

############### ref common files, format in manual ################
exon_file = ref_dir + 'exon_pos.txt'  # positions of numbered exons
gene_file = ref_dir + 'gene_pos.txt'  # positions of genes
###################################################################

################## parameters to read GWAS data ######################
GWAS_file = dta_dir + 'ebi_006719_download.txt'
title_line = True  # if title line generally
SNP_chr, SNP_name, SNP_bp, SNP_pv = 0, 1, 2, 6  # position of fields
inf_pv, sup_pv = 0, 1  # values limits of p_values
shape = 0.4538  # computed by SNP_exon_analyse.py
scale = 34094.3  # computed by SNP_exon_analyse.py
######################################################################

################ source common files, format in manual ###################
pvZ_file = dta_dir + source + '_pvZ.txt'  # valid SNP, p-value and Z-score
sed_file = dta_dir + source + '_sed.txt'  # distances between SNP and exons
gene_list = ldcor_dir + '_list_g.txt'
##########################################################################

################# result files ###########################
gene_Z_file_default = result_dir + source + '_Zgd.txt'
gene_Z_file_ldcorr = result_dir + source + '_Zgc.txt'
##########################################################

############## Linkage desiquilibrium correlation coefficient ##################
### values by default ###
def_LDcorr = 0.63  # computed by SNP_exon_corr.py
default_option = False  # use only default LD correlation coefficient above to compute Z-score of genes
################################################################################

### function reading Linkage desiquilibrium correlation coefficient ###
#  reading LD correlation files and return a dictionary[SNP1][SNP2] by group pf SNPs influencing a gene
#  may be modified if file format of LD correlation is different
def make_corr_dict(corr_file):
    snps_one = list()
    ld_cor = list()
    for cline in open(corr_file):
        cline = cline[0:-2]
        sln = cline.split(';')
        snps_one.append(sln[0])
        ld_cor.append(1.0)
        for i in range(1, len(sln)):
            try:
                ld_cor.append(float(sln[i]))
            except ValueError:
                ld_cor.append(def_LDcorr)
    corr_dict = dict()
    n = len(snps_one)
    for r in range(n):
        corr_dict[snps_one[r]] = dict()
    for r in range(n):
        for c in range(r, n):
            ind = r * n - r * (r + 1) // 2 + c
            corr_dict[snps_one[r]][snps_one[c]] = ld_cor[ind]
            corr_dict[snps_one[c]][snps_one[r]] = ld_cor[ind]
    return corr_dict
###########################################################################################

################### Labels of relative position cases ######################
cases = ['SNP inside exon',  # 0
         'SNP inside the same gene',  # 1
         'SNP inside the the 2 nearest genes',  # 2
         'SNP inside this gene but outside the other nearest gene',  # 3
         'SNP outside this gene but inside the other nearest gene',  # 4
         'SNP between exons and outside the 2 nearest genes ',  # 5
         'SNP at the extremities and outside all genes']  # 6
#############################################################################

################################# Chromosome size  ##########################################
chr_size = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
            138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
            83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415]
##############################################################################################
