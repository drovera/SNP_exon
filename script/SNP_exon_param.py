# daniel.rovera@gmail.com Institut Curie - Mines Paris Tech
# common parameters: directories, files, data sources, parametres of scripts, results ...


############### directories #################
your_root = 'C:/Users/gisel/Documents/SNP_exon/'
source = 'ebi_006719'
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006719
ref_dir = your_root + 'ref/'
dta_dir = your_root + source + '/dta/'
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
shape = 0.4538  # computed by SNP_exon_distance_fit.py
scale = 34094.3  # computed by SNP_exon_distance_fit.py
######################################################################

################ source common files, format in manual ###################
pvZ_file = dta_dir + source + '_pvZ.txt'  # valid SNP, p-value, Z-score and bp
sed_file = dta_dir + source + '_sed.txt'  # distances between SNP and exons
##########################################################################

################# result files ###########################
gene_Z_file_fitLD = result_dir + source + '_Zgd.txt'
gene_Z_file_LDcorr = result_dir + source + '_Zgc.txt'
##########################################################

############## Linkage desiquilibrium ##################
ldcor_dir = your_root + source + '/LD_corr/'
gene_list = ldcor_dir + '_list_g.txt'
### Linkage desiquilibrium correlation coefficient by default ###
def_LDcorr = 0.63  # computed by LD_corr_analyse_default_value.py
################################################################################

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



