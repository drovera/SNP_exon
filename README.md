# A way to analyse GWAS data using exons

Genome Wide Association Studies allow to analyse the link between
frequency of single nucleotide polymorphism and phenotype by
comparing to a reference population.

The first step of this analysis is evaluating how the p-values of SNPs
resulting from this comparison are transferred to genes.

The method described in analyse_GWAS_by_exon.pdf is based 
on the base position of exons and the statistical relation
between positions of SNPs and positions of exons.

This method sheds light on the link between SNPs and
their effect on gene production, what is coherent
with some biological knowledge. The usual method
based on chi-deux statistics does not shed this light.

The manual to use Python scripts is SNP_exon_manual.pdf.
The reference data are in ref and GWAS data are to be extracted
from catalogs as NHGRI-EBI


daniel.rovera@gmail.com Institut Curie - Mines Paris Tech


