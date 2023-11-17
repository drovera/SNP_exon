# Compute LD correlation with LDlinkR with log to manage error
# in _list_g.txt 0 to extract, 1 extracted or one SNP
# and negative number -1 by hand if error
# follow the progress with log you can complete with error of LDlinkR

library(LDlinkR)
library(lubridate)

your_root <- "C:/Users/gisel/Documents/SNP_exon/"
source <- "ebi_006719"
your_token <- "beeaf956cd42"

log_file <- "C:/Users/gisel/Documents/R/LD_corr/LDlinkR.log"

cat(format(now()),"\n", file=log_file, append = TRUE)

gFile <- paste0(root, source, "/_list_g$.txt")
genes <- read.table(gFile, header = TRUE, sep=';')

gene_nb <- nrow(genes)
for(r in 1:gene_nb){
  if(genes[r, 3] == '0'){
    sFile <- paste0(your_root, source,"/LD_corr/", genes[r, 1], "_g.txt")
    snpwz <- read.table(sFile, sep=';')
    snp <- snpwz[,1][!duplicated(snpwz[,1])]# all unique SNPs
    cat(r, genes[r, 1], genes[r, 2], "\n", file = log_file, append = TRUE)
    cat(genes[r, 1],"\n")
    LD <- LDmatrix(snps = snp,
                pop = c("YRI", "CEU"), r2d = "d",
                genome_build = "grch38",
                token = your_token)
    cFile <- paste0(your_root, source,"/LD_corr/", genes[r, 1], "_c.txt")
    cat(NULL, file = cFile)
    n <- length(LD) - 1
    for (r in 2:n) {
      line <- LD[r - 1, 1]
      for (c in r:n) {
        rc <- LD[c, r]
        line <- paste(line, rc, sep = ";")
      }
      cat(line, "\n", file = cFile, append = TRUE)
    }
    cat(LD[n, 1], "\n", file = cFile, append = TRUE)
  }
}