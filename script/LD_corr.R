# Compute LD correlation of one gene with LDlinkR

library(LDlinkR)
library(lubridate)
your_root <- "C:/Users/gisel/Documents/SNP_exon/"
source <- "ebi_006719"
your_token <- "beeaf956cd42"

# input gene and where to write
gene <-"TMEM167A_1"

# Output file where to write result
cFile <- paste0(your_root, source,"/LD_corr/", gene, "_c.txt")
cFile <- ""  # if you want to write in console

sFile <- paste0(your_root, source,"/LD_corr/", gene, "_g.txt")
snpwz <- read.table(sFile, sep=';')
snp <- snpwz[,1][!duplicated(snpwz[,1])] # all unique SNPs
cat(gene,"\n")
LD <- LDmatrix(snps = snp,
            pop = c("YRI", "CEU"), r2d = "d",
            genome_build = "grch38",
            token = your_token)
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
