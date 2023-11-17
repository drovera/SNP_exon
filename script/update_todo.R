# update temporary file _list_g$.txt
# which contains the list of genes to be processed
# thing of setting flag to -1 if not corrected error on gene_position

your_root <- "C:/Users/gisel/Documents/SNP_exon/"
source <- "ebi_006719"

path <- paste0(your_root, source)
gFile <- paste0(your_root, source, "/_list_g.txt")

genes <- read.table(gFile, header = TRUE, sep=';')

lf<-list.files(path=path,  pattern = "*_c.txt", full.names=FALSE)
cFiles <- substring(lf, 1, nchar(lf) - 6)
genes_to <- genes[genes$flag == 0,]
genes_todo <- subset(genes_to, !gene_int_spl %in% cFiles)

temp_gFile <- paste0(your_root, source, "/_list_g$.txt")
write.table(genes_todo,temp_gFile , append = FALSE, sep = ";",row.names = FALSE, col.names = TRUE)

cat("processed", length(cFiles), "\n")
cat("remains to be done", nrow(genes_todo), "\n")

