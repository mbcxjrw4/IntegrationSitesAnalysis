# make bed file containing gene coords and TPMs
fraietta_TPM <- read.table("resource/fraietta_TPM.tsv", sep = "\t", header = T, stringsAsFactors = F)

geneTPM <- merge(
		read.table("intermediate/gencode.v22.gene.ranges")[, c(1:3, 5)],
		fraietta_TPM[, c("ensembl_id", "mean_TPM")],
		by.x = "V5", by.y = "ensembl_id", all.x = F, all.y = T, sort = F
	)
names(geneTPM)[1:4] <- c("ensembl_id","chrom", "chromStart", "chromEnd")
geneTPM <- geneTPM[, c("chrom", "chromStart", "chromEnd", "ensembl_id", "mean_TPM")]
geneTPM$chrom <- factor(geneTPM$chrom, levels = unique(geneTPM$chrom))
geneTPM$chrom <- factor(geneTPM$chrom, levels = unique(geneTPM$chrom))
geneTPM <- geneTPM[order(geneTPM$chrom, geneTPM$chromStart), ] # bedtools assumes sorted
# write bed file
file_out <- file("intermediate/gencode.v22.gene.TPM.bed", "wb")
write.table(geneTPM, file_out, sep = "\t", col.names = F, row.names = F, quote = F, na = "")
close(file_out)
