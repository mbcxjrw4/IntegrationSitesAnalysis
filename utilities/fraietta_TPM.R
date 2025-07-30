# Script to calculate Gene expression (ex-vivo CD3+ cells)

# read the sizes
gencode.v22.gene.sizes <- read.table("intermediate/gencode.v22.gene.sizes.tsv", sep = "\t", header = T, stringsAsFactors = F)

# from [Fraietta *et al*., 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6117613/)
fname <- "resource/NIHMS981956-supplement-Supplementary_Table_5.xlsx"
fraietta_gene_counts <- openxlsx::read.xlsx(fname, sheet = 3, startRow = 4, colNames = F)
colnames(fraietta_gene_counts)[1] <- "gene_name"
fraietta_gene_counts <- fraietta_gene_counts[ fraietta_gene_counts$gene_name %in% gencode.v22.genes$gene_name, ]
fraietta_gene_counts <- merge(
  fraietta_gene_counts,
	data.frame(mcols(gencode.v22.genes)[, c("ensembl_id", "gene_name")]),
	by = "gene_name", all.x = T, all.y = F, sort = F
) # some gene_names have more than one id - some duplication occurs

# merge with counts
fraietta_gene_counts <- merge(
  fraietta_gene_counts,
  gencode.v22.gene.sizes,
  by = "ensembl_id", all.x = T, all.y = F, sort = F)

# deal with duplicated gene names (more than one ensembl ID) - use the biggest gene size
for (dup_gene_name in unique(fraietta_gene_counts$gene_name[duplicated(fraietta_gene_counts$gene_name)])) {
  dups <- fraietta_gene_counts[fraietta_gene_counts$gene_name == dup_gene_name, ]
  dups <- dups[order(dups$gene_exon_size, decreasing = T), ]
  fraietta_gene_counts <- fraietta_gene_counts[!fraietta_gene_counts$gene_name == dup_gene_name, ]
  fraietta_gene_counts <- rbind(fraietta_gene_counts, dups[1, ])
}

# calculate TPMs
RPK <- t(apply(fraietta_gene_counts, 1, function(x) {
  as.numeric(x[3:16]) * 1000 / as.numeric(x[17])
  }))  # reads per kb

TPM <- apply(RPK, 2, function(x) {x / (sum(x)/1000000)})

fraietta_gene_counts$mean_TPM <- apply( TPM, 1, function(x) round(mean(x, na.rm = T), 3) )
fraietta_TPM <- fraietta_gene_counts[ , c("ensembl_id", "gene_name", "gene_exon_size", "mean_TPM")]
# save the table
write.table(fraietta_TPM, "intermediate/fraietta_TPM.tsv", sep = "\t", col.names = T, row.names = F, quote = F, na = "")
