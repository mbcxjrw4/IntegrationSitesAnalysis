#!/usr/bin/R
library(GenomicFeatures)
library(pander)
suppressWarnings(gencode.v22.txdb <- GenomicFeatures::makeTxDbFromGFF('path_to_file/gencode.v22.annotation.gtf', format='gtf'))
ex <- GenomicFeatures::exonsBy(gencode.v22.txdb, 'gene') # extract exons
exSum <- sum(width(reduce(ex)))  # sum of non-overlapping exon sequence
geneSize <- data.frame(ensembl_id = names(exSum), gene_exon_size = as.integer(exSum), stringsAsFactors = F)
write.table(geneSize, 'intermediate/gencode.v22.gene.sizes.tsv', sep = '\t', col.names = T, row.names = F, quote = F, na = '')
pander::pander(sessionInfo())
