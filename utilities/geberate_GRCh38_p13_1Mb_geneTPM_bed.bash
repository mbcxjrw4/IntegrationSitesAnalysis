# Find overlapping genes and sum their TPMs
# make bed file containing gene coords and TPMs
Rscript generate_gencode_v22_gene_TPM_bed.R

# Sum gene TPMs within 1 Mbs
bedtools map -o sum -a resource/GRCh38.p13_1Mb.bed -b resource/gencode.v22.gene.TPM.bed > resource/GRCh38.p13_1Mb_geneTPM.bed
