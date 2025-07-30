# generate resource/gencode.v22.gene.ranges
# Comprehensive gene annotation (Gencode v22, hg38)
# get gene info from GTF file (gencode v22)
awk '$3 == "gene" { print $1,$4,$5,$7,$10,$16,$12 }' path_to_file/gencode.v22.annotation.gtf | tr -d ";\""  > ../resource/gencode.v22.gene.ranges

# generate resource/fraietta_TPM.tsv
# extract gene sizes from gencode.v22 gtf
Rscript "gencode.v22.gene.sizes.R"

# calculate Gene expression (ex-vivo CD3+ cells)
Rscript "fraietta_TPM.R"

# generate resource/GRCh38.p13_1Mb_geneTPM.bed
# Find overlapping genes and sum their TPMs
# make bed file containing gene coords and TPMs
Rscript "generate_gencode_v22_gene_TPM_bed.R"

# make 1 Mb regions
bedtools makewindows -g path_to_file/GRCh38.p13/STAR/chrNameLength.txt -w 1000000 | grep '^chr' > ../intermediate/GRCh38.p13_1Mb.bed

# Sum gene TPMs within 1 Mbs
bedtools map -o sum -a resource/GRCh38.p13_1Mb.bed -b resource/gencode.v22.gene.TPM.bed > resource/GRCh38.p13_1Mb_geneTPM.bed
