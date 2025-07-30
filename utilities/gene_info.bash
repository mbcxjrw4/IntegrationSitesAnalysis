# Comprehensive gene annotation (Gencode v22, hg38)
# get gene info from GTF file (gencode v22)
awk '$3 == "gene" { print $1,$4,$5,$7,$10,$16,$12 }' path_to_file/gencode.v22.annotation.gtf | tr -d ";\""  > ../intermediate/gencode.v22.gene.ranges
