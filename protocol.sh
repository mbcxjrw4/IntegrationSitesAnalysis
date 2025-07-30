#!/bin/bash
# Pipeline to compare MP’s viral intergration sites to public data
set -euo pipefail

# define the subfolder where all scripts are located
SCRIPT_DIR="R"

# define the subfolder where all intermediate files are stored 
INTERMEDIATE_DIR="intermediate"

# Define the log file name
LOG_FILE="log.txt"

# Create the log file
touch $LOG_FILE

echo "Loading the data..."

# data read and preparation
Rscript "${SCRIPT_DIR}/dataReadInAndPreparation.R" 2> "log.txt"

echo "Running information content analysis..."

# Information content
# Get sequences
samtools faidx "${HOME}/path_to_genome/GRCh38.p13.genome.fa" -r "${INTERMEDIATE_DIR}/20bp.plus.region" > "${INTERMEDIATE_DIR}/20bp.plus.fa"
samtools faidx "${HOME}/path_to_genome/GRCh38.p13.genome.fa" -r "${INTERMEDIATE_DIR}/20bp.minus.region" -i > "${INTERMEDIATE_DIR}/20bp.minus.fa"
Rscript "${SCRIPT_DIR}/informationContent.R"

echo "Information content analysis finished. Running gene analysis..."

# gene analysis
Rscript "${SCRIPT_DIR}/geneAnalysis.R" 2> "log.txt"

echo "Gene analysis finished. Running open chromatin analysis..."

# Open chromatin
Rscript "${SCRIPT_DIR}/openChromatin.R" 2> "log.txt"

echo "Open chromatin analysis finished. Running GC content analysis..."

# GC content
# Calculate GC content in 1 kb bins
samtools faidx "${HOME}/path_to_genome/GRCh38.p13.genome.fa" -r "${INTERMEDIATE_DIR}/1kb.region" | awk 'BEGIN { IGNORECASE=1; a=0; c=0; g=0;t=0;}{ if (NR==1) next; a += gsub("A", ""); c += gsub("C", ""); g += gsub("G", ""); t += gsub("T", "");} />/ { print (c+g)/(a+c+g+t); a=0; c=0; g=0; t=0; next;} END { print (c+g)/(a+c+g+t);}' > "${INTERMEDIATE_DIR}/1kb.region.gc"

Rscript "${SCRIPT_DIR}/gcContent.R" 2> "log.txt"

echo "GC content analysis finished. Running clonality analysis..."

# Clonality
Rscript "${SCRIPT_DIR}/clonality.R" 2> "log.txt"

echo "Clonality finished. Generating the report..."

# Report
R -e "rmarkdown::render('${SCRIPT_DIR}/report.Rmd')"

mv R/report.html results/
mv log.txt results/
rm -rf intermediate/
