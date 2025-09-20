#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Load config params
params.input         = "input/intSites.tsv"
params.outdir        = "results"
params.intermediate  = "intermediate"
params.genome_fasta  = "${HOME}/path_to_genome/GRCh38.p13.genome.fa"
params.resources     = "resource"
params.report_rmd    = "report.Rmd"

// Channels
Channel.fromPath(params.input).set { raw_input }

// Step 1: Data preparation
process DataReadInAndPreparation {
    publishDir "${params.intermediate}", mode: 'copy'
    conda 'environment.yml'

    input:
    path raw_data

    output:
    path "20bp.plus.region"
    path "20bp.minus.region"
    path "IS.csv"

    script:
    """
    Rscript R/dataReadInAndPreparation.R \
        --patient_id config/variables.R \
        --input ${raw_data} \
        --outdir ${params.intermediate}
    """
}

// Step 2: Sequence logo
process SequenceLogo {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path plus_region from DataReadInAndPreparation.out.collect().map{ it[0] }
    path minus_region from DataReadInAndPreparation.out.collect().map{ it[1] }
    path is_csv from DataReadInAndPreparation.out.collect().map{ it[2] }

    output:
    path "logo20bp.png"
    path "IS.csv"

    script:
    """
    samtools faidx ${params.genome_fasta} -r ${plus_region} > ${params.intermediate}/20bp.plus.fa
    samtools faidx ${params.genome_fasta} -r ${minus_region} > ${params.intermediate}/20bp.minus.fa

    Rscript R/informationContent.R \
        --plus ${params.intermediate}/20bp.plus.fa \
        --minus ${params.intermediate}/20bp.minus.fa \
        --is ${is_csv} \
        --output ${params.outdir}/plots/logo20bp.png
    """
}

// Step 3: Gene analysis
process GeneAnalysis {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path is_csv from SequenceLogo.out.collect().map{ it[1] }

    output:
    path "IS_gene.png"
    path "ISGR.rds"

    script:
    """
    Rscript R/geneAnalysis.R \
        --is ${is_csv} \
        --genes ${params.resources}/gencode.v22.gene.ranges \
        --tpm ${params.resources}/fraietta_TPM.tsv \
        --plot ${params.outdir}/plots/IS_gene.png \
        --out ${params.intermediate}/ISGR.rds
    """
}

// Step 4: Open chromatin
process OpenChromatin {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path isgr_rds from GeneAnalysis.out.collect().map{ it[1] }

    output:
    path "IS_dnase.png"
    path "ISGR.rds"
    path "1kb.region"

    script:
    """
    Rscript R/openChromatin.R \
        --isgr ${isgr_rds} \
        --dnase ${params.resources}/wgEncodeRegDnaseClustered.txt.gz \
        --plot ${params.outdir}/plots/IS_dnase.png \
        --out ${params.intermediate}
    """
}

// Step 5: GC content
process GCContent {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path isgr_rds from OpenChromatin.out.collect().map{ it[1] }
    path region from OpenChromatin.out.collect().map{ it[2] }

    output:
    path "GC_content.png"
    path "integrationSiteData.tsv"
    path "ISGR.rds"

    script:
    """
    samtools faidx ${params.genome_fasta} -r ${region} | \
    awk 'BEGIN {IGNORECASE=1;a=0;c=0;g=0;t=0;}
         {if(NR==1)next;a+=gsub("A","");c+=gsub("C","");g+=gsub("G","");t+=gsub("T","");}
         />/ {print (c+g)/(a+c+g+t);a=0;c=0;g=0;t=0;next;}
         END {print (c+g)/(a+c+g+t);}' > ${params.intermediate}/1kb.region.gc

    Rscript R/gcContent.R \
        --isgr ${isgr_rds} \
        --gc ${params.intermediate}/1kb.region.gc \
        --out ${params.outdir}/integrationSiteData.tsv \
        --plot ${params.outdir}/plots/GC_content.png \
        --save ${params.intermediate}/ISGR.rds
    """
}

// Step 6: Clonality
process Clonality {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path isgr_rds from GCContent.out.collect().map{ it[2] }

    output:
    path "IS_clonality.png"
    path "clonalityData.tsv"

    script:
    """
    Rscript R/clonality.R \
        --isgr ${isgr_rds} \
        --bed ${params.resources}/GRCh38.p13_1Mb_geneTPM.bed \
        --plot ${params.outdir}/plots/IS_clonality.png \
        --out ${params.outdir}/clonalityData.tsv
    """
}

// Step 7: Final report
process Report {
    publishDir "${params.outdir}", mode: 'copy'
    conda 'environment.yml'

    input:
    path clonality_data from Clonality.out.collect().map{ it[1] }

    output:
    path "report.html"

    script:
    """
    Rscript -e "rmarkdown::render('${params.report_rmd}', output_file='report.html')"
    """
}

// Workflow
workflow {
    raw_input | DataReadInAndPreparation \
        | SequenceLogo \
        | GeneAnalysis \
        | OpenChromatin \
        | GCContent \
        | Clonality \
        | Report
}
