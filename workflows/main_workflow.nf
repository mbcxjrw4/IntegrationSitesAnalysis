#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ----------------------
// Parameters
// ----------------------
params.config_file   = "config/config.yaml"
params.input         = "input/intSites.tsv"
params.outdir        = "results"
params.genome_fasta  = "${HOME}/path_to_genome/GRCh38.p13.genome.fa"
params.report_rmd    = "report.Rmd"

// ----------------------
// Channels
// ----------------------
Channel.fromPath(params.input).set { raw_input }

// ----------------------
// Processes
// ----------------------

// Step 1: Data preparation
process DataReadInAndPreparation {
    conda 'environment.yml'

    input:
    path raw_data

    output:
    path "20bp.plus.region",  emit: plus_region
    path "20bp.minus.region", emit: minus_region
    path "IS.csv",            emit: is_csv

    script:
    """
    Rscript R/dataReadInAndPreparation.R \
        --input ${raw_data} \
        --config ${params.config_file}
    """
}

// Step 2: Sequence logo
process SequenceLogo {
    conda 'environment.yml'

    input:
    path plus_region from DataReadInAndPreparation.out.plus_region
    path minus_region from DataReadInAndPreparation.out.minus_region
    path is_csv from DataReadInAndPreparation.out.is_csv

    output:
    path "logo20bp.png",   emit: logo
    path "IS.updated.csv", emit: is_csv_updated

    script:
    """
    samtools faidx ${params.genome_fasta} -r ${plus_region} > 20bp.plus.fa
    samtools faidx ${params.genome_fasta} -r ${minus_region} > 20bp.minus.fa

    Rscript R/informationContent.R \
        --plus 20bp.plus.fa \
        --minus 20bp.minus.fa \
        --is ${is_csv} \
        --output logo20bp.png \
        --update IS.updated.csv
    """
    publishDir "${params.outdir}/plots", mode: 'copy', pattern: "*.png"
}

// Step 3: Gene analysis
process GeneAnalysis {
    conda 'environment.yml'

    input:
    path is_csv from SequenceLogo.out.is_csv_updated

    output:
    path "IS_gene.png",   emit: gene_plot
    path "ISGR.rds",      emit: isgr_rds
    path "geneData.tsv",  emit: genedata_tsv

    script:
    """
    Rscript R/geneAnalysis.R \
        --is ${is_csv} \
        --config ${params.config_file} \
        --plot IS_gene.png \
        --out  geneData.tsv \
        --isgr ISGR.rds
    """
    publishDir "${params.outdir}/plots", mode: 'copy', pattern: "*.png" 
    publishDir "${params.outdir}",       mode: 'copy', pattern: "*.tsv"
}

// Step 4: Open chromatin
process OpenChromatin {
    conda 'environment.yml'

    input:
    path isgr_rds from GeneAnalysis.out.isgr_rds

    output:
    path "IS_dnase.png",  emit: is_dnase_plot
    path "ISGR.rds",      emit: isgr_rds_updated
    path "1kb.region",    emit: region

    script:
    """
    Rscript R/openChromatin.R \
        --isgr ${isgr_rds} \
        --config ${params.config_file} \
        --plot IS_dnase.png \
        --update ISGR.rds
    """
    publishDir "${params.outdir}/plots", mode: 'copy'
}

// Step 5: GC content
process GCContent {
    conda 'environment.yml'

    input:
    path isgr_rds from OpenChromatin.out.isgr_rds_updated
    path region   from OpenChromatin.out.region

    output:
    path "GC_content.png",          emit: gc_plot
    path "integrationSiteData.tsv", emit: is_data_tsv
    path "ISGR.rds",                emit: isgr_rds_updated2

    script:
    """
    samtools faidx ${params.genome_fasta} -r ${region} | \
    awk 'BEGIN {IGNORECASE=1;a=0;c=0;g=0;t=0;}
         {if(NR==1)next;a+=gsub("A","");c+=gsub("C","");g+=gsub("G","");t+=gsub("T","");}
         />/ {print (c+g)/(a+c+g+t);a=0;c=0;g=0;t=0;next;}
         END {print (c+g)/(a+c+g+t);}' > 1kb.region.gc

    Rscript R/gcContent.R \
        --isgr ${isgr_rds} \
        --gc 1kb.region.gc \
        --out integrationSiteData.tsv \
        --plot GC_content.png \
        --update ISGR.rds
    """
    publishDir "${params.outdir}/plots", mode: 'copy'
    publishDir "${params.outdir}",       mode: 'copy', pattern: "*.tsv"
}

// Step 6: Clonality
process Clonality {
    conda 'environment.yml'

    input:
    path isgr_rds from GCContent.out.isgr_rds_updated2

    output:
    path "IS_clonality.png", emit: is_clonality_plot
    path "clonalityData.tsv", emit: clonality_data_tsv

    script:
    """
    Rscript R/clonality.R \
        --isgr ${isgr_rds} \
        --config ${params.config_file} \
        --plot IS_clonality.png \
        --out clonalityData.tsv
    """
    publishDir "${params.outdir}/plots", mode: 'copy'
    publishDir "${params.outdir}",       mode: 'copy', pattern: "*.tsv"
}

// Step 7: Final report
process Report {
    conda 'environment.yml'

    input:
    // A single input declaration for a tuple containing all the files
    tuple path(logo_plot), path(gene_plot), path(dnase_plot), path(gc_content_plot), path(clonality_plot)

    output:
    path "report.html"

    script:
    """
    Rscript -e "
    rmarkdown::render(
        '${params.report_rmd}', 
        output_file='report.html', 
        params = list(
            logo20bp_png = ${logo_plot.baseName},
            IS_gene_png = ${gene_plot.baseName},
            IS_dnase_png = ${dnase_plot.baseName},
            GC_content_png = ${gc_content_plot.baseName},
            IS_clonality_png = ${clonality_plot.baseName}
        )
    )"
    """
    publishDir "${params.outdir}", mode: 'copy'
}

// ----------------------
// Workflow definition
// ----------------------
workflow {
    DataReadInAndPreparation(raw_input) \
        | SequenceLogo \
        | GeneAnalysis \
        | OpenChromatin \
        | GCContent \
        | Clonality 

    // Combine all outputs into a single channel containing a tuple
    report_inputs = SequenceLogo.out.logo
        .mix(GeneAnalysis.out.gene_plot)
        .mix(OpenChromatin.out.is_dnase_plot)
        .mix(GCContent.out.gc_plot)
        .mix(Clonality.out.is_clonality_plot)
        .collect()

    // Pass the single tuple to the Report process
    Report(report_inputs)
}

