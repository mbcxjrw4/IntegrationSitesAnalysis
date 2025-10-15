#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ----------------------
// Parameters
// ----------------------
params.config_file   = "config/config.yaml"
params.input         = "input/intSites.tsv"
params.outdir        = "results"
params.genome_fasta  = null  // Must be provided via command line or config
params.report_rmd    = "report.Rmd"

// Resource defaults
params.max_cpus      = 4
params.max_memory    = '8.GB'
params.max_time      = '2.h'

// ----------------------
// Input validation
// ----------------------
if (!params.genome_fasta) {
    error "ERROR: --genome_fasta must be specified"
}

// ----------------------
// Processes
// ----------------------

// Step 1: Data preparation
process DataReadInAndPreparation {
    tag "Data preparation"
    label 'process_low'
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
    tag "Sequence logo generation"
    label 'process_medium'
    conda 'environment.yml'

    input:
    path plus_region
    path minus_region
    path is_csv

    output:
    path "logo20bp.png",   emit: logo
    path "IS.updated.csv", emit: is_csv_updated

    script:
    """
    # Index genome if not already indexed
    if [ ! -f "${params.genome_fasta}.fai" ]; then
        samtools faidx ${params.genome_fasta}
    fi

    # Extract sequences
    samtools faidx ${params.genome_fasta} -r ${plus_region} > 20bp.plus.fa
    samtools faidx ${params.genome_fasta} -r ${minus_region} > 20bp.minus.fa

    # Generate sequence logo
    Rscript R/informationContent.R \
        --plus 20bp.plus.fa \
        --minus 20bp.minus.fa \
        --is ${is_csv} \
        --output logo20bp.png \
        --update IS.updated.csv
    """
}

// Step 3: Gene analysis
process GeneAnalysis {
    tag "Gene analysis"
    label 'process_medium'
    conda 'environment.yml'

    input:
    path is_csv

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
}

// Step 4: Open chromatin
process OpenChromatin {
    tag "Open chromatin analysis"
    label 'process_medium'
    conda 'environment.yml'

    input:
    path isgr_rds

    output:
    path "IS_dnase.png",       emit: is_dnase_plot
    path "ISGR_chromatin.rds", emit: isgr_rds_updated
    path "1kb.region",         emit: region

    script:
    """
    Rscript R/openChromatin.R \
        --isgr ${isgr_rds} \
        --config ${params.config_file} \
        --plot IS_dnase.png \
        --update ISGR_chromatin.rds
    """
}

// Step 5: GC content
process GCContent {
    tag "GC content analysis"
    label 'process_medium'
    conda 'environment.yml'

    input:
    path isgr_rds
    path region

    output:
    path "GC_content.png",          emit: gc_plot
    path "integrationSiteData.tsv", emit: is_data_tsv
    path "ISGR_gc.rds",             emit: isgr_rds_updated

    script:
    """
    # Calculate GC content for regions
    samtools faidx ${params.genome_fasta} -r ${region} | \\
    awk 'BEGIN {IGNORECASE=1;a=0;c=0;g=0;t=0;}
         {if(NR==1)next;a+=gsub("A","");c+=gsub("C","");g+=gsub("G","");t+=gsub("T","");}
         />/ {print (c+g)/(a+c+g+t);a=0;c=0;g=0;t=0;next;}
         END {print (c+g)/(a+c+g+t);}' > 1kb.region.gc

    # Generate GC content analysis
    Rscript R/gcContent.R \
        --isgr ${isgr_rds} \
        --gc 1kb.region.gc \
        --out integrationSiteData.tsv \
        --plot GC_content.png \
        --update ISGR_gc.rds
    """
}

// Step 6: Clonality
process Clonality {
    tag "Clonality analysis"
    label 'process_low'
    conda 'environment.yml'

    input:
    path isgr_rds

    output:
    path "IS_clonality.png",   emit: is_clonality_plot
    path "clonalityData.tsv",  emit: clonality_data_tsv

    script:
    """
    Rscript R/clonality.R \
        --isgr ${isgr_rds} \
        --config ${params.config_file} \
        --plot IS_clonality.png \
        --out clonalityData.tsv
    """
}

// Step 7: Final report
process Report {
    tag "Report generation"
    label 'process_low'
    conda 'environment.yml'

    input:
    path logo_plot
    path gene_plot
    path dnase_plot
    path gc_content_plot
    path clonality_plot
    path report_rmd

    output:
    path "report.html", emit: html_report

    script:
    """
    # Copy plots to expected location if needed
    mkdir -p plots
    cp ${logo_plot} plots/
    cp ${gene_plot} plots/
    cp ${dnase_plot} plots/
    cp ${gc_content_plot} plots/
    cp ${clonality_plot} plots/

    # Render report
    Rscript -e "
    rmarkdown::render(
        '${report_rmd}', 
        output_file='report.html', 
        params = list(
            logo20bp_png = 'plots/${logo_plot.name}',
            IS_gene_png = 'plots/${gene_plot.name}',
            IS_dnase_png = 'plots/${dnase_plot.name}',
            GC_content_png = 'plots/${gc_content_plot.name}',
            IS_clonality_png = 'plots/${clonality_plot.name}'
        )
    )"
    """
}

// ----------------------
// Workflow definition
// ----------------------
workflow {
    // Validate inputs
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    genome_ch = Channel.fromPath(params.genome_fasta, checkIfExists: true)
    report_rmd_ch = Channel.fromPath(params.report_rmd, checkIfExists: true)

    // Step 1: Data preparation
    prep_out = DataReadInAndPreparation(input_ch)

    // Step 2: Sequence logo
    logo_out = SequenceLogo(
        prep_out.plus_region,
        prep_out.minus_region,
        prep_out.is_csv
    )

    // Step 3: Gene analysis
    gene_out = GeneAnalysis(logo_out.is_csv_updated)

    // Step 4: Open chromatin
    chrom_out = OpenChromatin(gene_out.isgr_rds)

    // Step 5: GC content
    gc_out = GCContent(
        chrom_out.isgr_rds_updated,
        chrom_out.region
    )

    // Step 6: Clonality
    clone_out = Clonality(gc_out.isgr_rds_updated)

    // Step 7: Generate report
    Report(
        logo_out.logo,
        gene_out.gene_plot,
        chrom_out.is_dnase_plot,
        gc_out.gc_plot,
        clone_out.is_clonality_plot,
        report_rmd_ch
    )

    // Publishing outputs
    logo_out.logo.subscribe { 
        it.copyTo("${params.outdir}/plots/${it.name}")
    }
    
    gene_out.gene_plot.subscribe { 
        it.copyTo("${params.outdir}/plots/${it.name}")
    }
    gene_out.genedata_tsv.subscribe { 
        it.copyTo("${params.outdir}/${it.name}")
    }
    
    chrom_out.is_dnase_plot.subscribe { 
        it.copyTo("${params.outdir}/plots/${it.name}")
    }
    
    gc_out.gc_plot.subscribe { 
        it.copyTo("${params.outdir}/plots/${it.name}")
    }
    gc_out.is_data_tsv.subscribe { 
        it.copyTo("${params.outdir}/${it.name}")
    }
    
    clone_out.is_clonality_plot.subscribe { 
        it.copyTo("${params.outdir}/plots/${it.name}")
    }
    clone_out.clonality_data_tsv.subscribe { 
        it.copyTo("${params.outdir}/${it.name}")
    }
    
    Report.out.html_report.subscribe { 
        it.copyTo("${params.outdir}/${it.name}")
    }
}

// ----------------------
// Process labels for resource allocation
// ----------------------
process {
    withLabel: 'process_low' {
        cpus = 1
        memory = 2.GB
        time = 1.h
    }
    withLabel: 'process_medium' {
        cpus = 2
        memory = 4.GB
        time = 2.h
    }
    withLabel: 'process_high' {
        cpus = 4
        memory = 8.GB
        time = 4.h
    }
}
