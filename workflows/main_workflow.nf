#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

processOptions = { cpus = 1; memory = '4 GB' }

workflow {

    /*
     Parameters expected:
       params.input_intsites        -> input/intSites.tsv
       params.genome_fa            -> path to genome fasta ($HOME/path_to_genome/GRCh38.p13.genome.fa)
       params.intermediate_dir     -> intermediate/
       params.resource_dir         -> resource/
       params.output_dir           -> results/
       params.patient              -> patient id (optional)
    */

    // Create channels for the single primary input file
    Channel
        .fromPath( params.input_intsites ?: 'input/intSites.tsv' )
        .set { intsites_ch }

    /*
      Step 1: Data read in & preparation
      - Runs R/dataReadInAndPreparation.R
      - Produces:
        intermediate/20bp.plus.region
        intermediate/20bp.minus.region
        intermediate/IS.csv
    */
    process DataPrep {
        tag "DataPrep"

        publishDir "${params.output_dir}/intermediate", mode: 'copy'

        conda 'environment.yml'

        input:
        path intsites_file

        output:
        path "20bp.plus.region", emit: plus_region
        path "20bp.minus.region", emit: minus_region
        path "IS.csv", emit: is_csv

        script:
        """
        mkdir -p ${params.intermediate_dir ?: 'intermediate'}
        Rscript R/dataReadInAndPreparation.R \
          --patient '${params.patient ?: ""}' \
          --input ${intsites_file} \
          --outdir ${params.intermediate_dir ?: 'intermediate'}
        # expect script to write 20bp.plus.region 20bp.minus.region IS.csv into intermediate/
        ln -sfn ${params.intermediate_dir}/20bp.plus.region 20bp.plus.region
        ln -sfn ${params.intermediate_dir}/20bp.minus.region 20bp.minus.region
        ln -sfn ${params.intermediate_dir}/IS.csv IS.csv
        """
    }

    /*
      Step 2A & 2B: samtools getfa for plus and minus regions
      Produces: intermediate/20bp.plus.fa and intermediate/20bp.minus.fa
    */
    process SamtoolsGetFasta {
        tag "$region_file.simpleName"
        publishDir "${params.output_dir}/intermediate", mode: 'copy'

        conda 'environment.yml'

        input:
        path region_file
        val genome_fa = params.genome_fa

        output:
        path "*.fa", emit: fasta_out

        script:
        """
        mkdir -p ${params.intermediate_dir ?: 'intermediate'}
        # Samtools faidx -r expects region file format (chr:start-end or BED depending on samtools build)
        samtools faidx ${genome_fa} -r ${region_file} > ${region_file.simpleName}.fa
        """
    }

    /*
      Step 2C: Sequence logo R script
      - Reads: 20bp.plus.fa, 20bp.minus.fa, IS.csv
      - Writes: results/plots/logo20bp.png and may update IS.csv
    */
    process SequenceLogo {
        tag "SequenceLogo"
        publishDir "${params.output_dir}/plots", mode: 'copy'

        conda 'environment.yml'

        input:
        path plus_fa
        path minus_fa
        path is_csv

        output:
        path "logo20bp.png"
        path "IS.csv"   // if updated
        script:
        """
        mkdir -p ${params.output_dir}/plots
        Rscript R/informationContent.R \
          --plus ${plus_fa} --minus ${minus_fa} --is ${is_csv} \
          --outplot ${params.output_dir}/plots/logo20bp.png \
          --outis ${params.intermediate_dir ?: 'intermediate'}/IS.csv
        cp ${params.intermediate_dir ?: 'intermediate'}/IS.csv IS.csv || true
        cp ${params.output_dir}/plots/logo20bp.png logo20bp.png
        """
    }

    /*
      Step 3: Gene analysis
      Reads: intermediate/IS.csv, resource/gencode.v22.gene.ranges
      Writes: results/plots/IS_gene.png and intermediate/ISGR.rds
    */
    process GeneAnalysis {
        tag "GeneAnalysis"
        publishDir "${params.output_dir}/plots", mode: 'copy'
        publishDir "${params.output_dir}/intermediate", mode: 'copy'

        conda 'environment.yml'

        input:
        path is_csv
        path gene_ranges from file("${params.resource_dir ?: 'resource'}/gencode.v22.gene.ranges")
        path tpm from file("${params.resource_dir ?: 'resource'}/fraietta_TPM.tsv")

        output:
        path "ISGR.rds", emit: isgr_rds
        path "IS_gene.png"

        script:
        """
        mkdir -p ${params.output_dir}/plots
        Rscript R/geneAnalysis.R \
          --is ${is_csv} \
          --gene_ranges ${gene_ranges} \
          --tpm ${tpm} \
          --out_plot ${params.output_dir}/plots/IS_gene.png \
          --out_rds ${params.output_dir}/intermediate/ISGR.rds
        cp ${params.output_dir}/plots/IS_gene.png IS_gene.png
        cp ${params.output_dir}/intermediate/ISGR.rds ISGR.rds
        """
    }

    /*
      Step 4: Open chromatin analysis
      Reads: intermediate/ISGR.rds, resource/wgEncodeRegDnaseClustered.txt.gz
      Writes: results/plots/IS_dnase.png, and updates intermediate/ISGR.rds and creates intermediate/1kb.region
    */
    process OpenChromatin {
        tag "OpenChromatin"
        publishDir "${params.output_dir}/plots", mode: 'copy'
        publishDir "${params.output_dir}/intermediate", mode: 'copy'

        conda 'environment.yml'

        input:
        path isgr_rds
        path dnase_file from file("${params.resource_dir ?: 'resource'}/wgEncodeRegDnaseClustered.txt.gz")

        output:
        path "ISGR.rds", emit: isgr_rds_out
        path "1kb.region", emit: onekb_region
        path "IS_dnase.png"

        script:
        """
        mkdir -p ${params.output_dir}/plots
        Rscript R/openChromatin.R \
          --isgr ${isgr_rds} \
          --dnase ${dnase_file} \
          --out_rds ${params.output_dir}/intermediate/ISGR.rds \
          --out_region ${params.output_dir}/intermediate/1kb.region \
          --out_plot ${params.output_dir}/plots/IS_dnase.png

        cp ${params.output_dir}/plots/IS_dnase.png IS_dnase.png
        cp ${params.output_dir}/intermediate/ISGR.rds ISGR.rds
        cp ${params.output_dir}/intermediate/1kb.region 1kb.region
        """
    }

    /*
      Step 5A: GC content extraction via samtools faidx + awk -> intermediate/1kb.region.gc
    */
    process GCFaidxAndAWK {
        tag "GC_faidx"
        publishDir "${params.output_dir}/intermediate", mode: 'copy'
        conda 'environment.yml'

        input:
        path region_file
        val genome_fa = params.genome_fa

        output:
        path "1kb.region.gc", emit: gc_out

        script:
        """
        # region_file expected to be in ${params.intermediate_dir}/1kb.region or linked in work dir
        samtools faidx ${genome_fa} -r ${region_file} > 1kb.region.fa

        # Compute GC content per sequence using awk (as in original)
        awk 'BEGIN { IGNORECASE=1; a=0; c=0; g=0;t=0; } \
            { if (NR==1) next; a += gsub("A",""); c += gsub("C",""); g += gsub("G",""); t += gsub("T",""); } \
            />/ { if (a+c+g+t>0) print (c+g)/(a+c+g+t); else print "NA"; a=0; c=0; g=0; t=0; next; } \
            END { if (a+c+g+t>0) print (c+g)/(a+c+g+t); }' 1kb.region.fa > 1kb.region.gc

        cp 1kb.region.gc ${params.output_dir}/intermediate/1kb.region.gc
        """
    }

    /*
      Step 5B: GC content R script
      Reads: ISGR.rds and 1kb.region.gc
      Writes: results/integrationSiteData.tsv, results/plots/GC_content.png and updated ISGR.rds
    */
    process GCContentR {
        tag "GCContentR"
        publishDir "${params.output_dir}", mode: 'copy'
        publishDir "${params.output_dir}/plots", mode: 'copy'
        publishDir "${params.output_dir}/intermediate", mode: 'copy'

        conda 'environment.yml'

        input:
        path isgr_rds
        path gc_file

        output:
        path "integrationSiteData.tsv"
        path "GC_content.png"
        path "ISGR.rds"

        script:
        """
        Rscript R/gcContent.R \
          --isgr ${isgr_rds} \
          --gc ${gc_file} \
          --out_tsv ${params.output_dir}/integrationSiteData.tsv \
          --out_plot ${params.output_dir}/plots/GC_content.png \
          --out_rds ${params.output_dir}/intermediate/ISGR.rds

        cp ${params.output_dir}/integrationSiteData.tsv .
        cp ${params.output_dir}/plots/GC_content.png GC_content.png
        cp ${params.output_dir}/intermediate/ISGR.rds ISGR.rds
        """
    }

    /*
      Step 6: Clonality analysis
      Reads: ISGR.rds and resource/GRCh38.p13_1Mb_geneTPM.bed
      Writes: results/plots/IS_clonality.png and results/clonalityData.tsv
    */
    process Clonality {
        tag "Clonality"
        publishDir "${params.output_dir}/plots", mode: 'copy'
        publishDir "${params.output_dir}", mode: 'copy'

        conda 'environment.yml'

        input:
        path isgr_rds
        path geneTPM from file("${params.resource_dir ?: 'resource'}/GRCh38.p13_1Mb_geneTPM.bed")

        output:
        path "IS_clonality.png"
        path "clonalityData.tsv"

        script:
        """
        Rscript R/clonality.R \
          --isgr ${isgr_rds} \
          --gene_tpm ${geneTPM} \
          --out_plot ${params.output_dir}/plots/IS_clonality.png \
          --out_tsv ${params.output_dir}/clonalityData.tsv

        cp ${params.output_dir}/plots/IS_clonality.png IS_clonality.png
        cp ${params.output_dir}/clonalityData.tsv clonalityData.tsv
        """
    }

    /*
      Step 7: Final report.Rmd -> report.html
      Inputs: many files - we wire required outputs
    */
    process Report {
        tag "Report"
        publishDir "${params.output_dir}", mode: 'copy'

        conda 'environment.yml'

        input:
        path integrationSiteData  from file("${params.output_dir}/integrationSiteData.tsv")
        path logoPlot            from file("${params.output_dir}/plots/logo20bp.png")
        path genePlot            from file("${params.output_dir}/plots/IS_gene.png")
        path dnasePlot           from file("${params.output_dir}/plots/IS_dnase.png")
        path gcPlot              from file("${params.output_dir}/plots/GC_content.png")
        path clonalityPlot       from file("${params.output_dir}/plots/IS_clonality.png")
        path clonalityTsv        from file("${params.output_dir}/clonalityData.tsv")

        output:
        path "report.html"

        script:
        """
        Rscript -e "rmarkdown::render('report.Rmd', output_file='report.html',
                 params=list(
                   integration_tsv='${integrationSiteData}',
                   logo='${logoPlot}',
                   gene_plot='${genePlot}',
                   dnase_plot='${dnasePlot}',
                   gc_plot='${gcPlot}',
                   clonality_plot='${clonalityPlot}',
                   clonality_tsv='${clonalityTsv}'
                 ))"
        cp report.html ${params.output_dir}/report.html || true
        """
    }

    // =======================
    // Workflow wiring
    // =======================

    // Step 1: Data prep
    data_prep = intsites_ch | DataPrep

    // Step 2: samtools for plus/minus
    plus_region = data_prep.plus_region
    minus_region = data_prep.minus_region
    sam_plus = plus_region | SamtoolsGetFasta
    sam_minus = minus_region | SamtoolsGetFasta

    // Step 2: Sequence logo (wait for both fasta files and IS.csv)
    sequence_logo = Channel
        .combine( sam_plus.fasta_out.flatten(), sam_minus.fasta_out.flatten(), data_prep.is_csv )
        | SequenceLogo

    // Step 3: Gene analysis uses possibly updated IS.csv from SequenceLogo
    gene_analysis = sequence_logo.IS.csv.flatten() | GeneAnalysis

    // Step 4: Open chromatin (uses ISGR.rds)
    open_chrom = gene_analysis.isgr_rds | OpenChromatin

    // Step 5A: GC faidx + awk on generated 1kb.region
    gc_faidx = open_chrom.onekb_region | GCFaidxAndAWK

    // Step 5B: GC content R
    gc_content_r = open_chrom.isgr_rds.mix( gc_faidx.gc_out ) | GCContentR

    // Step 6: Clonality
    clonality = gc_content_r.isgr_rds | Clonality

    // Step 7: Report
    report = Channel
        .from([ params.output_dir + '/integrationSiteData.tsv' ])
        .combine(
            Channel.fromPath("${params.output_dir}/plots/logo20bp.png"),
            Channel.fromPath("${params.output_dir}/plots/IS_gene.png"),
            Channel.fromPath("${params.output_dir}/plots/IS_dnase.png"),
            Channel.fromPath("${params.output_dir}/plots/GC_content.png"),
            Channel.fromPath("${params.output_dir}/plots/IS_clonality.png"),
            Channel.fromPath("${params.output_dir}/clonalityData.tsv")
        )
        | Report
}
