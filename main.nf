#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
  Nextflow pipeline for IntegrationSitesAnalysis
  Mirrors the steps provided in the bash pipeline.
*/

include { workflow as wf } from './workflows/main_workflow.nf'
wf()
