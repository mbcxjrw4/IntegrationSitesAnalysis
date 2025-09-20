#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(optparse))

# Function to load YAML config
load_config <- function(config_file = "config/config.yaml") {
  cfg <- yaml::read_yaml(config_file)
  return(cfg)
}

# Command line option
option_list <- list(
  make_option(c("-c", "--config"), type="character", default="config/config.yaml",
              help="Path to YAML config file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cfg <- load_config(opt$config)

# Make cfg available globally
assign("cfg", cfg, envir = .GlobalEnv)

# Example usage inside analysis scripts:
# patient_id <- cfg$patient_id
# genome_fasta <- cfg$genome_fasta
