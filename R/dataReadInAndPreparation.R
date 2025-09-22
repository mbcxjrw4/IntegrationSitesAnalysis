#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))

# ---------------------------
# CLI options
# ---------------------------
option_list <- list(
  make_option(c("-i", "--input"),  type="character", help="Input integration sites file (.tsv)", metavar="file"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory for intermediate files", metavar="dir"),
  make_option(c("-c", "--config"), type="character", default="config/config.yaml", help="Path to YAML config file", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Missing required arguments: --input and/or --outdir", call.=FALSE)
}

# ---------------------------
# Load config
# ---------------------------
source("config/load_config.R")  # this defines global variable `cfg`

# ---------------------------
# Main script
# ---------------------------
message("Running Data Preparation...")
message("Input file: ", opt$input)
message("Output directory: ", opt$outdir)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Use config values
patient_id <- cfg$patient_id
message("Patient ID from config: ", patient_id)

# data read in and preparation
IS <- data.table::fread(opt$input)

IS <- IS[!grepl("_", IS$chromosome), ] # remove alternative and random sequences
IS <- IS[IS$subject %in% patient_id, ]

# Create 20 bp region files
IS[, "regions"] <- apply(IS, 1, function(x) paste0(x[1], ":", as.integer(x[2])-10, "-",  as.integer(x[2])+10))
IS$regions <- sapply(IS$regions, as.character)

if (!dir.exists("intermediate")) dir.create("intermediate")
if (!file.exists("intermediate/20bp.plus.region")) {
	writeLines(sapply(unique(IS[IS$strand == "+", "regions"]), as.character), file.path(opt$outdir, "20bp.plus.region"))
}
if (!file.exists("intermediate/20bp.minus.region")) {
	writeLines(sapply(unique(IS[IS$strand == "-", "regions"]), as.character), file.path(opt$outdir, "20bp.minus.region"))
}

# update IS file
data.table::fwrite(IS, file = file.path(opt$outdir, "IS.csv"))

message("Data preparation complete ✅")
