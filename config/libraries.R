# libraries
suppressPackageStartupMessages( {
	library(ggplot2)
	library(ggseqlogo)
	library(patchwork)
	library(GenomicRanges)
	library(GenomicFeatures)
	library(rtracklayer)
	library(openxlsx)
    library(data.table)
    library(gridExtra)
})

# variables
adaptPalette <- c(
  "#5062AF","#BF1E2D","#6D4685","#62AF50","#E1C000","#1EBFB0","#38457B","#AA86BF","#315728","#969696",
  "#29066B","#AF4BCE","#EB548C","#8D2D56","#FFA500","#1DE4BD","#820401","#DE542C","#EABD3B","#005D68",
  "#8A9A5B","#C7EA46","#50C878","#98FB98","#555555","#1AC9E6")
