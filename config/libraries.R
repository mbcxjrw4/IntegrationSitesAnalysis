# config/libraries.R

required_packages <- c(
  "ggplot2", "ggseqlogo", "data.table", "patchwork", "GenomicRanges", "GenomicFeatures", "rtracklayer",
  "openxlsx", "gridExtra"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Missing required package:", pkg))
  }
  library(pkg, character.only = TRUE)
}
