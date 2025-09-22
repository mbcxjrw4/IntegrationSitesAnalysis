#!/usr/bin/env Rscript
# GC content

suppressPackageStartupMessages(library(optparse))

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default="config/config.yaml"),
  make_option("--isgr", type="character", help="ISGR RDS"),
  make_option("--gc", type="character", help="GC content file"),
  make_option("--out", type="character", help="Output TSV"),
  make_option("--plot", type="character", help="Output plot path"),
  make_option("--update", type="character", help="Updated ISGR RDS")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# ---------------------------
# Load YAML config
# ---------------------------
source("config/load_config.R")  # defines `cfg`

message("Running GC content analysis...")

if (!is.null(opt$isgr)) {
  message("Input file: ", opt$isgr)
}
if (!is.null(opt$gc)) {
  message("Input file: ", opt$gc)
}
if (!is.null(opt$out)) {
  message("Output file/dir: ", opt$out)
}
if (!is.null(opt$plot)) {
  message("Plot path: ", opt$plot)
}

# ---------------------------
# Main analysis logic
# ---------------------------
# Calculate GC content in 1 kb bins
# Read IS GRanges object
IS.GR <- readRDS(file = opt$isgr)
IS.GR$GC <- as.numeric(readLines(opt$gc)) * 100
	
# save the IS table
if (!file.exists(opt$out)) {
	write.table(
		as.data.frame(IS.GR), opt$out,
		sep = "\t", col.names = T, row.names = F, quote = F, na = ""
	)
}
	
# Plot GC content
genomeGCmean <- 40.89 # [Piovesan *et al*., 2019](https://pubmed.ncbi.nlm.nih.gov/30813969/)
if (!file.exists(opt$plot)) {
	p1 <- ggplot(data.frame(GC = IS.GR$GC, stringsAsFactors = F), aes(x = "GC content", y = GC)) +
		  geom_jitter(size = 0.2, colour  = "#96969633", width = 0.38) +
		  geom_hline(yintercept = genomeGCmean, colour = "#8F1622", linewidth = 1.2, linetype = "dashed") +
		  geom_boxplot(outlier.shape = NA, colour = "#5062AF", size = 1.2) +
		  scale_y_continuous(limits = c(20, 80), breaks = 10*(2:8), labels = paste0(10*(2:8), "%")) +
		  theme_minimal(base_size = 12) +
		  theme(
			  panel.grid.major.x = element_blank(),
			  panel.grid.minor = element_blank(),
			  panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
			  panel.border = element_rect(fill = NA, colour = "grey30"),
			  axis.ticks = element_line(colour = "grey30"),
			  axis.ticks.length.x = unit(0, "cm"),
			  axis.text.x = element_blank(),
			  axis.title.x = element_blank(),
			  axis.text.y = element_text(colour = "grey50"),
			  axis.title.y = element_text(colour = "grey40"),
			  plot.title = element_text(hjust =0.5, colour ="grey40", size = 12)
		  )
	p2  <- ggplot(data.frame(GC = IS.GR$GC, stringsAsFactors = F)) +
		   geom_density(aes(y = GC), fill = "#5062AF55", colour = "#5062AF", size =1.4) +
		   geom_hline(yintercept = genomeGCmean, colour = "#8F1622", linewidth = 1.2, linetype = "dashed") +
		   scale_y_continuous(limits = c(20, 80)) +
		   theme_void()
	plotObject <- p1 + p2
	ggsave(opt$plot, plotObject, width = 12, height = 15, units = "cm", dpi = 300)
}
	
# Save IS GRanges object
saveRDS(IS.GR, file = opt$update)

message("GC content analysis complete ✅")
