#!/usr/bin/env Rscript
# Open chromatin

suppressPackageStartupMessages(library(optparse))

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default="config/config.yaml"),
  make_option("--isgr", type="character", help="ISGR RDS"),
  make_option("--plot", type="character", help="Output plot"),
  make_option("--out", type="character", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# ---------------------------
# Load YAML config
# ---------------------------
source("config/load_config.R")  # defines `cfg`

message("Running open chromatin analysis...")

if (!is.null(opt$isgr)) {
  message("Input file: ", opt$isgr)
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
# Read IS GRanges object
	IS.GR <- readRDS(file = opt$isgr)
	
# Read DNAseI clusters data
dnase <- read.table(cfg$dnase)
# table schema: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=wgEncodeRegDnaseClustered&hgta_table=wgEncodeRegDnaseClustered&hgta_doSchema=describe+table+schema
dnase <- GRanges(
	seqnames = dnase$V2,
    ranges = IRanges(start = dnase$V3, end = dnase$V4)
)
dnaseHits <- distanceToNearest(IS.GR, dnase, ignore.strand = T)
	
IS.GR$dnase_dist <- mcols(dnaseHits)$distance
IS.GR$dnase <- ifelse( IS.GR$dnase_dist <= 1000, "Within 1kb", "Further than 1kb")
IS.GR$dnase[IS.GR$dnase_dist == 0] <- "DNAseI HS"
IS.GR$dnase <- factor(IS.GR$dnase, levels = c("Further than 1kb", "Within 1kb", "DNAseI HS" ))
	
if (!file.exists(opt$plot)) {
	# plot
	plotObject <- ggplot(data.frame(mcols(IS.GR)), aes(x = "DNAseI", fill = dnase)) +
				  ggtitle("IS distribution: open chromatin") +
				  geom_bar(position = "fill", colour = "grey30") +
				  scale_y_continuous(breaks = 0.1*(0:10), labels = paste0(10*(0:10), "%")) +
				  scale_x_discrete(expand = c(0,0)) +
				  scale_fill_manual(values = adaptPalette[c(10, 1, 7)]) +
				  theme_minimal(base_size = 12) +
				  theme(
					  panel.grid = element_blank(),
					  axis.ticks = element_line(colour = "grey30"),
					  axis.ticks.length.x = unit(0, "cm"),
					  axis.text.x = element_blank(),
					  axis.text.y = element_text(colour = "grey50"),
					  axis.title = element_text(colour = "grey40"),
					  legend.margin=margin(0,0,0,0),
			          legend.box.margin=margin(-250,0,0,-5),
					  plot.title = element_text(hjust =0.5, colour ="grey40", size = 12)
				  ) +
				  labs(x = "", y = "Integration site proportion") +
				  guides(fill = guide_legend(title = ""))
	# save the plot
	ggsave(opt$plot, plotObject, width = 8, height = 15, units = "cm", dpi = 300)
}
	
# Save IS GRanges object
saveRDS(IS.GR, file = opt$isgr)
	
# Calculate GC content in 1 kb bins
if (!file.exists(file.path(opt$out, "1kb.region"))) {
	regions.1kb <- paste0(as.character(seqnames(IS.GR)), ":", as.integer(start(IS.GR))-499, "-",  as.integer(end(IS.GR))+500 )
	writeLines(regions.1kb, file.path(opt$out, "1kb.region"))
}

message("Open chromatin analysis complete ✅")
