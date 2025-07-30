# GC content
# get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# check if a log file name is provided
if(length(args) == 0){
	log_file <- "log.txt"
}else{
	log_file <- args[1]
}

# open log file connection
log_connection <- file(log_file, open = "at")

# Redirect message and warings to the log file
sink(log_connection, type = "message")

suppressWarnings( {
	# source the libraries
	source("config/libraries.R")
	
	# Calculate GC content in 1 kb bins
	# Read IS GRanges object
	IS.GR <- readRDS(file = "intermediate/ISGR.rds")
	IS.GR$GC <- as.numeric(readLines("intermediate/1kb.region.gc")) * 100
	
	# save the IS table
	if (!file.exists("results/integrationSiteData.tsv")) {
		write.table(
			as.data.frame(IS.GR), "results/integrationSiteData.tsv",
			sep = "\t", col.names = T, row.names = F, quote = F, na = ""
		)
	}
	
	# Plot GC content
	genomeGCmean <- 40.89 # [Piovesan *et al*., 2019](https://pubmed.ncbi.nlm.nih.gov/30813969/)
	if (!file.exists("results/plots/GC_content.png")) {
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
		ggsave("results/plots/GC_content.png", plotObject, width = 12, height = 15, units = "cm", dpi = 300)
	}
	
	# Save IS GRanges object
	saveRDS(IS.GR, file = "intermediate/ISGR.rds")
})

# stop redirection and close file
sink(type = "message")
close(log_connection)
