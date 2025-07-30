# Clonality
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
	
	# Read IS GRanges object
	IS.GR <- readRDS(file = "intermediate/ISGR.rds")
	
	# Read TPM sums and count IS
	hg38_1mb <- read.table("resource/GRCh38.p13_1Mb_geneTPM.bed", na = ".")
	hg38_1mb$V4[is.na(hg38_1mb$V4)] <- 0
	hg38_1mb <- GRanges(
		seqnames = hg38_1mb$V1,
		ranges = IRanges(start = hg38_1mb$V2, end = hg38_1mb$V3),
		TPM_sum = hg38_1mb$V4
	)
	# count IS in 1 Mb bins
	hg38_1mb$IS_count <- countOverlaps(hg38_1mb, IS.GR, ignore.strand = T)
	
	# Plot clonality vs expression
	if (!file.exists("results/plots/IS_clonality.png")) {
		p1 <- ggplot(data.frame(TPM_sum = hg38_1mb$TPM_sum, stringsAsFactors = F)) +
			  ggtitle("IS clonality correlates with gene expression") +
			  geom_density(aes(x = TPM_sum), fill = "#5062AF33", colour = "#5062AF88", size = 0.5) +
			  scale_x_continuous(
			  	trans = scales::pseudo_log_trans(base = 10),
			  	limits = c(NA, 12000), expand = c(0.01,0.01)
			  ) +
			  theme_void() +
			  theme(
				plot.title = element_text(hjust =0.5, vjust = 5, colour ="grey40", size = 12, face="bold")
			  )
	
		p2 <- ggplot(data.frame(mcols(hg38_1mb)), aes(x = TPM_sum, y = IS_count)) +
			  geom_point(colour = "#5062AF", size = 0.75) +
			  scale_x_continuous(
			  	name = "Sum of overlapping genes' expression, TPM",
			  	trans = scales::pseudo_log_trans(base = 10),
			  	breaks = 10^(0:5), limits = c(NA, 12000), expand = c(0.01,0.01)
			  ) +
			  scale_y_continuous(
			  	name = "Integration site count",
			  	trans = scales::pseudo_log_trans(base = 10),
			  	breaks = 10^(0:3), limits = c(NA, 1000), expand = c(0.01,0.01)
			  ) +
			  theme_minimal(base_size = 12) +
			  theme(
			  	panel.grid.minor = element_blank(),
			  	panel.grid.major = element_line(linetype = "dashed", colour = "grey80"),
			  	panel.border = element_rect(fill = NA, colour = "grey30"),
			  	axis.text = element_text(colour = "grey50"),
			  	axis.title = element_text(colour = "grey40"),
			  	plot.title = element_text(hjust =0.5, colour ="grey40", size = 12, face="bold")
			  )
	
		p3 <- ggplot(data.frame(IS_count = hg38_1mb$IS_count, stringsAsFactors = F)) +
			  geom_density(aes(y = IS_count), fill = "#5062AF33", colour = "#5062AF88", size =0.5) +
			  scale_y_continuous(
			  	trans = scales::pseudo_log_trans(base = 10),
			  	limits = c(NA, 1000), expand = c(0.01,0.01)
			  ) +
			  theme_void()
		
		plotObject <- p1 + plot_spacer() + p2 + p3 +
					  plot_layout(
					  	widths = c(0.9, 0.1), heights = c(0.1, 0.9),
					  	ncol = 2, nrow = 2
					  )
	
		ggsave("results/plots/IS_clonality.png", plotObject, width = 18, height = 15, units = "cm", dpi = 300)
	}
	
	if (!file.exists("results/clonalityData.tsv")) {
		write.table(
			as.data.frame(hg38_1mb), "results/clonalityData.tsv",
			sep = "\t", col.names = T, row.names = F, quote = F, na = ""
		)
	}
})

# stop redirection and close file
sink(type = "message")
close(log_connection)
