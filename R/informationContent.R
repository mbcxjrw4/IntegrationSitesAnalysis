#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default="config/config.yaml"),
  make_option("--plus", type="character", help="20bp plus fasta"),
  make_option("--minus", type="character", help="20bp minus fasta"),
  make_option("--is", type="character", help="Integration site CSV"),
  make_option("--output", type="character", help="Output plot PNG")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# ---------------------------
# Load YAML config
# ---------------------------
source("config/load_config.R")  # defines `cfg`

# ---------------------------
# Example usage
# ---------------------------
message("Running information content analysis...")
message("Config file: ", opt$config)
message("Patient ID: ", cfg$patient_id)

if (!is.null(opt$plus)) {
  message("Input file: ", opt$plus)
}

if (!is.null(opt$minus)) {
  message("Input file: ", opt$minus)
}

if (!is.null(opt$is)) {
  message("Input file: ", opt$is)
}

if (!is.null(opt$output)) {
  message("Output file/dir: ", opt$output)
}

# ---------------------------
# Main analysis logic
# ---------------------------
	
fasta <- readLines(opt$plus)
fasta <- c(fasta, readLines(opt$minus))
fasta <- t(matrix(fasta, nrow = 2)) # make long table of names and sequences
fasta[, 1] <- gsub("^>", "", fasta[, 1])
fasta <- cbind(
	gsub( ":.+$", "", fasta[ , 1] ),
	as.integer(gsub( "^.+:|-.+$", "", fasta[ , 1] )) + 10,
	fasta[ , 2]
) # extract original IS coordinate
fasta <- data.frame(fasta, stringsAsFactors = F)
names(fasta) <- c("chromosome", "position", "sequence")
	
fasta$position <- as.numeric(fasta$position)
	
IS <- data.table::fread(opt$is)
	
IS <- merge(IS, fasta, by = c("chromosome", "position"), all.x = T, all.y =F, sort = F)
IS$sequence <- substring(IS$sequence, 1, 20) # trim to 20 bp
	
#IS$timePoint <- factor(IS$timePoint, levels = c("M18", "M36", "M24"))
	
# Sequence logo plots	
if (!file.exists(opt$output)) {
	# Initialize an empty list to store individual plots
  	seqPlot <- list()
		
  	# Loop over each timepoint for the subject
  	for (i in 1:length(unique(IS$timePoint))) {
  		timepoint <- unique(IS$timePoint)[i]
	  	  
  		# Generate plot for the current timepoint
  		p <- ggplot(IS[IS$timePoint == timepoint, ]) + 
  	    	 ggtitle(paste("Logo Plot for", unique(IS$subject), "at Timepoint", timepoint)) +
  	    	 geom_logo(na.omit(IS$sequence), seq_type = "dna") + theme_logo() +
  	    	 scale_x_discrete(limits = factor(1:20), labels = -10:9) +
  	    	 scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2)) +
  	    	 theme(
				 plot.background = element_rect(colour = 'grey75'),
				 axis.line = element_line(colour="grey40", linewidth=0.75, linetype="solid"),
	  	    	 axis.ticks = element_line(colour="grey40", linewidth = 0.25, linetype="solid"),
	  	    	 axis.text.x = element_text(angle = 90, hjust = 1),
	  	    	 plot.title = element_text(hjust = 0.5, colour = "grey25")
			 )
		
		# Store the plot in the list with a unique identifier
	  	seqPlot[[paste("Timepoint", timepoint)]] <- p
	}
	  	
	# Combine plots with patchwork (if necessary)
	plotObject <- seqPlot[[1]]
	if (length(seqPlot) > 1) {
		for (j in 2:length(seqPlot)) {
			plotObject <- plotObject / seqPlot[[j]]
	  	}
	 }
	  
	# Save the combined plot
	ggsave(opt$output, plotObject, width = 30, height = 25, units = "cm")
}
	
# update IS file
data.table::fwrite(IS, file = opt$is)

message("Information content analysis complete ✅")
