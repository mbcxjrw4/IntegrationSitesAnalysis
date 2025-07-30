# Information content
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

suppressWarnings({
	# source the libraries
	source("config/libraries.R")
	
	fasta <- readLines("intermediate/20bp.plus.fa")
	fasta <- c(fasta, readLines("intermediate/20bp.minus.fa"))
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
	
	IS <- data.table::fread("intermediate/IS.csv")
	
	IS <- merge(IS, fasta, by = c("chromosome", "position"), all.x = T, all.y =F, sort = F)
	IS$sequence <- substring(IS$sequence, 1, 20) # trim to 20 bp
	
	#IS$timePoint <- factor(IS$timePoint, levels = c("M18", "M36", "M24"))
	
	# Sequence logo plots
	if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)
	
	if (!file.exists("results/plots/logo20bp.png")) {
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
	  	ggsave("results/plots/logo20bp.png", plotObject, width = 30, height = 25, units = "cm")
	}
	
	# update IS file
	data.table::fwrite(IS, file = "intermediate/IS.csv")
})

# stop redirection and close file
sink(type = "message")
close(log_connection)
