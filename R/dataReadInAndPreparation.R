# data read in and preparation
IS <- data.table::fread("input/intSites.tsv")

IS <- IS[!grepl("_", IS$chromosome), ] # remove alternative and random sequences
source("config/variables.R")
IS <- IS[IS$subject %in% patientId, ]

# Create 20 bp region files
IS[, "regions"] <- apply(IS, 1, function(x) paste0(x[1], ":", as.integer(x[2])-10, "-",  as.integer(x[2])+10))
IS$regions <- sapply(IS$regions, as.character)

if (!dir.exists("intermediate")) dir.create("intermediate")
if (!file.exists("intermediate/20bp.plus.region")) {
	writeLines(sapply(unique(IS[IS$strand == "+", "regions"]), as.character), "intermediate/20bp.plus.region")
}
if (!file.exists("intermediate/20bp.minus.region")) {
	writeLines(sapply(unique(IS[IS$strand == "-", "regions"]), as.character), "intermediate/20bp.minus.region")
}

# update IS file
data.table::fwrite(IS, file = "intermediate/IS.csv")
