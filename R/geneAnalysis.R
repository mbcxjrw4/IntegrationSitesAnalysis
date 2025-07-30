# gene analysis
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
	
	# Create IS GRanges object
	IS <- data.table::fread("intermediate/IS.csv")
	IS.GR <- unique(IS[, c("chromosome", "position", "sequence")])
	IS.GR <- GRanges(
		seqnames = IS.GR$chromosome,
		ranges = IRanges(start = IS.GR$position, end = IS.GR$position),
		sequence = IS.GR$sequence
	)
	
	# gencode.v22 gene info
	gencode.v22.genes <- read.table("resource/gencode.v22.gene.ranges")
	gencode.v22.genes <- GRanges(
		seqnames = gencode.v22.genes$V1,
	    ranges = IRanges(start = gencode.v22.genes$V2, end = gencode.v22.genes$V3),
		strand = gencode.v22.genes$V4, ensembl_id = gencode.v22.genes$V5,
		gene_name = gencode.v22.genes$V6, gene_type = gencode.v22.genes$V7
	)
	# calculate overlaps between IS and genes
	ISoverGenes <- findOverlaps( IS.GR, gencode.v22.genes, ignore.strand = T, select = "all" )
	IS.GR$gene_hit_count <- countQueryHits(ISoverGenes) # number of matched genes
	mcols(IS.GR)[ , c("ensembl_id","gene_name", "gene_type")] <- NA
	IS.GR$ensembl_id[queryHits(ISoverGenes)] <- gencode.v22.genes$ensembl_id[subjectHits(ISoverGenes)]
	IS.GR$gene_name[queryHits(ISoverGenes)] <- gencode.v22.genes$gene_name[subjectHits(ISoverGenes)]
	IS.GR$gene_type[queryHits(ISoverGenes)] <- gencode.v22.genes$gene_type[subjectHits(ISoverGenes)]
	# fewer gene groups
	IS.GR$gene_group <- ifelse(IS.GR$gene_type == "protein_coding", "Protein coding", "Other gene")
	IS.GR$gene_group[IS.GR$gene_hit_count == 0] <- "Not a gene"
	IS.GR$gene_group <- factor(IS.GR$gene_group, levels = rev(c("Protein coding", "Other gene", "Not a gene")))
	
	# Oncogenes
	# Oncoplex panel (2021-11-30) extracted from synonyms https://testguide.labmed.uw.edu/public/view/OPX
	# renamed MYCL1 to MYCL
	oncoPlex <- c(
		'ABL1','ABL2','ACVR1','AKT1','AKT2','AKT3','ALK','ANGPTL1','ANKRD26','APC','AR','ARAF','ARID1A','ARID1B','ASXL1','ASXL2',
		'ATM','ATR','ATRX','AURKA','AURKB','AXIN2','AXL','BABAM1','BAK1','BAP1','BARD1','BCL2','BCL2L11','BCOR','BCORL1','BCR',
		'BIRC3','BMPR1A','BRAF','BRCA1','BRCA2','BRIP1','BTK','C11orf95','CALR','CARD11','CBL','CBLB','CBLC','CCND1','CCND2',
		'CCNE1','CD19','CD274','CD33','CD74','CDC27','CDH1','CDK12','CDK4','CDK6','CDK8','CDK9','CDKN1A','CDKN2A','CDKN2B',
		'CEBPA','CHD1','CHEK1','CHEK2','CREBBP','CRLF2','CRX','CSF1R','CSF3R','CTCF','CTNNA1','CTNNB1','CUX1','DAXX','DDR2',
		'DDX41','DEPDC5','DICER1','DNAJB1','DNMT3A','DOCK7','DPYD','EBF1','EGFR','EIF3E','ELF1','EML4','EP300','EPAS1','EPCAM',
		'EPHA3','EPHA5','EPHB2','EPHB6','ERBB2','ERBB3','ERBB4','ERCC2','ERG','ESR1','ESR2','ETV6','EZH2','FAM175A','FANCA',
		'FANCB','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','FBXW7','FGFR1','FGFR2','FGFR3','FGFR4','FH','FKBP1A',
		'FLT1','FLT3','FLT4','FOXA1','FOXR2','GAB2','GALNT12','GATA1','GATA2','GATA3','GEN1','GLI1','GLTSCR1','GLTSCR2','GNA11',
		'GNAQ','GNAS','GREM1','GRIN2A','GRM3','H3F3A','H3F3B','HDAC4','HDAC9','HIF1A','HIST1H3B','HNF1A','HRAS','HSPH1','ID3',
		'IDH1','IDH2','IGF1R','IKZF1','IL7R','JAK1','JAK2','JAK3','KDM6A','KDR','KIF5B','KIT','KLF4','KMT2A','KMT2C','KMT2D',
		'KRAS','MAP2K1','MAP2K2','MAP2K4','MAPK1','MAX','MC1R','MCL1','MDM2','MDM4','MED12','MEGF6','MEN1','MET','MIOS','MITF',
		'MLH1','MLH3','MN1','MPL','MRE11A','MSH2','MSH6','MSLN','MTAP','MTOR','MUTYH','MYB','MYC','MYCL','MYCN','MYD88','MYOD1',
		'NAB2','NAT2','NBN','NF1','NF2','NKX2-1','NOTCH1','NOTCH2','NOTCH3','NPM1','NPRL2','NPRL3','NR4A3','NRAS','NT5C2','NTHL1',
		'NTRK1','NTRK2','NTRK3','NUDT15','PAK1','PALB2','PAX5','PBRM1','PDCD1LG2','PDGFRA','PDGFRB','PHF6','PHOX2B','PIK3CA',
		'PIK3CB','PIK3R1','PLCG2','PLK1','PLK2','PLK3','PLK4','PML','PMS2','POLD1','POLE','PPM1D','PRKAR1A','PRPF40B','PRPS1',
		'PTCH1','PTEN','PTPN11','PTPRD','QKI','RAC1','RAD21','RAD51B','RAD51C','RAD51D','RAF1','RARA','RB1','RECQL','RELA','RET',
		'RHEB','RICTOR','RINT1','RIT1','ROR1','ROS1','RPL10','RPS14','RPS15','RPS20','RPTOR','RRM1','RRM2','RSPO2','RSPO3',
		'RUNX1','SAMD9L','SDHA','SDHB','SDHC','SDHD','SETBP1','SETD2','SF1','SF3B1','SH2B3','SHH','SLX4','SMAD2','SMAD3','SMAD4',
		'SMARCA4','SMARCB1','SMC1A','SMC3','SMO','SPOP','SPRY4','SRC','SRP72','SRSF2','STAG2','STAT5B','STAT6','STK11','SUFU',
		'SUZ12','TACC3','TACSTD2','TCF3','TERC','TERT','TET1','TET2','TET3','TFE3','TFG','TGFBR2','TLX1','TMPRSS2','TP53','TP73',
		'TRAF7','TRRAP','TSC1','TSC2','TTYH1','TYMS','U2AF1','U2AF2','VHL','WRN','WT1','XRCC2','YAP1','ZBTB16','ZRSR2'
	)
	oncoPlex <- c(oncoPlex, "LMO2") # https://pubmed.ncbi.nlm.nih.gov/16084128/
	# counts of IS per gene
	gencode.v22.genes$IS_count <- countOverlaps(gencode.v22.genes, IS.GR, ignore.strand = T) # IS per gene
	gencode.v22.genes$IS_norm_count <- round(gencode.v22.genes$IS_count * 1000 / width(gencode.v22.genes), 2) # per kb of gene
	gencode.v22.genes$oncoPlex <- gencode.v22.genes$gene_name %in% oncoPlex
	gencode.v22.genes$oncoPlex <- ifelse(gencode.v22.genes$oncoPlex, "Oncogene", "Other")
	gencode.v22.genes$oncoPlex <- factor(gencode.v22.genes$oncoPlex, levels = c("Oncogene", "Other"))
	
	# Gene expression (ex-vivo CD3+ cells)
	fraietta_TPM <- read.table("resource/fraietta_TPM.tsv", sep = "\t", header = T, stringsAsFactors = F)
	
	# Oncogene comparisons within expression bins
	gencode.v22.geneDF <- merge(
		data.frame(mcols(gencode.v22.genes)), fraietta_TPM[, c("gene_name", "mean_TPM")],
		by = "gene_name", all.x = T, all.y = F, sort = F
	)
	gencode.v22.geneDF$expression <- cut(gencode.v22.geneDF$mean_TPM, breaks = c(-Inf, 1, 10, Inf), labels = c("Low expression", "Medium expression", "High expression"))
	
	# calculate p values using one-sided Kolmogorov-Smirnov test 
	# (alternative hypothesis: genes in OncoPlex panel have MORE integrations)
	KS_greater <- data.frame(
		expression = levels(gencode.v22.geneDF$expression),
		pval = as.numeric(
			by(gencode.v22.geneDF[, c("IS_norm_count", "oncoPlex")], gencode.v22.geneDF$expression, function(x) {
				suppressWarnings(ks.test(
					x = as.numeric(x[ x[, 2] == "Oncogene", 1 ]),
					y = as.numeric(x[x[ ,2] == "Other", 1]),
					alternative = "greater"
				)$p.value)
		})),
		stringsAsFactors = F
	)
	KS_greater$expression <- factor(KS_greater$expression, levels = unique(KS_greater$expression))
	# gene number data
	geneN <- data.frame(
		rbind(by(
			gencode.v22.geneDF[!is.na(gencode.v22.geneDF$mean_TPM), ],
			gencode.v22.geneDF[!is.na(gencode.v22.geneDF$mean_TPM), c("expression", "oncoPlex")],
			nrow
		))
	)
	geneN$expression <- rownames(geneN)
	geneN <- data.frame(
		tidyr::pivot_longer(geneN, cols = 1:2, names_to = "oncoPlex", values_to = "n"),
		stringsAsFactors = F
	)
	geneN$expression <- factor(geneN$expression, levels = levels(gencode.v22.geneDF$expression))
	
	# save the gene table
	if (!file.exists("results/geneData.tsv")) {
		write.table(
			gencode.v22.geneDF, "results/geneData.tsv",
			sep = "\t", col.names = T, row.names = F, quote = F, na = ""
		)
	}
	
	# Gene plotting
	if (!file.exists("results/plots/IS_gene.png")) {
		# Barplot of IS in/out of genes
		p1 <- ggplot(data.frame(mcols(IS.GR)), aes(x = "Gene", fill = gene_group)) +
			ggtitle("Integration site distribution") +
			geom_bar(position = "fill", colour = "grey30") +
			scale_y_continuous(breaks = 0.1*(0:10), labels = paste0(10*(0:10), "%")) +
			scale_x_discrete(expand = c(0,0)) +
			scale_fill_manual(values = adaptPalette[c(10, 8, 3)]) +
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
		# Boxplots of IS counts in expression bins (oncogene/not)
		p2 <- ggplot(gencode.v22.geneDF[!is.na(gencode.v22.geneDF$mean_TPM), ]) +
			ggtitle("Integration site frequencies within genes") +
			geom_jitter(
				aes(x = oncoPlex, y = IS_norm_count),
				colour = "grey70", size = 0.5, width = 0.35, shape = 16
			) +
			geom_boxplot(
				aes(x = oncoPlex, y = IS_norm_count, fill = oncoPlex),
				outlier.shape = NA, colour = "grey30"
			)  +
			geom_text(
				aes(x = 1.5, y = 17, label = paste0("p=",round(pval, 4))), data = KS_greater,
				hjust = 0.5, colour = "grey70", fontface = "italic", size = 3
			) +
			geom_text(
				aes(x = oncoPlex, y = 0, label = paste0("n=", n)), data = geneN,
				hjust = 0.5, vjust = 1.5, size = 3, colour = "grey40"
			) +
			geom_segment( aes(x=0.9, xend=2.1,y=16, yend=16), colour = "grey70", linewidth = 0.25 ) +
			geom_segment( aes(x=0.9, xend=0.9,y=15.5, yend=16), colour = "grey70", linewidth = 0.25 ) +
			geom_segment( aes(x=2.1, xend=2.1,y=16, yend=15.5), colour = "grey70", linewidth = 0.25 ) +
			facet_wrap(.~expression) +
			scale_y_continuous(
				trans = scales::pseudo_log_trans(base = 10),
				breaks = c(0:5, 5*(2:4))
			) +
			scale_fill_manual(values = adaptPalette[c(2, 1)]) +
			theme_minimal(base_size = 12) + 
			theme(
				strip.text = element_text(colour = "grey60", face = "bold", size = 10),
				panel.grid.major.x = element_blank(),
				panel.grid.minor = element_blank(),
				panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
				panel.border = element_rect(fill = NA, colour = "grey30"),
				axis.text.x = element_blank(),
				axis.text.y = element_text(colour = "grey50"),
				axis.title = element_text(colour = "grey40"),
				legend.position = "bottom",
				legend.margin=margin(0,0,0,0),
		        legend.box.margin=margin(-20,0,0,0),
				plot.title = element_text(hjust =0.5, colour ="grey40", size = 12)
			) +
			labs(x = "", y = "Normalised integration site count") +
			guides(fill = guide_legend(title = "Gene type:"))
		# combine
		plotObject <- p1 + plot_spacer() + p2 +
			plot_layout(widths = c(0.75, 0.25, 6))
		# save
		ggsave("results/plots/IS_gene.png", plotObject, width = 25, height = 15, units = "cm", dpi = 300)
	}
	
	# Save IS GRanges object
	saveRDS(IS.GR, file = "intermediate/ISGR.rds")
})

# stop redirection and close file
sink(type = "message")
close(log_connection)
