# Defined as a function, return value is table of differential fitness results (source into R/RStudio for interactive analyses)

TnSeqDESeq <- function(ctrl_pfx, ctrl_reps, test_pfx, test_reps, gff_pfx, out_pfx, to_trim, in_files) {
	# Read in sites files
	library(dplyr)
	sites <- data.frame(Pos=c(0)) 
	for (i in 1:length(in_files)) {
		newsites <- read.table(paste(paste(in_files[i], in_files[i], sep="/"), "sites.txt", sep="-")) 
		colnames(newsites) <- c(paste("V", i, sep=""), "Pos")
		newsites <- tail(newsites, n=-to_trim) %>% arrange(Pos)
		sites <- merge(sites, newsites, all=T) 
	}
	sites <- tail(sites, n=-1)
	sites[is.na(sites)] <- 0
	
	# OPTIONAL - perform site filtering. Example: only consider sites identified in 2 or more of 6 replicates
	#sites <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0) + (V3 > 0) + (V4 > 0)) %>% filter(numreps >= 2)
	#sites <- sites[-6]
	sites  <- sites %>% mutate(numreps = (V1 > 0 ) + (V2 > 0) + (V3 > 0) + (V4 > 0) + (V5 > 0) + (V6 > 0)) %>% filter(numreps >= 2)
	sites <- sites[-8]
	
	# LOESS smooth data
	for (i in 2:(length(sites))) { 
		counts.loess <- loess(sites[,i] ~ sites$Pos, span=1, data.frame(x=sites$Pos, y=sites[,i]), control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
		counts.predict <- predict(counts.loess, data.frame(x=sites$Pos))
		counts.ratio <- counts.predict/median(counts.predict)
	    sites[,i] <- sites[,i]/counts.ratio
	}
	
	# Normalize data by reads/site
	library(DESeq2)
	colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps), rep(test_pfx, test_reps)), condition=c(rep("untreated", ctrl_reps), rep("treated", test_reps)))
	sitescds <- sites[,2:length(sites)] %>% round %>% DESeqDataSetFromMatrix(colData = colData, design = ~ condition)
	sitescds <- estimateSizeFactors(sitescds)
	counts.norm <- counts(sitescds, normalized=T)
	rownames(counts.norm) <- sites$Pos
	counts.df <- data.frame(counts.norm)

	# Initialize the list of genes, read counts per gene, and number of independent Tn sites per gene
	gff <- read.delim(file=paste(gff_pfx, ".trunc.gff", sep=""), sep="\t", fill=TRUE, skip=2, header=FALSE, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att", "KO", "pathways"))
	gff <- tail(gff, n=-2)
	genomelength <- as.numeric(max(gff$end)*1.1)
	gff <- gff[(gff$feature=="CDS"),]
	controlreps <- 0
	testreps <- 0
	for (c in 1:length(counts.norm[1,])) {
		gff[,c+11] <- rep(1,length(gff[,1]))
		if (controlreps < ctrl_reps) {
			controlreps <- controlreps + 1
			colnames(gff)[c+11] <- paste(ctrl_pfx, controlreps, sep=".")
		}
		else {
			testreps <- testreps + 1
			colnames(gff)[c+11] <- paste(test_pfx, testreps, sep=".")
		}
	}
	
	# Output gene boundaries and read counts per Tn site for Perl binning script
	print("Binning read counts by gene boundaries")
	boundariesfile <- paste(out_pfx, ".boundaries.tsv", sep="")
	sitecountsfile <- paste(out_pfx, ".sitecounts.tsv", sep="")
	write.table(gff[,c(4,5, 12:length(gff))], boundariesfile, quote=FALSE, sep="\t", row.names=F)
	write.table(counts.df, sitecountsfile, quote=FALSE, sep="\t", row.names=T)
	system(paste("perl TnGeneBin.pl", boundariesfile, sitecountsfile))
	genecounts <- read.table(paste(boundariesfile, "out", sep="."), header=T)[,-c(1,2)]
	numsites <- read.table(paste(boundariesfile, "numsites.out", sep="."), header=T)[,-c(1,2)]
	system(paste("rm", boundariesfile,
		paste(boundariesfile, "out", sep="."),
		paste(boundariesfile, "numsites.out", sep=".")))

	# Uncomment this section if you have a kegg annotation description file of the genes and their products
	#genes <- read.delim(file=paste(gff_pfx, ".gene.products.kegg.txt", sep=""), sep="\t", header=TRUE)
	#rownames(genecounts) <- genes$Locus
	#write.table(genecounts, paste(out_pfx, ".genecounts.tsv", sep=""), quote=FALSE, sep="\t")

	# Uncomment this section if you DO NOT have a kegg annotation description file of the genes and their products
	genes <- matrix("", length(gff[,1]), 4)
	for (i in 1:length(gff[,1])) {
		genes[i,1] <- strsplit(grep("locus_tag",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
		genes[i,2] <- strsplit(grep("product",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
		genes[i,3] <- as.character(gff$KO[i])
		genes[i,4] <- as.character(gff$pathways[i])	
	}
	colnames(genes) <- c("id", "product", "KO", "pathways")
	write.table(genecounts, paste(out_pfx, ".genecounts.tsv", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
	
	# Perform differential fitness analysis and output results
	colnames(numsites) <- colnames(gff)[12:length(gff)]
	
	genescds <-DESeqDataSetFromMatrix(countData = round(genecounts), colData=colData, design = ~ condition) 
	genescds$sizeFactor <- rep(1, length(genecounts[1,])) # This is manually set as 1 because we normalized by site above
	#genescds <- estimateDispersions(genescds, fitType="local") # Use this if estimateDispersions fails
	genescds <- estimateDispersions(genescds)
	genescds <- nbinomWaldTest(genescds)
	res <- results(genescds, contrast = c("condition", "treated", "untreated"))
	print(head(res)) 
	#colnames(res)[3] <- paste(ctrl_pfx, "Mean", sep="")
	#colnames(res)[4] <- paste(test_pfx, "Mean", sep="")
	res <- cbind(res, genes[,1:4], numsites) # Uncomment if you have a kegg annotation
	#res <- cbind(res, genes[,1:2], numsites) %>% tbl_df # Uncomment if you do not have a kegg annotation
	write.table(res, file=paste(out_pfx, ".DESeq.tsv", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
	return(res)
}
Args <- commandArgs(TRUE)
TnSeqDESeq(Args[1], as.numeric(Args[2]), Args[3], as.numeric(Args[4]), Args[5], Args[6], as.numeric(Args[7]), Args[-(1:7)])
