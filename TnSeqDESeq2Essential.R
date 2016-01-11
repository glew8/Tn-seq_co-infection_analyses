# Defined as a function, return value is table of differential fitness results with essentiality calls (source into R/RStudio for interactive analyses)

TnSeqDESeqEssential <- function(ctrl_pfx, ctrl_reps, gff_pfx, out_pfx, to_trim, num_expected, in_files) {
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

	# OPTIONAL - perform site filtering. Example: only consider sites identified in both of 2 replicates
	sites <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps == 2)
	sites <- sites[-4]
	
	# LOESS smooth data
	for (i in 2:(length(sites))) { 
		counts.loess <- loess(sites[,i] ~ sites$Pos, span=1, data.frame(x=sites$Pos, y=sites[,i]), control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
		counts.predict <- predict(counts.loess, data.frame(x=sites$Pos))
		counts.ratio <- counts.predict/median(counts.predict)
	    sites[,i] <- sites[,i]/counts.ratio
	}

	# Normalize data by reads/site
	library(DESeq2)
	colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps)), condition = rep("untreated", ctrl_reps))
	sitescds <- sites[,2:length(sites)] %>% round %>% DESeqDataSetFromMatrix(colData = colData, design= ~ 1)
	sitescds <- estimateSizeFactors(sitescds)
	
	counts.norm <- counts(sitescds, normalized=T)
	rownames(counts.norm) <- sites$Pos
	
	# Initialize the list of genes, determine genome length
	# Expects a standard gff file appended with Kegg Orthology and Pathways information appended at the end of each line, separated by tabs
	gff <- read.delim(file=paste(gff_pfx, ".trunc.gff", sep=""), sep="\t", fill=TRUE, header=FALSE, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att")) 
	#colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att", "KO", "pathways")
	#genomelength <- as.numeric(strsplit(as.character(gff[3,1]), " ")[[1]][4])
	genomelength <- as.numeric(max(gff$end)*1.1)
	#print(head(gff))
	#gff <- tail(gff, n=-2)
	#change this value if your gff file contains anything besides genes
	gff <- gff[(gff$feature=="gene"),]
	#print(head(gff))
	# Generate pseudo-datasets with the same number of insertion sites and total reads mapping to those sites, randomly distributed across the genome
	print("Generating pseudo-datasets")
	counts.df <- data.frame(counts.norm)
	counts.df$Pos <- as.numeric(rownames(counts.df))
	numreads <- sum(counts.norm)/ctrl_reps
	numsites <- length(which(counts.norm>0))/ctrl_reps
	for (i in 1:num_expected) {
		expected <- data.frame(Pos=sample(1:genomelength, numsites), Exp=sample(sites$V1, numsites)) %>% arrange(Pos)
		colnames(expected)[2] <- paste("Expected", i, sep=".")
		counts.df <- merge(counts.df, expected, by="Pos", all=T) %>% arrange(Pos)
		counts.df[is.na(counts.df)] <- 0
	}
	rownames(counts.df) <- counts.df$Pos
	counts.norm <- as.matrix(counts.df[,(2:length(counts.df))])

	# Initialize the lists of read counts per gene and number of independent Tn sites per gene
	controlreps <- 0
	expreps <- 0
	for (c in 1:length(counts.norm[1,])) {
		gff[,c+9] <- rep(1,length(gff[,1]))
		if (controlreps < ctrl_reps) {
			controlreps <- controlreps + 1
			colnames(gff)[c+9] <- paste(ctrl_pfx, controlreps, sep=".")
		}
		else {
			expreps <- expreps + 1
			colnames(gff)[c+9] <- paste("Expected", expreps, sep=".")
		}
	}

	# Output gene boundaries and read counts per Tn site for Perl binning script
	print("Binning read counts by gene boundaries")
	boundariesfile <- paste(out_pfx, ".boundaries.tsv", sep="")
	sitecountsfile <- paste(out_pfx, ".sitecounts.tsv", sep="")
	write.table(gff[,c(4,5, 10:length(gff))], boundariesfile, quote=FALSE, sep="\t", row.names=F)
	write.table(counts.df, sitecountsfile, quote=FALSE, sep="\t", row.names=F)
	system(paste("perl ~/local/bin/TnGeneBin.pl", boundariesfile, sitecountsfile))
	genecounts <- read.table(paste(boundariesfile, "out", sep="."), header=T)[,-c(1,2)]
	numsites <- read.table(paste(boundariesfile, "numsites.out", sep="."), header=T)[,-c(1,2)]
	system(paste("rm", boundariesfile,
		paste(boundariesfile, "out", sep="."),
		paste(boundariesfile, "numsites.out", sep=".")))

	# Uncomment this section if you DO NOT have a kegg annotation description file of the genes and their products
	genes <- data.frame(id = rep("", length(gff[,1]), stringsAsFactors = FALSE))#, length(gff[,1]), 1)
	genes$id <- as.character(genes$id)
	
	for (i in 1:length(gff[,1])) {
		genes$id[i] <- strsplit(grep("locus_tag",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
	}
	
	colnames(genes) <- c("id")
	write.table(genecounts, paste(out_pfx, ".genecounts.tsv", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

	# Perform differential fitness analysis
	colnames(numsites) <- colnames(gff)[12:length(gff)]
	numsitesout <- data.frame(numsites[,(1:ctrl_reps)])
	numsitesout[,ctrl_reps+1] <- rowMeans(numsites[,-(1:ctrl_reps)])
	colnames(numsitesout)[ctrl_reps+1] <- "Expected"
	colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)), condition = c(rep(ctrl_pfx, ctrl_reps),rep("Expected", num_expected)))
	genescds <- DESeqDataSetFromMatrix(countData = round(genecounts), colData = colData, design = ~ condition)
	#genescds <- newCountDataSet(round(genecounts), c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)))
	genescds$sizeFactor <- rep(1, length(genecounts[1,])) # This is manually set as 1 because we normalized by site above
	genescds <- estimateDispersions(genescds)
	genescds <- nbinomWaldTest(genescds)
	res <- results(genescds, contrast = c("condition", ctrl_pfx, "Expected"))
	print(head(res)) 
	#colnames(res)[4] <- paste(ctrl_pfx, "Mean", sep="")
	#colnames(res)[3] <- "ExpectedMean"
	out <- cbind(res, genes$id, numsitesout) # Uncomment if you have a kegg annotation
	colnames(out)[6] <- "id"
	#out <- cbind(res, genes[,2:3], numsitesout) %>% tbl_df # Uncomment if you do not have a kegg annotation

	# Perform bimodal clustering and essentiality calling and output results
	library(mclust)
	fit <- Mclust(out$log2FoldChange, G=1:2)
	category <- rep("",length(out$id))
	for (i in 1:length(out$id)) {
		if (fit$classification[i] == 1 & out$log2FoldChange[i] < 0) {
			category[i] <- "Reduced"
		}
		else {
			category[i] <- "Unchanged"
		}
	}

	fit$uncertainty[which(out$log2FoldChange > 0)] <- 0
	print(head(category, 10))
	essentiality <- as.data.frame(cbind(category, fit$uncertainty))
	colnames(essentiality) <- c("Essentiality", "Uncertainty")
	out <- cbind(out, essentiality) 
	write.table(out, file=paste(out_pfx, ".DESeq.tsv", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
	return(out)
}
Args <- commandArgs(TRUE)
TnSeqDESeqEssential(Args[1], as.numeric(Args[2]), Args[3], Args[4], as.numeric(Args[5]), as.numeric(Args[6]), Args[-(1:6)])
