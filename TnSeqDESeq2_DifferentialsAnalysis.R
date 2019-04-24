
# return value is table of differential fitness results for Tn-seq dataset (source into R/RStudio for interactive analyses)
##################################################################################################################
#Clear environment
rm(list=ls())
##################################################################################################################
#manually set this
wd = "~/Desktop/Tn_Analysis"
#need to have a folder titled "DE" in this wd
setwd(wd)
filenames <- dir(pattern = "-sites.txt")
conditions <- sub("-sites.txt", "", filenames) #sub(pattern, replacement, x)
filenames
conditions
##################################################################################################################
#manually change each parameter
gff_pfx = "Aa624_sRNAs_reordered"
gff_pfx
#tell R where your gff file is stored
#standard protocol is to download gff from NCBI and make the following changes: (1) remove any header lines; (2) trim 3' 10% of each gene (currently done manually in excel); (3) name as .trunc.gff
gff_dir = "~/Desktop/Tn_Analysis/Aa624"
gff_dir
#set the name and number of replicates for your control condition
ctrl_pfx = paste(gff_pfx, "mono", sep = "_")
ctrl_pfx
ctrl_reps = 2
ctrl_reps
#set the name and number of replicates for your test condition
test_pfx = paste(gff_pfx, "Pg", sep = "_")
test_pfx
test_reps = 2
test_reps
#set the output file prefix
out_pfx = paste(ctrl_pfx, "vs", test_pfx, sep = "_")
out_pfx
#set the number of sites with the most reads to "trim" (remove) from the data
to_trim = 50
to_trim
#can put the "conditions" here, see variable above, just make sure the order is correct (controls first and then tests)
in_files = c("mono_1_cat", "mono_2_cat", "Pg-1", "Pg-2")
in_files
write(out_pfx, file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""))
write(paste( "\n", "gff: ", gff_pfx, sep=""), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste( "\n", "input files:", sep=""), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(in_files, file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("\n", "to trim: ", to_trim, sep=""), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)

##################################################################################################################
# Read in sites files
library(dplyr)
sites <- data.frame(Pos=c(0)) 
for (i in 1:length(in_files)) {
  newsites <- read.table(paste(paste(in_files[i], sep="/"), "sites.txt", sep="-")) 
  colnames(newsites) <- c(paste("V", i, sep=""), "Pos")
  newsites <- tail(newsites, n=-to_trim) #comment if to_trim=0
  newsites <- arrange(newsites,Pos)
  sites <- merge(sites, newsites, all=T) 
}
sites <- tail(sites, n=-1)
sites[is.na(sites)] <- 0

# OPTIONAL - look at share sites in control and test conditions
# ADJUST IF NOT 2 CONTROL and 2 TEST
sites_ctrl <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps >= 2)
write(paste("\n", paste("number of sites shared by at least 2 control samples:", nrow(sites_ctrl)), sep =""), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
sites_test <- sites %>% mutate(numreps = (V3 > 0) + (V4 > 0)) %>% filter(numreps >= 2)
write(paste("number of sites shared by at least 2 test samples:", nrow(sites_test)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
sites_all <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0) + (V3 > 0) + (V4 >0)) %>% filter(numreps >= 2)
write(paste("number of sites shared by at least 2 samples:", nrow(sites_all)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)


# OPTIONAL - perform site filtering. Example: only consider sites identified in 2 or more of 4 replicates
# ADJUST IF NOT 4 REPLICATES
sites  <- sites %>% mutate(numreps = (V1 > 0 ) + (V2 > 0) + (V3 > 0) + (V4 > 0)) %>% filter(numreps >= 2) %>% select(-numreps)

# Normalize data by reads/site
library(DESeq2)
colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps), rep(test_pfx, test_reps)), condition=c(rep("untreated", ctrl_reps), rep("treated", test_reps)))
sitescds <- sites[,2:length(sites)] %>% round %>% DESeqDataSetFromMatrix(colData = colData, design = ~ condition)
sitescds <- estimateSizeFactors(sitescds)
counts.norm <- counts(sitescds, normalized=F)
rownames(counts.norm) <- sites$Pos
counts.df <- data.frame(counts.norm)

#Make file with sites and counts per site
library('tibble')
counts.df.out <- counts.df
colnames(counts.df.out) <- in_files
counts.df.out <- counts.df.out %>% rownames_to_column('Tn_site')
write.table(counts.df.out, paste("./DE/", paste(out_pfx, ".sitecounts2.tsv", sep=""), sep=""), quote=FALSE, sep="\t", row.names=F, col.names=T)

# Initialize the list of genes, determine genome length
setwd(gff_dir)
#can add column names if your gff has extra columns
gff <- read.delim(file=paste(gff_pfx, ".trunc.gff", sep=""), sep="\t", fill=TRUE, header=FALSE, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att"))
setwd(wd)
print(head(gff))
#gff <- tail(gff, n=-6) #uncommentif gff has extra info at top
#change this value if your gff file contains anything besides genes that you want to read in
gff <- gff[(gff$feature=="gene"),]
print(head(gff))

# Initialize the lists of read counts per gene and number of independent Tn sites per gene
controlreps <- 0
testreps <- 0
for (c in 1:length(counts.norm[1,])) {
  gff[,c+9] <- rep(1,length(gff[,1]))
  if (controlreps < ctrl_reps) {
    controlreps <- controlreps + 1
    colnames(gff)[c+9] <- paste(ctrl_pfx, controlreps, sep=".")
  }
  else {
    testreps <- testreps + 1
    colnames(gff)[c+9] <- paste(test_pfx, testreps, sep=".")
  }
}

# Output gene boundaries and read counts per Tn site for Perl binning script
print("Binning read counts by gene boundaries")
boundariesfile <- paste("./DE/", paste(out_pfx, ".boundaries.tsv", sep=""), sep="")
sitecountsfile <- paste("./DE/", paste(out_pfx, ".sitecounts.tsv", sep=""), sep="")
write.table(gff[,c(4,5, 10:length(gff))], boundariesfile, quote=FALSE, sep="\t", row.names=F)
write.table(counts.df, sitecountsfile, quote=FALSE, sep="\t", row.names=T)

#Use perl script to bin sites into genes
#You will need to change the line below to where your TnGeneBin.pl script is kept
system(paste("perl ~/Desktop/Tn_Analysis/TnGeneBin.pl", boundariesfile, sitecountsfile))        
genecounts <- read.table(paste(boundariesfile, "out", sep="."), header=T)[,-c(1,2)]
numsites <- read.table(paste(boundariesfile, "numsites.out", sep="."), header=T)[,-c(1,2)]
system(paste("rm", boundariesfile,
             paste(boundariesfile, "out", sep="."),
             paste(boundariesfile, "numsites.out", sep=".")))

genes <- data.frame(id = rep("", length(gff[,1]), stringsAsFactors = FALSE))#, length(gff[,1]), 1)
genes$id <- as.character(genes$id)
for (i in 1:length(gff[,1])) {
  genes[i,1] <- strsplit(grep("locus_tag",strsplit(as.character(gff$att[i]),";")[[1]], value=T),"=")[[1]][2]
}
colnames(genes) <- c("id")
write.table(genecounts, paste("./DE/", paste(out_pfx, ".genecounts.tsv", sep=""), sep=""), quote=FALSE, sep="\t", row.names=FALSE)

# Perform differential fitness analysis and output results
colnames(numsites) <- colnames(gff)[10:length(gff)] #edit if extra columns in gff

###betaPrior=FALSE
#When betaPrior=FALSE in DESeq2 step "nbinomWaldTest," log2fold changes are NOT shrunk towards zero, even when counts are low, dispersion is high, or degrees of freedom is low.
#We consider betaPrior=FALSE as the less conservative analysis. It is more likely to identify genes with a large l2fc simply because the data is noisy 
genescds_F <-DESeqDataSetFromMatrix(countData = round(genecounts), colData=colData, design = ~ condition) 
genescds_F <- estimateSizeFactors(genescds_F)
genescds_F <- estimateDispersions(genescds_F) #can set fitType="local" if getting flag at this step for some samples
genescds_F <- nbinomWaldTest(genescds_F, betaPrior = FALSE)

res_F <- results(genescds_F, contrast = c("condition", "treated", "untreated"))
print(head(res_F))

#Count number of significant genes
write(paste(paste("\n", "number of genes with pvalue <= 0.05 in FALSE:", sep=""), sum(res_F$pvalue<=.05, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with padj <= 0.05 in FALSE:", sum(res_F$padj<=.05, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with pvalue <= 0.01 in FALSE:", sum(res_F$pvalue<=.01, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with padj <= 0.01 in FALSE:", sum(res_F$padj<=.01, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)

res_F <- cbind(res_F, genes[,1], genecounts, numsites)
colnames(res_F)[2] <- "log2FoldChange_FALSE"
colnames(res_F)[3] <- "lfcSE_FALSE"
colnames(res_F)[4] <- "stat_FALSE"
colnames(res_F)[5] <- "pvalue_FALSE"
colnames(res_F)[6] <- "padj_FALSE"
colnames(res_F)[7] <- "id"
write.csv(res_F, file=paste("./DE/", paste(out_pfx, ".FALSE.DESeq.csv", sep=""), sep=""), row.names=T)

###betaPrior=TRUE
#When betaPrior=TRUE in DESeq2 step "nbinomWaldTest," log2fold changes ARE shrunk towards zero when counts are low, dispersion is high, or degrees of freedom is low.
#We consider betaPrior=TRUE as the more conservative analysis. It is less likely to identify genes with a large l2fc simply because the data is noisy
#betaPrior=TRUE is also the legacy DESeq2 option prior to version 1.16 (November 2016)
genescds_T <-DESeqDataSetFromMatrix(countData = round(genecounts), colData=colData, design = ~ condition)
genescds_T <- estimateSizeFactors(genescds_T)
genescds_T <- estimateDispersions(genescds_T) #can set fitType="local" if getting flag at this step for some samples
genescds_T <- nbinomWaldTest(genescds_T, betaPrior = TRUE)

res_T <- results(genescds_T, contrast = c("condition", "treated", "untreated"))
print(head(res_T)) 

#Count number of significant genes
write(paste(paste("\n", "number of genes with pvalue <= 0.05 in TRUE:", sep=""), sum(res_T$pvalue<=.05, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with padj <= 0.05 in TRUE:", sum(res_T$padj<=.05, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with pvalue <= 0.01 in TRUE:", sum(res_T$pvalue<=.01, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of genes with padj <= 0.01 in TRUE:", sum(res_T$padj<=.01, na.rm=TRUE)), file = paste("./DE/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)

res_T <- cbind(res_T, genes[,1], genecounts, numsites)
colnames(res_T)[2] <- "log2FoldChange_TRUE"
colnames(res_T)[3] <- "lfcSE_TRUE"
colnames(res_T)[4] <- "stat_TRUE"
colnames(res_T)[5] <- "pvalue_TRUE"
colnames(res_T)[6] <- "padj_TRUE"        
colnames(res_T)[7] <- "id"
write.csv(res_T, file=paste("./DE/", paste(out_pfx, ".TRUE.DESeq.csv", sep=""), sep=""), row.names=T)

#Combine betaPrior=FALSE and betaPrior=TRUE data
out_combo <- cbind(res_T[,1:6], res_F)
write.csv(out_combo, file=paste("./DE/", paste(out_pfx, ".combo.DESeq.csv", sep=""), sep=""), row.names=T)
