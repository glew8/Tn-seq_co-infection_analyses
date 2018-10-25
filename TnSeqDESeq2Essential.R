
#Clear environment
rm(list=ls())
##################################################################################################################
#manually set this
wd = "C:/Users/Gina/Documents/Whiteley_Lab/AA-co-culture-TN-Seq/Data_Analysis/CoInfTnSeq_combo/results_20170426"
#need to have a folder titled "Essentials_output" in this wd
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
gff_dir = "C:/Users/Gina/Documents/Whiteley_Lab/AA-co-culture-TN-Seq/Data_Analysis/Aa624/GCF_001594265.1_ASM159426v1 9-12-17/GCF_001594265.1_ASM159426v1_genomic.gff/"
gff_dir
genomelength = 6537648 #only matters for non-mariner transposons
genomelength
ctrl_pfx = paste(gff_pfx, "Pg", sep = "_")
ctrl_pfx
ctrl_reps = 2
ctrl_reps
num_expected = 100
num_expected
out_pfx = paste(ctrl_pfx, "Essential_100", sep = "_")
out_pfx
to_trim = 50
to_trim
#can put the "conditions" here, see variable above, the order does not matter
in_files = c("Pg_1_cat", "Pg-2")
in_files

write(out_pfx, file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""))
write(paste( "\n", "input files:", sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(in_files, file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste( "\n", "gff: ", gff_pfx, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste( "\n", "genome length: ", genomelength, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("\n", "to trim: ", to_trim, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
write(paste("number of expected datasets: ", num_expected, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)

##################################################################################################################
        # Read in sites files
        library(dplyr)
        sites <- data.frame(Pos=c(0)) 
        for (i in 1:length(in_files)) {
                newsites <- read.table(paste(paste(in_files[i], sep="/"), "sites.txt", sep="-")) 
                colnames(newsites) <- c(paste("V", i, sep=""), "Pos")
                newsites <- tail(newsites, n=-to_trim) #uncomment if to_trim=0
                newsites <- arrange(newsites,Pos)
                sites <- merge(sites, newsites, all=T) 
        }
        sites <- tail(sites, n=-1)
        sites[is.na(sites)] <- 0

        #output number of total and shared sites across samples. Need to edit based on number of samples.
        sites_total <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps >= 1)
        write(paste("\n","number of sites:", nrow(sites_total)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        sites_shared <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps >= 2)
        write(paste("number of sites shared by at least 2 samples:", nrow(sites_shared)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        # OPTIONAL - perform site filtering. Example: only consider sites identified in both of 2 replicates. Need to edit based on number of samples.
        sites <- sites %>% mutate(numreps = (V1 > 0) + (V2 > 0)) %>% filter(numreps >= 2) %>% select(-numreps)
        
        # LOESS smooth data - Usually skip, not sure if this works
        #for (i in 2:(length(sites))) { 
        #        counts.loess <- loess(sites[,i] ~ sites$Pos, span=1, data.frame(x=sites$Pos, y=sites[,i]), control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
        #        counts.predict <- predict(counts.loess, data.frame(x=sites$Pos))
        #        counts.ratio <- counts.predict/median(counts.predict)
        #    sites[,i] <- sites[,i]/counts.ratio
        #}
        
        # Normalize data by reads/site
        library(DESeq2)
        colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps)), condition = rep("untreated", ctrl_reps))
        sitescds <- sites[,2:length(sites)] %>% round %>% DESeqDataSetFromMatrix(colData = colData, design= ~ 1)
        sitescds <- estimateSizeFactors(sitescds)
        #Output the normalized counts
        counts.norm <- counts(sitescds, normalized=F)
        rownames(counts.norm) <- sites$Pos
        
        # Initialize the list of genes, determine genome length
        setwd(gff_dir)
        #can add column names if your gff has extra columns
        gff <- read.delim(file=paste(gff_pfx, ".trunc.gff", sep=""), sep="\t", fill=TRUE, header=FALSE, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "att"))
        setwd(wd)
        print(head(gff))
        #gff <- tail(gff, n=-6) #if gff has extra info at top
        #change this value if your gff file contains anything besides genes that you want to read in
        gff <- gff[(gff$feature=="gene"),]
        print(head(gff))
        
        # Initialize read counts per gene and number of independent Tn sites per gene
        # Generate pseudo-datasets with the same number of insertion sites and total reads mapping to those sites, randomly distributed across the genome
        print("Generating pseudo-datasets")
        counts.df <- data.frame(counts.norm)
        counts.df$Pos <- as.numeric(rownames(counts.df))
        numreads <- sum(counts.norm)/ctrl_reps
        numsites <- length(which(counts.norm>0))/ctrl_reps
        #generates a new sites data frame and calculates the mean counts for each site in a new column. This column (containing the mean) is then used to sample from to generate the expected datasets.
        sites2 <- sites[2:(ctrl_reps+1)]
        sites2$mean <-rowMeans(sites2) %>% round
        #For mariner transposon, define Possible_sites
        Possible_sites <- read.csv("C:/Users/Gina/Documents/Whiteley_Lab/AA-co-culture-TN-Seq/Data_Analysis/scripts/Aa624_TAsites.csv") #only need to do if mariner transposon
        #compute expected dataset
        for (i in 1:num_expected) {
          expected <- data.frame(Pos=sample(Possible_sites$sites, numsites), Exp=sample(sites2$mean, numsites)) %>% arrange(Pos) #for mariner, may need to change column name depending on sites file
          #expected <- data.frame(Pos=sample(1:genomelength, numsites), Exp=sample(sites2$mean, numsites)) %>% arrange(Pos) #for non-mariner transposons
                colnames(expected)[2] <- paste("Expected", i, sep=".")
                counts.df <- merge(counts.df, expected, by="Pos", all=T) %>% arrange(Pos)
                counts.df[is.na(counts.df)] <- 0
        }
        rownames(counts.df) <- counts.df$Pos
        counts.norm <- as.matrix(counts.df[,(2:length(counts.df))])
        
        write(paste("\n", "numreads: ", numreads, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("numsites: ", numsites, sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)

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
        boundariesfile <- paste("./Essentials_output/", paste(out_pfx, ".boundaries.tsv", sep=""), sep="")
        sitecountsfile <- paste("./Essentials_output/", paste(out_pfx, ".sitecounts.tsv", sep=""), sep="")
        write.table(gff[,c(4,5, 10:length(gff))], boundariesfile, quote=FALSE, sep="\t", row.names=F)
        write.table(counts.df, sitecountsfile, quote=FALSE, sep="\t", row.names=F)
        
        #Use perl script to bin sites into genes
        #You will need to change the line below to where your TnGeneBin.pl script is kept
        #system(paste("perl ~/Desktop/Scripts_for_Denmark/TnGeneBin.pl", boundariesfile, sitecountsfile))
        system(paste("perl TnGeneBin.pl", boundariesfile, sitecountsfile))        
        genecounts <- read.table(paste(boundariesfile, "out", sep="."), header=T)[,-c(1,2)]
        genecounts2 <- read.table(paste(boundariesfile, "out", sep="."), header=T)
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
        write.table(genecounts2, paste("./Essentials_output/", paste(out_pfx, ".genecounts2.tsv", sep=""), sep=""), quote=FALSE, sep="\t", row.names=FALSE)

        # Perform differential fitness analysis
        colnames(numsites) <- colnames(gff)[10:length(gff)] #change this number based on number of columns in gff
        numsitesout <- data.frame(numsites[,(1:ctrl_reps)])
        numsitesout[,ctrl_reps+1] <- rowMeans(numsites[,-(1:ctrl_reps)])
        colnames(numsitesout)[ctrl_reps+1] <- "Expected_numsites"
        colData <- data.frame(c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)), condition = c(rep(ctrl_pfx, ctrl_reps),rep("Expected", num_expected)))
        numcountsout <- data.frame(genecounts2[,1:(ctrl_reps+2)])
        numcountsout[,ctrl_reps+3] <- rowMeans(genecounts2[,-(1:(ctrl_reps+2))])
        colnames(numcountsout)[ctrl_reps+3] <- "Expected_numreads"
        
        ###betaPrior=FALSE
        genescds_F <- DESeqDataSetFromMatrix(countData = round(genecounts), colData = colData, design = ~ condition)
        #genescds_F <- newCountDataSet(round(genecounts), c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)))
        #genescds_F$sizeFactor <- rep(1, length(genecounts[1,])) # This is manually set as 1 because we normalized by site above
        genescds_F <- estimateSizeFactors(genescds_F)
        genescds_F <- estimateDispersions(genescds_F) #can set fitType="local" if getting flag at this step for some samples
        genescds_F <- nbinomWaldTest(genescds_F, betaPrior=FALSE)
        res_F <- results(genescds_F, cooksCutoff = FALSE, contrast = c("condition", ctrl_pfx, "Expected"))
        print(head(res_F)) 
        #colnames(res_F)[4] <- paste(ctrl_pfx, "Mean", sep="")
        #colnames(res_F)[3] <- "ExpectedMean"
        out_F <- cbind(res_F, genes$id, numcountsout, numsitesout) # Uncomment if you have a kegg annotation
        #out_F <- cbind(res_F, genes[,2:3], numsitesout) %>% tbl_df # Uncomment if you do not have a kegg annotation
        colnames(out_F)[2] <- "log2FoldChange_FALSE"
        colnames(out_F)[3] <- "lfcSE_FALSE"
        colnames(out_F)[4] <- "stat_FALSE"
        colnames(out_F)[5] <- "pvalue_FALSE"
        colnames(out_F)[6] <- "padj_FALSE"
        colnames(out_F)[7] <- "id"
        
        #Count number of significant genes; can adjust pvalue and padj here based on output(s) that you are interested in
        write(paste("\n", "Analyses with betaPrior = FALSE (no l2fc shrinkage):", sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste(paste("\n", "number of genes with pvalue <= 0.01 in FALSE:", sep=""), sum(res_F$pvalue<=.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("number of genes with padj <= 0.01 in FALSE:", sum(res_F$padj<=.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        # Perform bimodal clustering and essentiality calling and output results
        library(mclust)
        fit_F <- Mclust(out_F$log2FoldChange_FALSE, G=1:2, modelNames = "V")
        summary(fit_F, parameters = TRUE)
        category_F <- rep("",length(out_F$id))
        for (i in 1:length(out_F$id)) {
          if (fit_F$classification[i] == 1 & out_F$log2FoldChange_FALSE[i] < 0) {
            category_F[i] <- "Reduced"
          }
          else {
            category_F[i] <- "Unchanged"
          }
        }
        
        density_F <- densityMclust(out_F$log2FoldChange_FALSE, G=1:2, modelNames = "V")
        summary(density_F)
        pdf(paste("./Essentials_output/", paste(out_pfx, "FALSE_histogram.pdf", sep="_"), sep=""), width = 10, height = 4)
        par(mfrow = c(1,2))
        plot(density_F, what = "density", data = out_F$log2FoldChange_FALSE, breaks=100)
        plot(fit_F, what = "uncertainty")
        dev.off()

        #Count number of unchanged and reduced genes, l2fc means, and l2fc variances
        write(paste(paste("\n", "Mclust Unchanged FALSE:", sep=""), sum(category_F=="Unchanged")), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust Reduced FALSE:", sum(category_F=="Reduced")), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("\n", "Mclust_reduced_mean_False: ", fit_F$parameters$mean[1], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust_unchanged_mean_False: ", fit_F$parameters$mean[2], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste(paste("\n", "Mclust_reduced_variance_False: ", sep=""), fit_F$parameters$variance$sigmasq[1], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust_unchanged_variance_False: ", fit_F$parameters$variance$sigmasq[2], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        fit_F$uncertainty[which(out_F$log2FoldChange_FALSE > 0)] <- 0
        print(head(category_F, 10))
        essentiality_F <- as.data.frame(cbind(category_F, fit_F$uncertainty))
        colnames(essentiality_F) <- c("Essentiality_FALSE", "Uncertainty_FALSE")
        out_F <- cbind(out_F, essentiality_F) 
        
        #Count number of significant genes with reduced Mclust; can adjust padj (or pvalue) and uncertaintly values here based on output(s) that you are interested in
        write(paste("\n", "Reduced and Uncertainty <=0.01 if FALSE:", sum(out_F$Essentiality_FALSE=="Reduced" & as.numeric(as.character(out_F$Uncertainty_FALSE))<=0.01)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Reduced and Uncertainty <=0.01 and padj <=0.01 if FALSE:", sum(out_F$Essentiality_FALSE=="Reduced" & out_F$padj_FALSE<=.01 & as.numeric(as.character(out_F$Uncertainty_FALSE))<=0.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("\n", "Reduced and Uncertainty <=0.05 if FALSE:", sum(out_F$Essentiality_FALSE=="Reduced" & as.numeric(as.character(out_F$Uncertainty_FALSE))<=0.05)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Reduced and Uncertainty <=0.05 and padj <=0.05 if FALSE:", sum(out_F$Essentiality_FALSE=="Reduced" & out_F$padj_FALSE<=.05 & as.numeric(as.character(out_F$Uncertainty_FALSE))<=0.05, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        write.csv(out_F, file=paste("./Essentials_output/", paste(out_pfx, ".FALSE.DESeq.csv", sep=""), sep=""), row.names=T)
        
        ###betaPrior=TRUE
        genescds_T <- DESeqDataSetFromMatrix(countData = round(genecounts), colData = colData, design = ~ condition)
        #genescds_T <- newCountDataSet(round(genecounts), c(rep(ctrl_pfx, ctrl_reps), rep("Expected", num_expected)))
        #genescds_T$sizeFactor <- rep(1, length(genecounts[1,])) # This is manually set as 1 because we normalized by site above
        genescds_T <- estimateSizeFactors(genescds_T)
        genescds_T <- estimateDispersions(genescds_T) #can set fitType="local" if getting flag at this step for some samples
        genescds_T <- nbinomWaldTest(genescds_T, betaPrior=TRUE)
        res_T <- results(genescds_T, cooksCutoff = FALSE, contrast = c("condition", ctrl_pfx, "Expected"))
        print(head(res_T)) 
        #colnames(res_T)[4] <- paste(ctrl_pfx, "Mean", sep="")
        #colnames(res_T)[3] <- "ExpectedMean"
        out_T <- cbind(res_T, genes$id, numcountsout, numsitesout) # Uncomment if you have a kegg annotation
        #out_T <- cbind(res_T, genes[,2:3], numsitesout) %>% tbl_df # Uncomment if you do not have a kegg annotation
        colnames(out_T)[2] <- "log2FoldChange_TRUE"
        colnames(out_T)[3] <- "lfcSE_TRUE"
        colnames(out_T)[4] <- "stat_TRUE"
        colnames(out_T)[5] <- "pvalue_TRUE"
        colnames(out_T)[6] <- "padj_TRUE"
        colnames(out_T)[7] <- "id"

        #Count number of significant genes; can adjust pvalue and padj here based on output(s) that you are interested in
        write(paste("\n", "Analyses with betaPrior = TRUE (l2fc shrinkage):", sep=""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste(paste("\n", "number of genes with pvalue <= 0.01 in TRUE:", sep=""), sum(res_T$pvalue<=.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("number of genes with padj <= 0.01 in TRUE:", sum(res_T$padj<=.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        # Perform bimodal clustering and essentiality calling and output results
        library(mclust)
        fit_T <- Mclust(out_T$log2FoldChange_TRUE, G=1:2, modelNames = "V")
        summary(fit_T, parameters = TRUE)
        category_T <- rep("",length(out_T$id))
        for (i in 1:length(out_T$id)) {
          if (fit_T$classification[i] == 1 & out_T$log2FoldChange_TRUE[i] < 0) {
            category_T[i] <- "Reduced"
          }
          else {
            category_T[i] <- "Unchanged"
          }
        }
        
        density_T <- densityMclust(out_T$log2FoldChange_TRUE, G=1:2, modelNames = "V")
        summary(density_T)
        pdf(paste("./Essentials_output/", paste(out_pfx, "TRUE_histogram.pdf", sep="_"), sep=""), width = 10, height = 4)
        par(mfrow = c(1,2))
        plot(density_T, what = "density", data = out_T$log2FoldChange_TRUE, breaks=100)
        plot(fit_T, what = "uncertainty")
        dev.off()
        
        #Count number of unchanged and reduced genes, l2fc means, and l2fc variances
        write(paste(paste("\n", "Mclust Unchanged TRUE:", sep=""), sum(category_T=="Unchanged")), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust Reduced TRUE:", sum(category_T=="Reduced")), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste(paste("\n","Mclust_reduced_mean_TRUE: ", sep=""), fit_T$parameters$mean[1], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust_unchanged_mean_TRUE: ", fit_T$parameters$mean[2], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste(paste("\n", "Mclust_reduced_variance_TRUE: ", sep=""), fit_T$parameters$variance$sigmasq[1], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Mclust_unchanged_variance_TRUE: ", fit_T$parameters$variance$sigmasq[2], sep =""), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        fit_T$uncertainty[which(out_T$log2FoldChange_TRUE > 0)] <- 0
        print(head(category_T, 10))
        essentiality_T <- as.data.frame(cbind(category_T, fit_T$uncertainty))
        colnames(essentiality_T) <- c("Essentiality_TRUE", "Uncertainty_TRUE")
        out_T <- cbind(out_T, essentiality_T) 
        
        #Count number of significant genes with reduced Mclust; can adjust padj (or pvalue) and uncertaintly values here based on output(s) that you are interested in
        write(paste("\n", "Reduced and Uncertainty <=0.01 if TRUE:", sum(out_T$Essentiality_TRUE=="Reduced" & as.numeric(as.character(out_T$Uncertainty_TRUE))<=0.01)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Reduced and Uncertainty <=0.01 and padj <=0.01 if TRUE:", sum(out_T$Essentiality_TRUE=="Reduced" & out_T$padj_TRUE<=.01 & as.numeric(as.character(out_T$Uncertainty_TRUE))<=0.01, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("\n", "Reduced and Uncertainty <=0.05 if TRUE:", sum(out_T$Essentiality_TRUE=="Reduced" & as.numeric(as.character(out_T$Uncertainty_TRUE))<=0.05)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        write(paste("Reduced and Uncertainty <=0.05 and padj <=0.05 if TRUE:", sum(out_T$Essentiality_TRUE=="Reduced" & out_T$padj_TRUE<=.05 & as.numeric(as.character(out_T$Uncertainty_TRUE))<=0.05, na.rm=TRUE)), file = paste("./Essentials_output/", paste(out_pfx, "_stats.txt", sep=""), sep=""), append = TRUE)
        
        write.csv(out_T, file=paste("./Essentials_output/", paste(out_pfx, ".TRUE.DESeq.csv", sep=""), sep=""), row.names=T)
        
        #Combine betaPrior=TRUE and betaPrior=FALSE data
        out_combo <- cbind(out_F[,1:6], out_F[,(ctrl_reps+ctrl_reps+12):(ctrl_reps+ctrl_reps+13)], out_T[,1:6], out_T[,(ctrl_reps+ctrl_reps+12):(ctrl_reps+ctrl_reps+13)], out_T[,7:(ctrl_reps+ctrl_reps+11)])
        write.csv(out_combo, file=paste("./Essentials_output/", paste(out_pfx, ".combo.DESeq.csv", sep=""), sep=""), row.names=T)
        
        gc()