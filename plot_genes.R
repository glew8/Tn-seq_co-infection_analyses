##Scripting for making pretty Gviz charts of genes
#Necessary libraries
library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
options(ucscChromosomeNames=FALSE)
#Read in wkb2 genome from ncbi gff file and make transcript db
setwd("~/Dropbox/Moran/genomes/")
txdb <- makeTxDbFromGFF("wkb2_ncbi.gff")

##put your favorite gene here!  Use the "new" locus tag from ncbi gff file
## maybe can update?
gene = "SALWKB2_RS10635"
name = "Lipopolysaccharide synthesis protein RfaS"
chromosome = "NZ_CP007446.1"
#cdsBy will extract a list of granges from txdb.  These can be used for the annotation track

gene_list <- cdsBy(txdb, use.names=TRUE)
gen <- genome(txdb)


###Read in data from tn-seq
##Read in locations of insertions.  These files were created with awk
#awk -v OFS="\t" '$1=$1' EXP1-sites.txt > EXP1-sites.tsv
input1 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/CONT1-sites.tabbed.tsv")
input2 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/CONT2-sites.tsv")
input3 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/CONT3-sites.tsv")
exp1 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/EXP1-sites.tsv")
exp2 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/EXP2-sites.tsv")
exp3 <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/EXP2-sites.tsv")
names(input1) <- c("abundance", "location")
names(input2) <- c("abundance", "location")
names(input3) <- c("abundance", "location")
names(exp1) <- c("abundance", "location")
names(exp2) <- c("abundance", "location")
names(exp3) <- c("abundance", "location")

site_counts <- read.table("/Users/seanleonard/Dropbox/Moran/TnSeq/diff_output_3.sitecounts.tsv", header=TRUE)



plot_tn_all <- function(gene, name, expand = 0){
    #use res to get start and stop if necessary
    res <- select(txdb, gene, columns(txdb), keytype="GENEID")
    start_id <- res$CDSSTART - expand
    stop_id <- res$CDSEND + expand
    gene_track <- AnnotationTrack(gene_list[[gene]], name = name, id = name, fill = "darkgray")
    gtrack <- GenomeAxisTrack()
    #plotTracks(list(gtrack, gene_track))
    data <- subset(site_counts, location >= start_id  & location <= stop_id)
    input1_test <- DataTrack(data = data$input1, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Input 1", type = "h", lwd = 2,  ylim = c(0, 100))
    input2_test <- DataTrack(data = data$input2, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Input 2", type = "h",  lwd = 2,  ylim = c(0, 100))
    input3_test <- DataTrack(data = data$input3, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Input 3", type = "h",  lwd = 2,  ylim = c(0, 100))
    exp1_test <- DataTrack(data = data$exp1, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Exp 1", type = "h",  lwd = 2,  ylim = c(0, 100))
    exp2_test <- DataTrack(data = data$exp2, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Exp 2", type = "h",  lwd = 2,  ylim = c(0, 100))
    exp3_test <- DataTrack(data = data$exp3, start = data$location, width = 1, genome = genome(txdb), chr = chromosome, name = "Exp 3", type = "h",  lwd = 2,  ylim = c(0, 100))
    plotTracks(list(gtrack,gene_track, input1_test, input2_test, input3_test, exp1_test, exp2_test, exp3_test), from = start_id , to = stop_id, featureAnnotation="id", sizes = c(1, 1.5, 1, 1, 1, 1, 1, 1))
}
#To plot individual genes
plot_tn_all(gene, name, expand = 1500)

#to plot gene region
expand <- 10000
res <- select(txdb, gene, columns(txdb), keytype="GENEID")
start_id <- 857381
stop_id <- 865074
whole_chromosome <- GeneRegionTrack(txdb, chromosome = chromosome, start = start_id - expand, end = stop_id + expand, shape = "arrow", fill = "light blue")
ht <- HighlightTrack(trackList = list(whole_chromosome), start = start_id, end = stop_id)
plotTracks(list(ht), stacking = "squish")

library(dplyr)
##Average of input and output and normalization
site_counts_avg <- site_counts %>% mutate(input_mean = (input1 + input2 + input3)/3, exp_mean = (exp1 + exp2 + exp3)/ 3)

plot_tn_avg <- function(gene, name, expand = 0, type = "h"){
    #use res to get start and stop if necessary
    res <- select(txdb, gene, columns(txdb), keytype="GENEID")
    start_id <- res$CDSSTART - expand
    stop_id <- res$CDSEND + expand
    gene_track <- AnnotationTrack(gene_list[[gene]], name = name, id = name, fill = "darkgray")
    gtrack <- GenomeAxisTrack()
    #plotTracks(list(gtrack, gene_track))
    data <- subset(site_counts_avg, location >= start_id & location <= stop_id)
    input <- DataTrack(data = data$input_mean, start = data$location, width = 2, genome = genome(txdb), chr = chromosome, name = "Input Mean", type = type, ylim = c(0, 100), lwd = 2)
    exp <- DataTrack(data = data$exp_mean, start = data$location, width = 2, genome = genome(txdb), chr = chromosome, name = "Exp Mean", type = type, ylim = c(0, 100), lwd = 2)
    plotTracks(list(gtrack,gene_track, input, exp), from = start_id, to = stop_id, featureAnnotation="id", sizes = c(1, 1, 1, 1))
    
}

plot_tn_avg_range <- function(start_loc, stop_loc, name, expand = 0, type = "h"){
    #use res to get start and stop if necessary
    res <- select(txdb, gene, columns(txdb), keytype="GENEID")
    start_id <- start_loc - expand
    stop_id <- stop_loc + expand
    gene_track <- AnnotationTrack(gene_list[[gene]], name = name, id = name, fill = "darkgray")
    whole_chromosome <- GeneRegionTrack(txdb, chromosome = chromosome, start = start_id, end = stop_id, shape = "arrow", fill = "dark grey")
    gtrack <- GenomeAxisTrack()
    #plotTracks(list(gtrack, gene_track))
    data <- subset(site_counts_avg, location >= start_id & location <= stop_id)
    input <- DataTrack(data = data$input_mean, start = data$location, width = 2, genome = genome(txdb), chr = chromosome, name = "Input Mean", type = type, ylim = c(0, 100), lwd = 2)
    exp <- DataTrack(data = data$exp_mean, start = data$location, width = 2, genome = genome(txdb), chr = chromosome, name = "Exp Mean", type = type, ylim = c(0, 100), lwd = 2)
    plotTracks(list(gtrack,whole_chromosome, input, exp), from = start_id, to = stop_id, featureAnnotation="id", sizes = c(1, 1, 1, 1))
    
}
plot_tn_avg(gene, name, type = "h", expand = 4000)
plot_tn_avg_range(2320000, 2350000, name = "LPS Syntehsis")
?DataTrack
##Make Real plots
setwd("~/Dropbox/talks/2015_sep_lab_meetings/")
#DNA gyrase A
gene = "SALWKB2_RS10385"
name = "DNA Gyrase A"
pdf("dna_gyraseA_all.pdf", height=8, width = 8)
plot_tn_all(gene, name, expand = 1000)
dev.off()
pdf("dna_gyraseA_avg.pdf", height=8, width = 8)
plot_tn_avg(gene, name, expand = 1000)
dev.off()
#alkaline phosphatase no change
gene = "SALWKB2_RS01075"
name = "Alkaline Phosphatase"
pdf("alkaline_phosphatase_all.pdf", height=8, width = 8)
plot_tn_all(gene, name, expand = 1000)
dev.off()
pdf("alkaline_phosphatase_avg.pdf", height=8, width = 8)
plot_tn_avg(gene, name, expand = 1000)
dev.off()

#differential fitness bee gut
#pyruvate dehydrogenase (aceE)
gene = "aceE"
name = "pyruvate dehydrogenase"
pdf("aceE_all.pdf", height=8, width = 8)
plot_tn_all(gene, name, expand = 1000)
dev.off()
pdf("ace_avg.pdf", height=8, width = 8)
plot_tn_avg(gene, name, expand = 1000)
dev.off()
#plot pyruvate genes
plot_tn_avg_range(857381, 865074, "Pyruvate genes")

gene = "SALWKB2_RS03550"
plot_tn_avg(gene, name, expand = 1000)
