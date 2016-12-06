#script to count AT/TA regions in genome


library(Biostrings, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)

setwd("~/Dropbox/articles/2015_tnseq_snod/analysis_rerun/")
genome <- readDNAStringSet("wkB2/wkB2.fna")

#Match TA pattern
ta_sites <- vmatchPattern("TA", genome)

#Make a single sorted vector of all locations
ta_vect <- unname(unlist(start(ta_sites)))
#at_vect <- unname(unlist(start(at_sites)))
#all_vect <- sort(c(ta_vect, at_vect))
ta_df <- data_frame(sites = ta_vect)
write_tsv(ta_df, "wkB2_TAsites.tsv", col_names = TRUE)
ta_df$abundance <- 1
ta_df <- ta_df %>% dplyr::select(abundance, sites)
write_tsv(ta_df, "ta_wkB2.tsv", col_names = FALSE)

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
options(ucscChromosomeNames=FALSE)

#Read in wkb2 genome from ncbi gff file and make transcript db
setwd("~/Dropbox/articles/2015_tnseq_snod/analysis_rerun/")
txdb <- makeTxDbFromGFF("wkB2/wkB2.gff")

exbygene <- exonsBy(txdb, "gene")
genebygene <- genes(txdb)
transcripts <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
lengths <- width(transcripts)
ta_sites <- as(ta_sites, "GRanges")
seqlevels(ta_sites) <- as.character(seqnames(transcripts))[1]
test <- summarizeOverlaps(transcripts, ta_sites, ignore.strand=TRUE)
hist(assay(test))
head(assay(test))
library(rtracklayer)
test <- import("wkB2/wkB2.gff")
text_txdb <- makeTxDbFromGRanges(test)
columns(text_txdb)
txdb_dump <- as.list(text_txdb)
txdb_dump[["transcripts"]]
test_df <- as.data.frame(test)
locus_tags <- dplyr::filter(test_df, type == "gene") %>% dplyr::select(locus_tag, Name)

rownames(test) <- rowRanges(test)$gene_id
count_data <- as.data.frame(assay(test))
count_data$locus_tag <- rownames(count_data)
count_data$length <- width(test)

colnames(count_data) <- c("TA_sites", "locus_tag", "length")
count_data <- inner_join(count_data, locus_tags, by = c("locus_tag" = "Name"))
count_data <- count_data %>% dplyr::select(TA_sites, locus_tag.y, length) %>% rename(locus_tag = locus_tag.y )
write_tsv(count_data, "TA_sites_wkB2.tsv")
