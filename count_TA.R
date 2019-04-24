#script to count AT/TA regions in genome

library(Biostrings, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)

setwd("~/Desktop/staph_ref_genome")
getwd()
genome <- readDNAStringSet("Bcenocepacia_J2315_modJV_cat.fa")

#Match TA pattern
ta_sites <- vmatchPattern("TA", genome)
at_sites <- vmatchPattern("AT", genome)

#Make a single sorted vector of all locations
ta_vect <- unname(unlist(start(ta_sites)))
at_vect <- unname(unlist(start(at_sites)))
all_vect <- sort(c(ta_vect, at_vect))
ta_df <- data_frame(sites = all_vect)
write_tsv(ta_df, "Bcen_TAsites.tsv", col_names = TRUE)
