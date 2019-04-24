Tn-seq
===========
Custom scripts for analyzing (parsing, mapping, and tallying) Tn-seq reads and
determining differentially abundant transposon insertion mutants.

Copyright (c) 2014 Keith H. Turner, Jake Everett, Urvish Trivedi, Kendra P.
Rumbaugh, and Marvin Whiteley

The scripts contained herein can be used to automatically analyze high-throughput
(Illumina) sequencing reads derived from transposon-genome junctions. First, each
individual dataset is analyzed with TnSeq2.sh. Then two different anlayses can be done. TnSeqDESeq2_DifferentialsAnalysis.R can be used to compare transposon mutant abundance between a control and
test condition. Alternatively, TnSeqDESeq2Essential_randomTn.R or TnSeqDESeq2Essential_mariner.R can be used to perform an essentials analysis (see https://www.pnas.org/content/112/13/4110) to identify genes required for growth under a certain condition.
See below for specific usage details and software dependencies. Please direct any
questions on their use or construction to Keith H. Turner (khturner at utexas.edu).

As of May 2015, this version of scripts is maintained by Sean Leonard (sean.p.leonard at utexas.edu).
Scripts were updated to use DEseq2 for analyses.

As of April 2019, these version of scripts have been updated by Carolyn Ibberson (carolyn.ibberson at biosci.gatech.edu) and Gina Lewin (gina.lewin at biosci.gatech.edu) in the Whiteley lab at Georgia Tech.

TnSeq2.sh
===========
This script takes one FASTQ file specifying read 1 single-end sequencing run done
on a Tn-seq library. This sequencing library should have been prepared with an
end-specific method (i.e. the transposon end should be on read 1), and should
contain a "tag" or "IR" sequence derived from the end of the transposon to identify
which reads are transposon-derived. These reads are found, trimmed of sequence that
will not map to the genome (both transposon- and sequencing adapter-derived), and
mapped to your genome with bowtie2. Finally, insertion site locations and read
counts are tallied. All results and run information is put in a directory named for
your files, and these results can be fed directly into subsequent analyses (see below).

Usage: ./TnSeq2.sh [-i <IR seq>] [-g <path to genome files>] <pfx>

Arguments:

<IR seq>     - The sequence of the transposon end sequence remaining (for junction
   authentication)

<path to genome files>   - The location of the bowtie2 index you're using. You should have already built the bowtie2 indices files.

<pfx>        - the file prefix for your sequence files (If your sequence file is
   named condition1_R1.fastq, the prefix is "condition1")

Dependencies:

-cutadapt (https://cutadapt.readthedocs.io/en/stable/)

-bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
