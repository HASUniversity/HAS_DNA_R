library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)


d <- "E:/Seq/S4B/ITS/" #path to directory with the FASTQ files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- head(sort(list.files(d, pattern="_R1.fastq.gz", full.names = TRUE)), -1) #laatste weglaten
fnRs <- head(sort(list.files(d, pattern="_R2.fastq.gz", full.names = TRUE)), -1)

FWD <- "GAACGCAGCRAAIIGYGA" #NB: loopt vast op I
REV <- "TCCTCCGCTTATTGATATGC"

# Check primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Remove sequences with Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# count primer number of times primers appear.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(
  FWD.ForwardReads = sapply(
    FWD.orients, 
    primerHits, 
    fn = fnFs.filtN[[1]]
    ), 
  FWD.ReverseReads = sapply(
    FWD.orients,
    primerHits, 
    fn = fnRs.filtN[[1]]
    ), 
  REV.ForwardReads = sapply(
    REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), 
  REV.ReverseReads = sapply(
    REV.orients, 
    primerHits, 
    fn = fnRs.filtN[[1]])
  )

# NB cutadapt installeren op linux-pc!
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

