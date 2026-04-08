library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)


d <- "E:/Seq/S4B/ITS/" #path to directory with the FASTQ files
d <- "/mnt/Bravo/ITS"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- head(sort(list.files(d, pattern="_R1.fastq.gz", full.names = TRUE)), -1) #laatste weglaten
fnRs <- head(sort(list.files(d, pattern="_R2.fastq.gz", full.names = TRUE)), -1)

# FWD <- "GAACGCAGCRAAIIGYGA" #NB: loopt vast op I
FWD <- "GAACGCAGCRAANNGYGA" #NB: loopt vast op I
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
fnFs.filtN <- file.path(d, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(d, "filtN", basename(fnRs))
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

# Remove primers

# cutadapt <- "/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
cutadapt <- "C:/Users/SmMa/AppData/Local/Packages/PythonSoftwareFoundation.Python.3.13_qbz5n2kfra8p0/LocalCache/local-packages/Python313/Scripts"
system2(cutadapt, args = "--version") # Run shell commands from R


path.cut <- file.path(d, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(
  FWD.ForwardReads = sapply(
    FWD.orients, 
    primerHits, 
    fn = fnFs.cut[[1]]
    ), 
  FWD.ReverseReads = sapply(
    FWD.orients,
    primerHits, 
    fn = fnRs.cut[[1]]
    ), 
  REV.ForwardReads = sapply(
    REV.orients, 
    primerHits,
    fn = fnFs.cut[[1]]
    ), 
  REV.ReverseReads = sapply(
    REV.orients, 
    primerHits, 
    fn = fnRs.cut[[1]]
    )
  )

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: **i5SAMPLENAME_R**
sample.names <- str_match(
  fnFs, "i5.(.*)-ITS"
  )[, 2]
# # Extract sample names, assuming filenames have format:
# get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
# sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(
  cutFs, 
  filtFs, 
  cutRs, 
  filtRs, 
  maxN = 0, 
  maxEE = c(2, 2), 
  truncQ = 2,
  minLen = 50, 
  rm.phix = TRUE, 
  compress = TRUE, 
  multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# Learning the error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# Sample inference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

# Check reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
unite.ref <- "~/Documents/DADA2_analyses/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#Use phyloseq to visualize
theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
# subject <- sapply(strsplit(samples.out, "\\."), `[`, 2)
# head(subject)


#load metadata
Metadata <- readxl::read_excel("/mnt/Charlie/users/Lectoraat_Bodem/S4B/DNA extractie_codering.xlsx") %>% 
  rename(sample = `Nr epje`)
# filter  to sample.names
Metadata <- Metadata %>% filter(sample %in% sample.names)
# reorder to sample.names
Metadata <- Metadata[match(sample.names, Metadata$sample),]
# add row names
rownames(Metadata) <- Metadata$sample

# create data frame
S4B_ITS <- phyloseq(
  otu_table(
    seqtab.nochim, 
    taxa_are_rows = FALSE
  ),
  # sample_data(
  #   Metadata
  #   ),
  tax_table(
    taxa
  )
)

dna <- Biostrings::DNAStringSet(taxa_names(S4B_ITS))
names(dna) <- taxa_names(S4B_ITS)
ps <- merge_phyloseq(S4B_ITS, dna)
taxa_names(S4B_ITS) <- paste0("ASV", seq(ntaxa(S4B_ITS)))
S4B_ITS

# opslaan
saveRDS(ps, file = "S4B_ITS.rds")
