library(tidyverse)
library(dada2)
library(phyloseq)

d <- "E:/Seq/S4B/16S/" #path to directory with the FASTQ files
# d <- "/mnt/Bravo/16S"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(
  list.files(
    d, 
    pattern="_R1.fastq.gz", 
    full.names = TRUE)
  )[-153] #laatste weglaten
fnRs <- sort(
  list.files(
    d, 
    pattern="_R2.fastq.gz", 
    full.names = TRUE
    )
  )[-153]
# fnFs <- "E:/Seq/S4B/16S/NX.VH01519_270.001.FLD_ill_050_i7---FLD_ill_0242_i5.S50_R1.fastq.gz"
# fnRs <- "E:/Seq/S4B/16S/NX.VH01519_270.001.FLD_ill_050_i7---FLD_ill_0242_i5.S50_R2.fastq.gz"

# Extract sample names, assuming filenames have format: **i5SAMPLENAME_R**
sample.names <- str_match(
  fnFs, "i5.(.*)_R"
  )[, 2]

#Inspect read quality profiles
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(
  d, 
  "filtered", 
  paste0(
    sample.names, 
    "_F_filt.fastq.gz"
    )
  )

filtRs <- file.path(
  d, 
  "filtered", 
  paste0(
    sample.names, 
    "_R_filt.fastq.gz"
    )
  )

names(filtFs) <- sample.names
names(filtRs) <- sample.names

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT" 

out <- filterAndTrim(
  fnFs, 
  filtFs, 
  fnRs, 
  filtRs, 
  truncLen=c(280,240),
  # trimLeft = c(FWD, REV), # primers still presented in the raw data!
  trimLeft = 21, # primers still presented in the raw data!
  maxN=0, 
  maxEE=c(2,2), 
  truncQ=2, 
  rm.phix=TRUE,
  compress=TRUE, 
  multithread=TRUE
  )

#NB: version without trimming out primers

# out <- filterAndTrim(fnFs, 
#                      filtFs, 
#                      fnRs, 
#                      filtRs, 
#                      truncLen=c(280,240),
#                      maxN=0, 
#                      maxEE=c(2,2), 
#                      truncQ=2, 
#                      rm.phix=TRUE,
#                      compress=TRUE, 
#                      multithread=TRUE)
#first number is truncLen for Forward reads; second number is truncLen of Reverse reads)
head(out)

#learning the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#sample inference 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #Chimera's make up about 11 % of the merge sequence variants, is that too much?

#Track read through the pipeline
getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/Lectoraat_Bodem/Documents/DADA2_analyses/Silva_db/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/Lectoraat_Bodem/Documents/DADA2_analyses/Silva_db/silva_v138.2_assignSpecies.fa.gz")

#inspect assignments
taxa.print <- taxa # Removing sequence rownames for display only
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
S4B_16S <- phyloseq(
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
  
dna <- Biostrings::DNAStringSet(taxa_names(S4B_16S))
names(dna) <- taxa_names(S4B_16S)
ps <- merge_phyloseq(S4B_16S, dna)
taxa_names(S4B_16S) <- paste0("ASV", seq(ntaxa(S4B_16S)))
S4B_16S

# opslaan
saveRDS(ps, file = "S4B_16S.rds")
