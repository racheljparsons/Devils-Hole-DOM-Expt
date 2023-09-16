library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(DECIPHER)
path <- '/home/ubuntu/data/DH/'

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

#FWD <- "GTGYCAGCMGCCGCGGTAA"
#REV <- "GGACTACNVGGGTWTCTAAT"

# Need to reverse the FWD and REV primers as OSU seems to use a different sequencing protocol that allows the adapters to bind to the other end.
FWD <- "GGACTACNVGGGTWTCTAAT"
REV <- "GTGYCAGCMGCCGCGGTAA"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- '/home/ubuntu/miniconda3/envs/cutadapt/bin/cutadapt'

path.cut <- file.path(path, "cutadapt")
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

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


path <- '/home/ubuntu/data/DH/cutadapt'
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <-gsub('.*-(DH[0-9]+)_.*', '\\1', basename(fnFs)) ## Samples DH1-96

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Check quality
plotQualityProfile(fnFs[1:4])

#I would probably trim forward reads at 230 here.

plotQualityProfile(fnRs[1:4])

#The quality in the reverse reads is pretty low and more of a steady decline than a crash, so it's hard to pick where we should
#cut. Let's pick a point where they drop below 30 (say position 175) 


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,175),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#Now we need to learn the error rates for dada2
#here we are looking for two things:

# 1. Does the estimated error rate (black line map well to the points)
# 2. Does the error rate drop with increasing quality?
# The red line is the expected error rate just based on Q-value.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplication for identifying unique sequences

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Application of dada2 pipeline.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#merge our forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


#construct ASV table
seqtab <- makeSequenceTable(mergers)

#check read lengths
table(nchar(getSequences(seqtab)))

#plot distribution
plot(table(nchar(getSequences(seqtab))))

#looking at the distribution, let's toss anything that is <252 and more than 256
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(252,256)]

#Now remove chimeras
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#assign taxa to our nonchimeric reads, try reverse complementing the reads too.
taxa <- assignTaxonomy(seqtab2.nochim, "~/data/silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, "~/data/silva/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.table(taxa.print, file="BATS_OMZ.csv",sep=",",row.names=F)
