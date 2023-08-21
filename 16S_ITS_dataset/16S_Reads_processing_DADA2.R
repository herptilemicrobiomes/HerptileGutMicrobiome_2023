#R Libraries needed

library(dada2)
library(DECIPHER)
library(ggplot2)
library(ShortRead)
library(Biostrings) 

#Establish working directory:

setwd("~/16S_cryo")


#Establish the directory where the input and unzipped files are placed

input_files <- "/0_Raw_reads" 
list.files(input_files)

#Read the names of the fastq files and extract the names
fnFs <- sort(list.files(input_files, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(input_files, pattern="_R2_001.fastq", full.names = TRUE))

#Identify the primers (515F and 806R designed by Caporaso):
#Do not forget to change the primers sequences

FWD <- "GTGCCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"


### VERIFYING THE PRESENCE AND ORIENTATION OF THE PRIMERS ON THE READS.

#Function to know all the primers orientations

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD) #Creating a vector with all the primers orientations
REV.orients <- allOrients(REV)


#Creating a sub directory where the "modified" reads will be stored
## CHANGE ME to the directory containing the fastq files.
path1 <- "/nfs1/BPP/Spatafora_Lab/VargasGL/16S_ITS_2022/Cryo/16S_cryo/2_Cutadapt"

fnFs.filtN <- file.path(path1, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path1, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Count the number of times that the primers appear on the reads:

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
count_primers <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
                       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


### Removing primers with Cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html)

#"Calling" cutadapt as shell command from R

cutadapt <- "/local/cluster/cutadapt/bin/cutadapt"
system2(cutadapt, args = "--version")

# Creating a new directory to store the cutadapt output files

path.cut <- file.path(path1, "cutadapt_output")

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


## Counting again the presence of primers

count_primers_AfterCut <-rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
                               FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
                               REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
                               REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#Saving files

write.table(FWD.orients, file="FWD.orients.tsv", sep="\t", quote=F, col.names=NA)
write.table(REV.orients, file="REV.orients.tsv", sep="\t", quote=F, col.names=NA)

write.table(count_primers, file="count_primers_BEFORE_cutadapt.tsv", sep="\t", quote=F, col.names=NA)
write.table(count_primers_AfterCut, file="count_primers_AFTER_cutadapt.tsv", sep="\t", quote=F, col.names=NA)

# DADA2 -16S processing - AFTER CUTADAPT

#    Target region: V4 region of 16S SSU rRNA (~250 - 390 bp).

### Indicating where the new files to use are located, in this case, where the cutadapt resulting files are.

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))


#Extract sample names, assuming filenames have a format:

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 1)
sample.names

#Inspecting read quality profiles:

ForwardExamplePlot<-plotQualityProfile(fnFs[3:2])
ReverseExamplePlot<-plotQualityProfile(fnRs[3:2])

#Saving Quality plots

ggsave("Forward_plots_example.pdf", ForwardExamplePlot, device="pdf")
ggsave("Reverse_plots_example.pdf", ReverseExamplePlot, device="pdf")


#Creating a new folder to store the new filtered files:

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))



#Filter and trim
#cutFS and cutRs: where the fastq files (after cut adapt) are stored
#filtFs and filtRs: where the filtered files will be stored after running the filterAndTrim command
#truncLen: To set this command, you need to look into the quality scores of your data. This command should be set where in your forward and reverse reads the average base quality steeply drops in quality. In our case, the Forward reads have better quality than the Reverse reads, and the Phared Score drops between 240 and 250 bp, but for Reverse redas
#maxEE: Refers to the amount of expected errors in reads, forward and reverse. The higher the value, more reads will pass


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(230,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)

## LEARN THE ERROR RATES - The parametric error model.

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Plotting and saving plot in PDF format

errorplots <- plotErrors(errF, nominalQ=TRUE)
ggsave("Error_plots_example.pdf", errorplots, device="pdf")


#Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


##Merge Pair ends and editing seqs by length

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)



# Inspect the merger data.frame from the first sample

head(mergers[[1]])
seqtable <- makeSequenceTable(mergers)
dim(seqtable)

table(nchar(getSequences(seqtable)))

length_16S<-hist(nchar(getSequences(seqtable)), main="Distribution of sequence lengths")
jpeg('length_16S.jpg')
plot(length_16S)
dev.off()

#Filter those seqs that the length is below 250

seqtab2 <- seqtable[,nchar(colnames(seqtable)) %in% 250:350]
table(nchar(getSequences(seqtab2)))

## REMOVE CHIMERAS

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#To know the percentage of seqs retained after chimera filtering.

sum(seqtab.nochim)/sum(seqtab2)


##Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#Assign taxonomy with DADA2

taxa16S <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, minBoot=80)



### SAVE DATA ### 

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")


for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "outputs/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "outputs/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
#ASV taxonomy
asv_tax_dadaassign <- taxa16S
row.names(asv_tax_dadaassign) <- sub(">", "", asv_headers)

#merging abundance and tax table
asv_tax_dadaassign <- merge(asv_tab, asv_tax_dadaassign, by=0)

write.table(asv_tax_dadaassign, "outputs/ASVs_taxonomy_dadaassing.tsv", sep = "\t", quote=F, col.names=NA)
