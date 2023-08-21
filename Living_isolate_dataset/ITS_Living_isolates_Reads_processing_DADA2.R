
#DADA2 pipeline - ITS analyses

#First, we are analyzing the ITS data from BASIDIOBOLUS CULTURES collection to document The Basidiobolus ITS diversity from single cultures

#Target region: ITS1 region (~250 - 390 bp)

#Primers:

#Earth Microbiome Primers:
#ITS1F (Gardes & Bruns 1993) 5' CTTGGTCATTTAGAGGAAGTAA 3'
#ITS2 (White et al. 1990) 5' GCTGCGTTCTTCATCGATGC 3'


#Primers references:
##Gardes M, Bruns TD. 1993. ITS primers with enhanced specificity for basidiomycetes – application
#	to the identification of mycorrhizas and rusts. Mol. Ecol. 2: 113-118.

## White TJ, Bruns TD, Lee S, Taylor J. 1990. Amplification and direct sequencing of fungal 
#	ribosomal RNA genes for phylogenetics. In: Innis MA, Gelfand DH. (eds). PCR Protocols: 
#	A Guide to Methods and Applications. Academic Press: London, pp. 315 – 322.


#Principal pipeline: https://benjjneb.github.io/dada2/ITS_workflow.html

#Call the needed libraries:

library(dada2)
library(DECIPHER)
library(ggplot2)
library(ShortRead)
library(Biostrings) 

#Establish working directory:

setwd("~/ITS_Basidiobolus_results")

#Establish the directory where the input and unzipped files are placed

input_files <- "~/0_RawSeqs/ITS_Basidiobolus_ITS_diversity/"
list.files(input_files)

#Read the names of the fastq files and extract the names
fnFs <- sort(list.files(input_files, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_files, pattern="_R2_001.fastq.gz", full.names = TRUE))

#Identify the primers (515F and 806R designed by Caporaso):

FWD <- "CTTGGTCATTTAGAGGAAGTAA"

REV <- "GCTGCGTTCTTCATCGATGC"


### VERIFYING THE PRESENCE AND ORIENTATION OF THE PRIMERS ON THE READS.
#Note: It could be common to find the forward primer not only at the beginning of the read, but also in other position a long the read, and this occurs when the reads are longer than the target region.

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
FWD.orients
REV.orients


#Creating a sub directory where the "modified" reads will be stored

path1 <- "~/ITS_Basidiobolus_results/cutadapt"  ## CHANGE ME to the directory containing the fastq files.


fnFs.filtN <- file.path(path1, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path1, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Count the number of times that the primers appear on the reads:

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))



### Removing primers with Cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html)

#"Calling" cutadapt as. shell command from R

cutadapt <- "/local/cluster/cutadapt/bin/cutadapt"
system2(cutadapt, args = "--version")

# Creating a new directory to store the cutadapt output files

path.cut <- file.path(path1, "cutadapt")

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



#Count again the number of times that the primers appear on the reads:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# RUNNING ITSXPRESS (SEE ITS_itsxpress.sh)> Conda environment (Running in bash)


### Indicating where the new files to use are located, in this case, where the itsxpress (filtered files are) resulting files are.

path<-"~/cutadapt_first/itsxpress_AfterCutadapt"

filtITSFs <- sort(list.files(path, pattern = "filtered.fastq.gz", full.names = TRUE))
#cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have a format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(filtITSFs, get.sample.name))
head(sample.names)

sample.names <- sapply(strsplit(basename(filtITSFs), "_"), `[`, 1)
sample.names

#Inspecting read quality profiles:
ForwardExamplePlot<-plotQualityProfile(filtITSFs[3:2])
#ReverseExamplePlot<-plotQualityProfile(fnRs[3:2])

#Saving Quality plots
ggsave("Forward_plots_example.pdf", ForwardExamplePlot, device="pdf")
#ggsave("Reverse_plots_example.pdf", ReverseExamplePlot, device="pdf")


#Creating a new folder to store the new filtered files:
filtFs <- file.path(path, "filtered", basename(filtITSFs))
#filtRs <- file.path(path.cut, "filtered", basename(cutRs))


#Filter and trim
#cutFS and cutRs: where the fastq files (after cut adapt) are stored
#filtFs and filtRs: where the filtered files will be stored after running the filterAndTrim command
#truncLen: To set this command, you need to look into the quality scores of your data. This command should be set where in your forward and reverse reads the average base quality steeply drops in quality. In our case, the Forward reads have better quality than the Reverse reads, and the Phared Score drops between 240 and 250 bp, but for Reverse redas
#maxEE: Refers to the amount of expected errors in reads, forward and reverse. The higher the value, more reads will pass


out <- filterAndTrim(filtITSFs, filtFs, minLen = 140, maxLen = 200,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


head(out)



## LEARN THE ERROR RATES - The parametric error model.

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Plotting and saving plot in PDF format
errorplots <- plotErrors(errF, nominalQ=TRUE)
ggsave("Error_plots_example.pdf", errorplots, device="pdf")


#Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")
dadaFs[[1]]


# Inspect the merger data.frame from the first sample
head(dadaFs[[1]])
seqtable <- makeSequenceTable(dadaFs)
dim(seqtable)

table(nchar(getSequences(seqtable)))

## REMOVE CHIMERAS

seqtab.nochim <- removeBimeraDenovo(seqtable, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtable)

##Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign Taxonomy#

unite.ref <- "sh_general_release_dynamic_all_10.05.2021.fasta"  # CHANGE ME to location on your machine


set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, outputBootstraps=TRUE)

#try bootstrap of minBoot=80
taxa_80 <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE, outputBootstraps=TRUE, minBoot=80)

head(unname(taxa))


#### SAVE DATA ### 

## Extracting the standard goods from DADA2 ##

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")


for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_ITS1.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)	

# taxa table
asv_tax_UNITE <- taxa
row.names(asv_tax_UNITE) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_taxa.tsv", sep="\t", quote=F, col.names=NA)	