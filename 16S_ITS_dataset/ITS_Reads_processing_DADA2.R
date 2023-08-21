#DADA2 pipeline - ITS analyses

#First, we are analyzing the ITS data from CRYOPRESERVED collection to document fungal gut diversity across herptiles
#Target region: ITS1 region (~250 - 390 bp)

#Primers:

## Original primers:
# ITS1FI2 (Schmidt et al. 2013) 5′GAACCWGCGGARGGATCA 3′
# ITS2 (White et al. 1990) 5’ GCTGCGTTCTTCATCGATGC ‘3
#Primers references:
## Schmidt, P. A., Bálint, M., Greshake, B., Bandow, C., Römbke, J., & Schmitt, I. (2013).
#	Illumina metabarcoding of a soil fungal community. Soil Biology and Biochemistry, 65, 128-132.

## White TJ, Bruns TD, Lee S, Taylor J. 1990. Amplification and direct sequencing of fungal 
#	ribosomal RNA genes for phylogenetics. In: Innis MA, Gelfand DH. (eds). PCR Protocols: 
#	A Guide to Methods and Applications. Academic Press: London, pp. 315 – 322.



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

#Run R into de cluster by typing R in the shell

#Call the needed libraries:

library(dada2)
library(DECIPHER)
library(ggplot2)
library(ShortRead)
library(Biostrings) 

#Establish working directory:

setwd("/nfs1/BPP/Spatafora_Lab/VargasGL/16S_ITS_2022/ITS_Basidiobolus_results")

#Establish the directory where the input and unzipped files are placed

input_files <- "/nfs1/BPP/Spatafora_Lab/VargasGL/16S_ITS_2022/ITS_Basidiobolus_results/2_ITSxpress"
list.files(input_files)

#Read the names of the fastq files and extract the names
fnFs <- sort(list.files(input_files, pattern="_filtered.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_files, pattern="_R2_001.fastq", full.names = TRUE))

#Identify the primers (515F and 806R designed by Caporaso):

FWD <- "CTTGGTCATTTAGAGGAAGTAA"

REV <- "GCTGCGTTCTTCATCGATGC"


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
FWD.orients
REV.orients


#Creating a sub directory where the "modified" reads will be stored

#path1 <- "/nfs1/BPP/Spatafora_Lab/VargasGL/16S_ITS_2022/ITS_Basidiobolus_results/4_Reads_orientation"  ## CHANGE ME to the directory containing the fastq files.


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



rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))


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


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen = 50,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


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

## REMOVE CHIMERAS

seqtab.nochim <- removeBimeraDenovo(seqtable, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtable)

##Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign Taxonomy#

unite.ref <- "sh_general_release_dynamic_all_10.05.2021.fasta"  # CHANGE ME to location on your machine

taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

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
write(asv_fasta, "ASVs_for_ITSX.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_NOITSX.tsv", sep="\t", quote=F, col.names=NA)	

asv_tax_UNITE <- taxa
row.names(asv_tax_UNITE) <- sub(">", "", asv_headers)
write.table(asv_tax_UNITE, "asv_tax_UNITE.tsv", sep="\t", quote=F, col.names=NA)
