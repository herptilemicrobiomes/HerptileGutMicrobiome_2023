#########################################################################################################

##Reloading data into Phyloseq

#Loading packages

library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyverse)
library(stats)
library(vegan)
library(rstatix)
library(FSA)
library(ggpubr)
library(data.table)
library(patchwork)
library(packcircles)
library(ggforce)
library(dplyr)
library(BiocGenerics)
library(ampvis2)
library(kableExtra)
library(tibble)
library(fantaxtic)
library(cowplot)

#Establishing working directory

setwd("~/Cryo_2022")
set.seed(1000)

# UPLOADING ITS DATA
its_asv_rare <- read.table("ITS_Cryo/Phyloseq_rarefied_data/220913_CryoITS_rarified.dataset.csv", sep=',', header=T, row.names=1)
its_meta <- read.csv("ITS_Cryo/Phyloseq_rarefied_data/ITS_cryo_metadata.csv", header=T, row.names=1, check.names=F)
its_tax <- read.csv("ITS_Cryo/Phyloseq_rarefied_data/220913_CryoITS_taxonomyfile.csv", colClasses = "character", row.names = 1)

# Create matrix from data frame
# ASV=ASV frequency table and BAC=Bacteria / TAX=Taxonomy information
ITSASVTAB<- data.matrix(its_asv_rare)
ITSTAX<- as.matrix(its_tax)


#Tell phyloseq what is what
ITSASVTAB = otu_table(ITSASVTAB, taxa_are_rows=TRUE)
ITSTAX = tax_table(ITSTAX)
ITSMETA = sample_data(its_meta)

#Integrating all data into a single phyloseq object
## ALL=all data integrated - BAC=Bacteria
ITS_phy_rare = phyloseq(ITSASVTAB,ITSTAX, ITSMETA)
#Visualazing the Phyloseq object
ITS_phy_rare


#Exploring
## Number of samples
nsamples(ITS_phy_rare)
## Number of taxa and rank names
ntaxa(ITS_phy_rare)
rank_names(ITS_phy_rare)
## Listing sample metadata variables
sample_variables(ITS_phy_rare)

#number of samples per host_animal
summary(as.factor(as.data.frame(phyloseq::sample_data(ITS_phy_rare))$host_animal))

#number of samples per host_genus
summary(as.factor(as.data.frame(phyloseq::sample_data(ITS_phy_rare))$host_genus))


#Exploring phyloseq object
summarize_phyloseq(ITS_phy_rare)


# UPLOADING ITS DATA
cryo16S_asv_rare <- read.table("16S_Cryo/Phyloseq_rarefied_data/220913_Cryo16s_rarified.dataset.csv", sep=',', header=T, row.names=1)
cryo16S_meta <- read.csv("16S_Cryo/Phyloseq_rarefied_data/16S_cryo_metadata.csv", header=T, row.names=1, check.names=F)
cryo16S_tax <- read.csv("16S_Cryo/Phyloseq_rarefied_data/220913_Cryo16s_taxonomyfile.csv", colClasses = "character", row.names = 1)


# Create matrix from data frame
# ASV=ASV frequency table and BAC=Bacteria / TAX=Taxonomy information
cryo16SASVTAB<- data.matrix(cryo16S_asv_rare)
cryo16STAX<- as.matrix(cryo16S_tax)


#Tell phyloseq what is what
cryo16SASVTAB = otu_table(cryo16SASVTAB, taxa_are_rows=TRUE)
cryo16STAX = tax_table(cryo16STAX)
cryo16SMETA = sample_data(cryo16S_meta)

#Integrating all data into a single phyloseq object
## ALL=all data integrated - BAC=Bacteria
cryo16S_phy_rare = phyloseq(cryo16SASVTAB,cryo16STAX, cryo16SMETA)
#Visualazing the Phyloseq object
cryo16S_phy_rare


#Exploring
## Number of samples
nsamples(cryo16S_phy_rare)
## Number of taxa and rank names
ntaxa(cryo16S_phy_rare)
rank_names(cryo16S_phy_rare)
## Listing sample metadata variables
sample_variables(cryo16S_phy_rare)

# COUNTING NUMBER OF ASVs BY TRANSFORMING THE FINAL TABLE TO PRESENCE/ABSENCE DATA
ITS_PresAbs_rare <- transform_sample_counts(ITS_phy_rare, function(abund) 1*(abund>0))
sample_sums(ITS_PresAbs_rare)

C16S_PresAbs_rare <- transform_sample_counts(cryo16S_phy_rare, function(abund) 1*(abund>0))
sample_sums(C16S_PresAbs_rare)

#merge_by_type
merged16S_PA = merge_samples(cryo16S_phy_rare, "host_animal")
merged16S_PA_2 <- transform_sample_counts(merged16S_PA, function(abund) 1*(abund>0))

sample_sums(merged16S_PA_2)

#number of samples per host_animal
summary(as.factor(as.data.frame(phyloseq::sample_data(cryo16S_phy_rare))$host_animal))

#number of samples per host_genus
summary(as.factor(as.data.frame(phyloseq::sample_data(cryo16S_phy_rare))$host_genus))

#Exploring phyloseq object
summarize_phyloseq(cryo16S_phy_rare)