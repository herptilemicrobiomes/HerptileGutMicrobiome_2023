#TESTS GEO AND HOSTS - selected 99 OTUs


setwd("~/Ancestral_rec_R/geo")

metadata <- read.csv("metadata.csv", header=T, row.names=1, check.names=F)

distmatrix  <- read.csv("Patristics_distances_STP_UHM.csv", header=T, row.names=1, check.names=F)
distmatrix[upper.tri(distmatrix)] <- NA
m <- as.dist(distmatrix)



set.seed(1234)
Host_EcoIII <- adonis2(m ~ Host*ecoregion_III, data = metadata, permutations = 999)
set.seed(1234)
Host_EcoIV <- adonis2(m ~ Host*ecoregion_IV, data = metadata, permutations = 999)


set.seed(1234)
Host_Genus_EcoIII <- adonis2(m ~ Genus*ecoregion_III, data = metadata, permutations = 999)
set.seed(1234)
Host_Genus_EcoIV <- adonis2(m ~ Genus*ecoregion_IV, data = metadata, permutations = 999)



#### Betadisper analysis

#Hosts
set.seed(1234)
betadisp_host <- betadisper(m, metadata$Host)
betadisp_host 

#plotting host:
box_betadisper_host<-boxplot(betadisp_host, xlab = "", las = 1, cex.axis = 0.8)
set.seed(1234)
permutest.betadisper.host<-permutest(betadisp_host, pairwise=TRUE,
                                     permutations=999, parallel=getOption("mc.cores"),)

plot(betadisp_host)


#Ecoregion
set.seed(1234)
betadisp_EcoIII <- betadisper(m, metadata$ecoregion_III)
betadisp_EcoIII

#plotting Ecoregion:
box_betadisper_ecoIII<-boxplot(betadisp_EcoIII, xlab = "", las = 1, cex.axis = 0.7)

set.seed(1234)
permutest.betadisper.ecoIII<-permutest(betadisp_EcoIII, pairwise=TRUE,
                                     permutations=999, parallel=getOption("mc.cores"),)

plot(betadisp_EcoIII)


#hosts Genera

set.seed(1234)
betadisp_Genera <- betadisper(m, metadata$Genus)
betadisp_Genera

#plotting genera:
box_betadisper_genera<-boxplot(betadisp_Genera, xlab = "", las = 1, cex.axis = 0.7)

set.seed(1234)
permutest.betadisper.genera<-permutest(betadisp_Genera, pairwise=TRUE,
                                       permutations=999, parallel=getOption("mc.cores"),)

plot(betadisp_Genera)



