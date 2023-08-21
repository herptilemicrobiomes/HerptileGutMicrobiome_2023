


## load packages
library(ggtree)
library(ape)
library(phylotools)
library(phytools)

setwd("~/Ancestral_rec_R/geo")

## set seed for reproducibility
set.seed(100)

TREE_99 <- read.nexus("tree_only_STP_UHM.nexus")
TREE_99 

metdata_99_OTUs <- read.csv("metadata_only_STP_UHM.csv", row.names=1,stringsAsFactors=TRUE)


#Metadata hosts
Hosts<-setNames(metdata_99_OTUs$Host,
                      rownames(metdata_99_OTUs))

head(Hosts)


#Metadata ecoregions III 
ecoregionIII<-setNames(metdata_99_OTUs$ecoregion_III,
                      rownames(metdata_99_OTUs))

head(ecoregionIII)



# Ancestral model

#Hosts
fitHosts <- rerootingMethod(TREE_99, Hosts, model="SYM")
fitHosts
fitHosts$loglik


#Ecoregion III
fitEcoIII <- rerootingMethod(TREE_99, ecoregionIII, model="SYM")
fitEcoIII
fitEcoIII$loglik

rotateNodes<-function(tree,nodes,polytom=c(1,2),...){
  n<-length(tree$tip.label)
  if(nodes[1]=="all") nodes<-1:tree$Nnode+n
  for(i in 1:length(nodes))
    tree<-rotate(tree,nodes[i],polytom)
  if(hasArg(reversible)) reversible<-list(...)$reversible
  else reversible<-TRUE
  if(reversible){
    ii<-which(tree$edge[,2]<=n)
    jj<-tree$edge[ii,2]
    tree$edge[ii,2]<-1:n
    tree$tip.label<-tree$tip.label[jj]
  }
  return(tree)
}


#Rotation of nodes
#tree_99_root<-rotateNodes(TREE_99,"all")


#HOSTS
cols_Hosts <-  hcl.colors(length(levels(Hosts)), "Temps")
plotTree(TREE_99, fsize=1, ftype="reg",lwd=1.5, offset = 2) #, node.numbers=TRUE)
nodelabels(pie = fitHosts$marginal.anc, piecol = cols, cex=0.6)
tiplabels(pie=to.matrix(Hosts[TREE_99$tip.label],
                        levels(Hosts)),piecol=cols,cex=0.1)
legend("topleft",levels(Hosts),pch=16,
       col=cols,bty="n",cex=0.8,
       pt.cex=1.5)



#Ecoregion III
cols_EcoIII <-  hcl.colors(length(levels(ecoregionIII)), "Spectral")
plotTree(TREE_99, fsize=1, ftype="reg",lwd=1.5, offset = 2) #, node.numbers=TRUE)
nodelabels(pie = fitEcoIII$marginal.anc, piecol = cols_EcoIII, cex=0.6)
tiplabels(pie=to.matrix(ecoregionIII[TREE_99$tip.label],
                        levels(ecoregionIII)),piecol=cols_EcoIII,cex=0.1)
legend("topleft",levels(ecoregionIII),pch=16,
       col=cols_EcoIII,bty="n",cex=0.8,
       pt.cex=1.5)




