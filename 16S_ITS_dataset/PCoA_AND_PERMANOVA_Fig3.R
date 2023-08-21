
#USE SAME PHYLOSEQ OBJECTS
#PCoA
#Bray and Curtis distance
ITS_dist = phyloseq::distance(ITS_phy_rare, method="bray")
ordination_ITS = ordinate(ITS_phy_rare, method="PCoA", distance=ITS_dist)
PCoA_ITSB <- plot_ordination(ITS_phy_rare, ordination_ITS) + 
  geom_point(aes(color = state, shape = host_animal), size=2) +
  ggtitle("ITS")+
  stat_ellipse(aes(linetype=host_animal)) +
  scale_linetype_manual(values= c("F1","dashed", "dotted"))+
  theme(text = element_text(size = 12)) + labs(color = "State", shape = "Host") + 
  scale_color_manual(values = c("#E63946", "#1F968BFF", "#A7C957", "#386641", "#8338EC", "#FFB703", "#FB8500", "#4CC9F0"))

#Centroids ITS-PCoA

#Extract sample data from phyloseq object

ord.meta_ITS <- data.frame(sample_data(ITS_phy_rare))

bc.envfit.pcoa_ITS <- envfit(ordination_ITS$vectors[,1:2] ~ host_animal,
                         data = ord.meta_ITS,
                         display = 'sites',
                         na.rm = T)

centroids.pcoa_ITS <- as.data.frame(bc.envfit.pcoa_ITS$factors$centroids)

centroids.pcoa_ITS$Label <- c('Frog', 'Lizard', 'Salamander')

PCoA_ITSB_C <- PCoA_ITSB +
  geom_label(aes(label = Label,
                 x = Axis.1,
                 y = Axis.2,
                 alpha = 0.90),
             show.legend = F,
             data = centroids.pcoa_ITS,
             inherit.aes = F) +
  theme_bw() +
  background_grid(major = 'none', minor = 'none')



 


#16S PCoA and centroids

c16S_dist = phyloseq::distance(cryo16S_phy_rare, method="bray")
ordination_16S = ordinate(cryo16S_phy_rare, method="PCoA", distance=c16S_dist)
PCoA_16SB <- plot_ordination(cryo16S_phy_rare, ordination_16S) + 
  geom_point(aes(color = state, shape = host_animal), size=2) +
  ggtitle("16S")+
  stat_ellipse(aes(linetype=host_animal)) +
  scale_linetype_manual(values= c("F1","dashed", "dotted"))+
  theme(text = element_text(size = 12)) + labs(color = "State", shape = "Host") + 
  scale_color_manual(values = c("#E63946", "#1F968BFF", "#A7C957", "#386641", "#8338EC", "#FFB703", "#FB8500", "#4CC9F0"))

#Centroids

ord.meta_16s <- data.frame(sample_data(cryo16S_phy_rare))

bc.envfit.pcoa.16S <- envfit(ordination_16S$vectors[,1:2] ~ host_animal,
                         data = ord.meta_16s,
                         display = 'sites',
                         na.rm = T)

centroids.pcoa.16S <- as.data.frame(bc.envfit.pcoa.16S$factors$centroids)

centroids.pcoa.16S$Label <- c('Frog', 'Lizard', 'Salamander')


PCoA_16SB_C <- PCoA_16SB +
  geom_label(aes(label = Label,
                 x = Axis.1,
                 y = Axis.2,
                 alpha = 0.90),
             show.legend = F,
             data = centroids.pcoa.16S,
             inherit.aes = F) +
  theme_bw() +
  background_grid(major = 'none', minor = 'none')

PCoA_ITSB_C + PCoA_16SB_C + plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect")


################################################################################################################################################



#PERMANOVA/ADONIS
metadata_ITS <- data.frame(sample_data(ITS_phy_rare))
test.adonis_ITS <- adonis(ITS_dist ~ host_animal*state, data = metadata_ITS)
test.adonis_ITS <- as.data.frame(test.adonis_ITS$aov.tab)
test.adonis_ITS

fungal_model <- adonis2(ITS_dist ~ host_animal*state, data = metadata_ITS, permutations=999, strata = metadata_ITS$state)


metadata_16S <- data.frame(sample_data(cryo16S_phy_rare))
test.adonis_16S <- adonis(c16S_dist ~ host_animal, data = metadata_16S)
test.adonis_16S <- as.data.frame(test.adonis_16S$aov.tab)
test.adonis_16S


bacterial_model <- adonis2(c16S_dist ~ host_animal*state, data = metadata_16S, permutations=999, strata = metadata_16S$state)
test.adonis_16S <- as.data.frame(bacterial_model$aov.tab)




metadata_ITS <- data.frame(sample_data(ITS_phy_rare))
test.adonis_ITSJ <- adonis(ITS_distJ ~ host_animal, data = metadata_ITS)
test.adonis_ITSJ <- as.data.frame(test.adonis_ITSJ$aov.tab)
test.adonis_ITSJ


metadata_16S <- data.frame(sample_data(cryo16S_phy_rare))
test.adonis_16SJ <- adonis(c16S_distJ ~ host_animal, data = metadata_16S)
test.adonis_16SJ <- as.data.frame(test.adonis_16SJ$aov.tab)
test.adonis_16SJ

#PAIRWISE PERMANOVA

cbn_ITS <- combn(x=unique(metadata_ITS$host_animal), m = 2)
p <- c()

for(i in 1:ncol(cbn_ITS)){
  ps.subs <- subset_samples(ITS_phy_rare, host_animal %in% cbn_ITS[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ host_animal, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn_ITS), p=p, p.adj=p.adj)
p.table


cbn_16S <- combn(x=unique(metadata_16S$host_animal), m = 2)
p <- c()

for(i in 1:ncol(cbn_16S)){
  ps.subs <- subset_samples(cryo16S_phy_rare, host_animal %in% cbn_16S[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "jaccard") ~ host_animal, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table16 <- cbind.data.frame(t(cbn_16S), p=p, p.adj=p.adj)
p.table16
