
#Loading packages
library(phyloseq)
library(ggplot2)
library(fantaxtic)

#Same phyloseq objects
########### ITS abundance tables_ most abundant Genera  ######
genus_ITS_rare<-tax_glom(ITS_phy_rare, taxrank = "Genus")


##library(fantaxtic)
ITS_asv <- top_taxa(genus_ITS_rare, n_taxa = 20)
Genus20_ITS_rare<-plot_nested_bar(ps_obj = ITS_asv$ps_obj,
                                  top_level = "Phylum",
                                  nested_level = "Genus",
                                  palette = c(p__Basidiobolomycota = "#0A9396", 
                                              p__Mortierellomycota = "#AE2012",
                                              p__Ascomycota = "#D9ED92",
                                              p__Basidiomycota = "#FFFAE5")) +
  facet_wrap(~host_animal, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 12, 
                                  face = "bold"),
        legend.position = "right",
        legend.key.size = unit(12, "points"),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Genus20_ITS_rare +
  theme(text = element_text("Arial"))



########### 16S abundance tables_ most abundant Genera  ######

genus_16S_rare<-tax_glom(cryo16S_phy_rare, taxrank = "genus")

C16S_asv <- top_taxa(genus_16S_rare, n = 20)
Genus20_16S_rare<-plot_nested_bar(ps_obj = C16S_asv$ps_obj,
                                  top_level = "phylum",
                                  nested_level = "genus",
                                  palette = c(Proteobacteria = "#FFF75E",
                                              Campylobacterota = "#F77F00",
                                              Fusobacteriota = "#FF006E")) + 
  facet_wrap(~host_animal, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 12, 
                                  face = "bold"),
        legend.position = "right",
        legend.key.size = unit(12, "points"),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Genus20_16S_rare +
  theme(text = element_text("Arial"))



Genus20_ITS_rare + Genus20_16S_rare + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 1, nrow = 2)

###################################################

State_ITS_20Genus<-plot_nested_bar(ps_obj = ITS_asv$ps_obj,
                                  top_level = "Phylum",
                                  nested_level = "Genus",
                                  palette = c(p__Basidiobolomycota = "#0A9396", 
                                              p__Mortierellomycota = "#AE2012",
                                              p__Ascomycota = "Green",
                                              p__Basidiomycota = "#FFFAE5")) +
  facet_grid(~host_animal + state, scales = "free_x", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 12, 
                                  face = "bold"),
        legend.position = "right",
        legend.key.size = unit(12, "points"),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




##
State_16S_20Genus<-plot_nested_bar(ps_obj = C16S_asv$ps_obj,
                                  top_level = "phylum",
                                  nested_level = "genus",
                                  palette = c(Proteobacteria = "#FFF75E",
                                              Campylobacterota = "#F77F00",
                                              Fusobacteriota = "#FF006E")) + 
  facet_grid(~host_animal + state, scales = "free_x", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 12, 
                                  face = "bold"),
        legend.position = "right",
        legend.key.size = unit(12, "points"),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


State_ITS_20Genus + State_16S_20Genus + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 1, nrow = 2)


######## Extract data

write.table(genus_ITS_rare %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(host_animal, state, Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "ITS.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(ITS_phy_rare %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "ITS.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)



######## Extract data from the 20 most abundant genera

write.table(ITS_asv$ps_obj %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(host_animal, state, Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "DATA_ITS.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


            
c16S_genus <- C16S_asv$ps_obj %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>%
              psmelt() %>% select(OTU, genus, Sample, Abundance, host_animal, state) %>% 
              spread(Sample, Abundance) 


write.table(c16S_genus, file = "DATA_16S_relative_abundance_genus.tsv", sep = "\t", quote = F, 
                        row.names = F, col.names = T)
            
            
            
            
            
            
            
            
            