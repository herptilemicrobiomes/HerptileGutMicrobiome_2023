#load packages
library(tidyverse)
library(ggtext)
library(gdata)
library(RColorBrewer)
library(shadowtext)
library(vegan)
library(ggpubr)

#load all relevant data 
meta <- read.csv("G:/My Drive/Walker Lab/Basidobolus/Sampling Metadata/220526_MECO_cryo_metadata.csv")
tax <- read.csv("G:/My Drive/Walker Lab/Basidobolus/Sampling Metadata/220922_CryoITS_taxonomyfile.csv")
seqs <- read.csv("G:/My Drive/Walker Lab/Basidobolus/Sampling Metadata/220922_CryoITS_rarified.dataset.csv")
asvs <- read.csv("G:/My Drive/Walker Lab/Basidobolus/ASV analysis/ASVcount_and_Tax_ITSBasi.csv", row.names = 1)

#format ASV dataframe
asvs  %>% 
  filter(tax.Genus == 'g__Basidiobolus') %>% 
  select(-starts_with("tax")) %>% 
  rename_all(~str_extract(.x, '^.*?(?=\\_)')) %>% 
  select(contains("UHM") | contains('stp')) -> asvs.f

#get count of all basidiobolus ASVs for all strains
asvs.f %>% 
  mutate(across(everything(), ~ifelse(.x > 1, 1, .x))) %>% 
  colSums(.) -> asvs.count

#get count of abundant basidiobolus ASVs (> 1% relative abundance)
asvs.abund.count <- asvs.f %>%
  t() %>%
  as.data.frame() %>%
  mutate(total.reads = rowSums(.)) %>%
  mutate(across(starts_with('ASV'), ~ifelse(.x/total.reads > .01, 1, 0))) %>%
  select(-total.reads) %>%
  rowSums(.) %>%
  subset(., . > 0)

#generate dataframe of ASVs count for downstream plot 
data.frame(strain = toupper(names(asvs.count)),
           total.asvs = asvs.count,
           abund.asvs = asvs.abund.count,
           y = 1.003,
           label2 = "*",
           y2 = c( 0.15, 0.15, 0.425, 0.6, 0.4, 0.4)) %>% 
  mutate(label = paste(abund.asvs, "/", total.asvs)) -> df.asvcounts

#plot of ITS ASV abundance 
asvs.f %>%
  t() %>%
  as.data.frame() %>%
  mutate(total.reads = rowSums(.)) %>% 
  mutate(across(-total.reads, ~.x/total.reads)) %>%
  select(-total.reads) %>% 
  rownames_to_column("strain") %>% 
  pivot_longer(starts_with("ASV"))  %>% 
  mutate(name = factor(str_replace(name, "_", " "), levels = gtools::mixedsort(str_replace(rownames(asvs), "_", " "))),
         strain = toupper(strain)) %>% 
  filter(value > .01) %>% 
  ggplot(aes(x = strain, y = value)) + 
  geom_col(aes(fill = name), color = "black", linewidth = 0.75) + 
  geom_text(data = df.asvcounts, aes(x = strain, y = y, label = label), size = 10, vjust = -0.5) +
  geom_text(data = df.asvcounts, aes(x = strain, y = y2, label = label2), size = 25, color = 'white') +  
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.06)) +
  scale_fill_manual(values = rev(colorRampPalette(brewer.pal(11, "Spectral"))(13))) +
  theme_classic() +
  labs(x = "<i>Basidiobolus</i> Strain",
       y = "ITS ASVs<br><span style = 'font-size:12pt'>(Relative Abundance)</span>") +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_markdown(margin = margin(12.5,0,5,0)),
        axis.title = element_markdown(size = 15),
        axis.title.y = element_markdown(margin = margin(0,12.5,0,5)),
        plot.margin = margin(t = -50)) -> asv.plot 

#get list of Basidiobolus ASVs from the cryo dataset 
basid.asvs <- tax %>% 
  filter(genus == "g__Basidiobolus") %>% 
  pull(ASV)

#generate dataframe of just basidiobolus sequences from cryo dataset 
seqs.basid <- seqs %>%
  select(any_of(basid.asvs))

#generate dataframe of basidiobolus richness and relative abundance for downstream plot 
df.basid.cryo <- data.frame(
  sample.ID = str_extract(seqs$sample.id, '^.*?(?=\\.)'),
  basid.richness = rowSums(decostand(seqs.basid, 'pa')),
  basid.reads = rowSums(seqs.basid),
  depth = rowSums(seqs[,-c(1:2)])
) %>% 
  mutate(basid.ra = basid.reads/depth) %>% 
  left_join(meta) %>% 
  rename(sample.id = 'sample.ID') %>% 
  select(sample.id, host_animal, depth, basid.reads, basid.ra, basid.richness) %>% 
  mutate(host_animal = str_to_sentence(host_animal)) 

#generate plot of basidiobolus richness by host taxon 
df.basid.cryo %>% 
  ggplot(aes(x = basid.richness, y = after_stat(scaled))) + 
  geom_density(aes(fill = host_animal), alpha = 0.65, linewidth = 0.75) + 
  labs(x = "*Basidiobolus* ITS<br><span style = 'font-size:10pt'>(n ASVs)</span>") + 
  scale_fill_manual(values = c('#fdaa89', '#ffe363', '#bce27f')) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_markdown(margin = margin(12.5,0,0,0), size = 15),
        axis.title.y = element_blank(),
        plot.margin = margin(5, 5, -100, 5)) -> basidrich.plot

#generate plot of basidiobolus relative abundance by host taxon 
df.basid.cryo %>% 
  ggplot(aes(x = basid.ra, y = after_stat(scaled))) + 
  geom_density(aes(fill = host_animal), alpha = 0.65, linewidth = 0.75) + 
  labs(x = "*Basidiobolus* Relative Abundance",
       y = "Scaled Density") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('#fdaa89', '#ffe363', '#bce27f')) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_markdown(margin = margin(12.5,0,0,0)),
        axis.title = element_markdown(size = 15),
        axis.title.y = element_markdown(margin = margin(0,12.5,0,5)),
        plot.margin = margin(5, 5, -100, 5)) -> basidabund.plot

#generate multipanel figure of basidiobolus plots 
ggarrange(basidabund.plot, basidrich.plot,
          common.legend = T,
          legend = "right",
          align = 'hv',
          labels = c("A", "B"),
          label.y = 1.025,
          label.x = 0.015) + 
  theme(plot.margin = margin(b = -100, t = 10)) -> basid.plots

#generate final multipanel 
ggarrange(basid.plots, asv.plot, 
         nrow = 2,
         align = "h",
         heights = c(0.75, 1),
         labels = c("", "C"),
         label.y = 1.025,
         label.x = 0.015
         ) 

#save high-resolution file of final multipanel
ggsave(plot = last_plot(),
       filename = "perspectives.supp1.jpg",
       device = "jpg",
       width = 3000,
       height = 3000,
       units = "px",
       path = "G:/My Drive/Walker Lab/Basidobolus/Sampling Metadata")
  

