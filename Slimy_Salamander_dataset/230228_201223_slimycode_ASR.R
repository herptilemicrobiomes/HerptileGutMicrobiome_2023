# Prepare coding environment & load data #####

#load packages 
suppressPackageStartupMessages({
  library(gridExtra)
  library(vegan)
  library(ape)
  library(gdata)
  library(indicspecies)
  library(robCompositions)
  library(zCompositions)
  library(betapart)
  library(indicspecies)
  library(patchwork)
  library(matrixStats)
  library(ggtext)
  library(cowplot)
  library(RColorBrewer)
  library(scales)
  suppressMessages(library(tidyverse))
})

#set random seed
set.seed(8765)

#set working directory
setwd("C:/Users/aromer/Downloads/201223_analysis")

#prioritize tidyverse functions 
select <- dplyr::select

#read in bacterial microbiome data and subset gut bacteria
bact.gut <- fread("bact_full.csv", header=TRUE) %>% 
    filter(sample == "gut") %>%  
    column_to_rownames("group") 

# Prepare 16s data for downstream analyses #####

#parse metadata and OTU data into seperate objects 
meta <- select(bact.gut, -starts_with("Otu"))
otus <- select(bact.gut, starts_with("Otu")) %>%
  #Ensure that otu matrix doesn't contain any zero sum columns 
  select_if(~sum(.) > 0) 

#calculate richness for OTU table 
otus %>% 
  decostand("pa") %>% # convert to presence/absence
  rowSums() -> meta$richness.bac #export richness to metadata 
summary(meta$richness.bac) #examine summary statistics

#Create OTU table where minimum read count for an OTU is 11
otus.abund <- otus %>% 
  select_if(~sum(.) > 10) 
#calculate richness for this OTU table 
otus.abund %>% 
  decostand("pa") %>% # convert to presence/absence
  rowSums() -> meta$richness.bac.abund #export richness to metadata 
summary(meta$richness.bac.abund) #examine summary statistics 

# Load and process ITS data ####

#Import Fungal OTU (FOtu) table (e.g., shared file)
otu.fungi <- read.csv("ITSstability.subsample.shared_FungalOTUs.csv", header=TRUE) %>% 
  column_to_rownames("group") %>% #store sample ids as rownames
  select(starts_with("FOtu")) #retain only FOtu columns 

#Import Taxonomy of each OTU (e.g., taxonomy file)
tax.fungi <- read.csv("ITSstability.cons.taxonomy.csv", header=TRUE) %>% 
  mutate(OTU = paste("F", OTU, sep = ""), #use same naming convention for fungal OTUs (FOtu)
         taxonomy = str_extract_all(Taxonomy,  "(?<=__)[A-Za-z_]*"), #neatly extract taxonomic data 
         .keep = "unused") %>% 
  #take taxonomy data and store in tidy format
  unnest(cols = taxonomy) %>% 
  mutate(tax.level = rep(str_to_lower(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
                         length.out = nrow(.))
  ) %>% 
  pivot_wider(names_from = tax.level,
              values_from = taxonomy) %>% 
  #only retain entries for OTUs that are in OTU matrix 
  filter(OTU %in% names(otu.fungi)) %>% 
  #remove size column and store Otu names as rownames 
  select(-Size) %>% 
  column_to_rownames("OTU")

#read in fungal metadata
meta.fungi <- read.csv("ITSstability.subsample.shared_FungalOTUs.csv", header=TRUE) %>% 
  column_to_rownames("group") %>% #store sample ids as rownames
  select(-c(starts_with("FOtu"), label, numOtus)) #remove all FOtu columns 

#Generate OTU matrices with Basidiobolus OTUs
#Include all Basidiobolus OTUs
basid.fouts <- tax.fungi %>% 
  filter(genus == "Basidiobolus") %>% 
  rownames()
fungdat.basid <- otu.fungi %>%
  select(any_of(basid.fouts))
#remove basidiobolus OTUs less with than 10 reads
fungdat.basid.abund <- fungdat.basid %>% 
  select_if(~sum(.) > 10)

#Report basic summary statistics about dataset
paste(
paste(ncol(otu.fungi), "total fungal OTUs"),
paste(ncol(fungdat.basid), "total basidiobolus OTUs"),
paste(ncol(fungdat.basid.abund), "abundant basidiobolus OTUs"),
sep = "\n"
) %>% cat()

# Indicator Otus based on host clade ####

#convert matrix of abundant 16s OTUs to presence/absence
otus.abund.pa <- decostand(otus.abund, "pa") 


#indicator sp. analysis using relative abundance data with species as grouping variable and all combinations considered 
otus.abund.indic <- multipatt(otus.abund, meta$species, func = "IndVal.g", duleg = FALSE, control = how(nperm=999))
summary(otus.abund.indic)#124 bacterial OTUs selected - use these results in comparison of bacterial to Basidiobolus with 'indpower' below

#indicator sp. analysis using relative abundance data using species as grouping variable and no combinations considered 
otus.abund.nocombo_indic <- multipatt(otus.abund, meta$species, func = "IndVal.g", duleg = TRUE, control = how(nperm=999))
summary(otus.abund.nocombo_indic)



#Assess if bacterial OTUs have preferences for all possible salamander groupings - based on presence absence  - 
#function "r.g" - correlation based approach that differs from indicator species analysis

#all combinations considered
otus.abund.pa.rg <- multipatt(otus.abund.pa, bact.gut$species, func = "r.g", duleg = FALSE, control = how(nperm=999)) 
summary(otus.abund.pa.rg)

#no combinations considered
otus.abund.pa.nocombo_rg <- multipatt(otus.abund.pa, meta$species, func = "r.g", duleg = TRUE, control = how(nperm=999)) 
summary(otus.abund.pa.nocombo_rg)


# Indicator Basidiobolus FOtus based on host clade ####

#convert relative abundances of abundant Basidiobolus to presence/absence 
fungdat.basid.abund.pa <- decostand(fungdat.basid.abund, "pa") 


#indicator sp. analysis using relative abundance data with species as grouping variable and all combinations considered 
basid.abund.indic <- multipatt(fungdat.basid.abund, meta.fungi$species, func = "IndVal.g", duleg = FALSE, control = how(nperm=999))
summary(basid.abund.indic) 

#indicator sp. analysis using relative abundance data using species as grouping variable and no combinations considered 
basid.abund.nocombo_indic <- multipatt(fungdat.basid.abund, meta.fungi$species, func = "IndVal.g", duleg = TRUE, control = how(nperm=999))
summary(basid.abund.nocombo_indic) 



#Assess if Basidiobolus FOtus have preferences for salamander host species - based on presence absence  - 
#function "r.g" - correlation based approach that differs from indicator species analysis

#all combinations considered
basid.abund.pa.rg <- multipatt(fungdat.basid.abund.pa, meta.fungi$species, func = "r.g", duleg = FALSE, control = how(nperm=999)) 
summary(basid.abund.pa.rg)

#no combinations considered
basid.abund.pa.nocombo_rg <- multipatt(fungdat.basid.abund.pa, meta.fungi$species, func = "r.g", duleg = TRUE, control = how(nperm=999)) 
summary(basid.abund.pa.nocombo_rg)


# correlation heat map ####

#Select top 100 most abundant 16S Otus and sort by abundance
bact.100 <- otus.abund[,names(sort(colSums(otus.abund), decreasing = TRUE))] %>%
  names() %>%
  .[1:100]
otus.abund.100 <- otus.abund %>% #generate dataframe of just these Otus
  select(any_of(bact.100))

#merge 100 top abundant 16s Otus and Abundant Basidiobolus FOtus by row name
abund.otus.fouts <- merge(otus.abund.100, fungdat.basid.abund, by = 0) %>%
  column_to_rownames("Row.names")

#convert to presence/absence
abund.otus.fouts.pa <- decostand(abund.otus.fouts, "pa")

#perform indicator power analysis using methodology from Halme et al. 2009 (https://doi.org/10.1111/j.1523-1739.2009.01206.x)
gut.bact.100.basid.pa.indpower <- indpower(abund.otus.fouts.pa)
#clean-up out distance matrix
diag(gut.bact.100.basid.pa.indpower) <- NA
gut.bact.100.basid.pa.indpower %>%
  as.data.frame() %>%
  select(-contains(".Otu")) %>%
  rownames_to_column("id") %>%
  filter(!grepl(".FOtu", id)) %>%
  mutate(id = str_replace_all(id, "i.", ""),
         id = str_replace_all(id, "(?<=Otu)0*", " "),) %>%
  column_to_rownames("id") %>%
  drop_na() %>%
  rename_all(.funs = ~str_replace(.x, "t.", "Basid. ")) %>%
  rename_all(.funs = ~str_replace(.x, "(?<=Otu)0*", " ")) %>%
  as.matrix() %>%
  t() ->  gut.bact.100.basid.pa.indpower

#generate correlation matrix for basidiobolus vs abundant 16s OTUs 
otus.abund.100.trim <- otus.abund.100[rownames(otus.abund.100) %in% rownames(fungdat.basid.abund),] 
cor.matrix <- cor(as.matrix(otus.abund.100.trim), as.matrix(fungdat.basid.abund))
cor.df <- cor.matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("bact") %>% 
  pivot_longer(cols = -bact, names_to = "basid", values_to = "cor") %>%  
  mutate(across(1:2, ~as.factor(.x)))

#using Hierarchical Clustering to arrange samples in heatmap
otu.levels <- rownames(cor.matrix)[hclust(dist(cor.matrix))$order]
fotu.levels <- colnames(cor.matrix)[hclust(dist(t(cor.matrix)))$order]

#color code x-axis labels according to basidiobolus species 
tax.fungi %>%
  filter(rownames(.) %in% names(fungdat.basid.abund)) %>% 
  mutate(species = str_extract(species, "(?<=_)[a-z]*")) %>% 
  mutate(species = fct_recode(species, `sp.` = "sp", `sp.` = "unclassified"),
         basid.species = paste(genus, species)) %>% 
  arrange(basid.species) %>% 
  rownames_to_column("fotu") %>% 
  mutate(label.col = ifelse(species == "heterosporus", "dodgerblue2",
                            ifelse(species == "ranarum", "red",
                                   ifelse(species == "sp.", 'gray40', NA))),
         plot.label = paste("<span style = 'color:", label.col, "'>", fotu, "</span>", sep = ""),
         plot.label = str_replace(plot.label, "FOtu[0]+", "Otu ")) %>% 
  select(fotu, plot.label) -> fotu.labels

#color code y-axis labels according to 16s OTU class 
fread("Bacteria.stability.an.cons.taxonomy.csv") %>% 
  filter(OTU %in% colnames(otus.abund.100.trim)) %>% 
  select(-Size) %>% 
  separate(Taxonomy,
           sep = ";",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           extra = "drop") %>% 
  rename(otu = OTU) %>% 
  mutate(across(.cols = -otu,
                .fns = ~ifelse(.x == "unclassified",
                               "unclassified",
                               map_chr(str_extract_all(.x, "[a-zA-Z_1-9\\s]*(?=[\\(\"])"), ~ str_c(.x, collapse=""))))) %>% 
  select(otu, class = Class) -> otu.abund.class
class.col <- otu.abund.class %>%
  group_by(class) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  ungroup() %>% 
  mutate(col.label = c(rep(T, 3), rep(F, nrow(.)-3)),
         col = ifelse(col.label == F, "gray40", brewer.pal(3, "Set1"))) %>% 
  select(class, col)
left_join(otu.abund.class, class.col) %>%
  mutate(otu.label = paste("<span style = 'color:", col, "'>", otu, "</span>", sep = "")) %>% 
  select(bact = otu, otu.label) %>% 
  left_join(cor.df, .) -> cor.df
cor.df %>% 
  left_join(fotu.labels, by = c("basid" = "fotu")) %>% 
  mutate(basid = factor(basid, levels = fotu.levels),
         bact = factor(bact, levels = otu.levels),
         otu.label = fct_reorder(str_replace_all(otu.label, "Otu[0]*", "OTU "), as.numeric(bact)),
         plot.label = fct_reorder(str_replace_all(plot.label, "Otu", "OTU"), as.numeric(basid))) %>% 
  ggplot(aes(x = plot.label, y = otu.label, fill = cor)) + 
  geom_tile(color = "gray50") + 
  scale_fill_gradientn(colors = c("red", "white", "royalblue"),
                       values = rescale(c(min(cor.matrix), 0, max(cor.matrix))),
                       guide = guide_colorbar(
                         frame.colour = "black",
                         ticks = TRUE, 
                         nbin = 10,
                         label.position = "top",
                         title.vjust = 0.25,
                         barwidth = 13,
                         barheight = 1.3, 
                         direction = 'horizontal'
                       )) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() + 
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 12.5),
        axis.text.y = element_markdown(hjust = 0, size = 7.5, color = "black"),
        axis.title.x = element_markdown(margin = margin(10,0,5,0), size = 15),
        axis.title.y = element_markdown(margin = margin(0,10,0,5), size = 15),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 15),
        legend.position = "top") + 
  labs(x = "*Basidiobolus* OTUs",
       y = "16S OTUs",
       fill = "Correlation") -> cor.hm

# Examine indicator power by bacterial taxonomy #####

#generate dataframe for plot
gut.bact.100.basid.pa.indpower %>%
as.data.frame() %>%
rownames_to_column("basid") %>%
pivot_longer(cols = -basid, names_to = "bacteria", values_to = "indpower") -> indic.df

#generate mean value indicator power for each bacterial OTU
indic.df %>%
  group_by(bacteria) %>% 
  summarize(indpower = mean(indpower)) %>%
  ungroup() -> mean.indic.df

#load 16s taxonomy data 
fread("Bacteria.stability.an.cons.taxonomy.csv") %>% 
  select(-Size) %>% 
  separate(Taxonomy,
           sep = ";",
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           extra = "drop") %>% 
  mutate(otu = str_replace_all(OTU, "(?<=Otu)0*", " "),
         .keep = "unused",
         .before = 1) %>%
  mutate(across(.cols = -otu,
                .fns = ~ifelse(.x == "unclassified",
                               "unclassified",
                               map_chr(str_extract_all(.x, "[a-zA-Z_1-9\\s]*(?=[\\(\"])"), ~ str_c(.x, collapse=""))))) %>% 
  filter(otu %in% str_replace_all(names(otus), "(?<=Otu)0*", " ")) -> tax.bac

#join indicator taxa with taxonomy data 
left_join(mean.indic.df,
          tax.bac,
          by = c('bacteria' = "otu")) -> tax.mean.indic.df

#generate figure of average indicator power of bacterial OTU per class
tax.mean.indic.df %>%
  group_by(Class) %>%
  mutate(indpower.mean = mean(indpower),
         se = sd(indpower)/sqrt(length(indpower)),
         n = n(),
         n.y = ifelse(n == 1,
                      indpower.mean + .015,
                      indpower.mean + se + .012)) %>%
  ungroup() %>%
  rename(class = Class) %>% 
  left_join(class.col) %>% 
  mutate(class = fct_recode(class, Unclassified = "unclassified", `Spirochaetia` = "Spirochaetes")) %>% 
  mutate(class.label = paste("<span style = 'color:", col, "'>", class, "</span>", sep = "")) %>% 
  mutate(class = fct_reorder(class, indpower.mean, .desc = T),
         class.label = fct_reorder(class.label, as.numeric(class))) %>%
  ggplot(aes(x = class.label, y = indpower.mean)) +
  geom_errorbar(aes(ymin = indpower.mean - se,
                    ymax = indpower.mean + se),
                width = 0.15,
                alpha = 0.75) +
  geom_point(aes(fill = indpower.mean), shape = 21, color = "black", size = 5) +
  scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(4, "Blues"))) +
  scale_y_continuous(limits = c(0.30, 0.5)) +
  geom_text(aes(label = n, y = n.y),
            size = 7.5,
            check_overlap = T) +
  theme_classic() +
  labs(y = "Indicator Power<br><span style = 'font-size:12pt;'>(Average per Class &plusmn; SE)</span>",
       x = "Bacterial Class") +
  theme(legend.position = "none",
        axis.title.y = element_markdown(size = 18, margin = margin(r = 10)),
        plot.margin = margin(20,85,20,20),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_markdown(angle = 35,
                                       hjust= 1,
                                       size = 15,
                                       margin = margin(b = 10, r = 10)),
        axis.text.y = element_text(size = 14)
  ) -> bacotu.indpower
bacotu.indpower

# Examine total indicator power of abundant of basidiobolus OTUs ####

# Store TIP values and their variance in a dataframe
data.frame(
  fotu = fotu.labels$fotu,
  tip = rowMeans(gut.bact.100.basid.pa.indpower, na.rm=TRUE),
  var = rowVars(gut.bact.100.basid.pa.indpower, na.rm=TRUE)
) %>%
  remove_rownames() -> tip.df

#store label data in a dataframe
label.df <- data.frame(
  label = "<span style = 'color:red'>*B. ranarum*</span>
  <br><span style = 'color:dodgerblue2'>*B. heterosporus*</span>
  <br><span style = 'color:gray40'>*Basidiobolus* spp.</span>",
  x = 4.5,
  y = .9
)

#Generate plot of TIP Values
tip.df %>%
  left_join(fotu.labels) %>%
  mutate(plot.label = str_replace_all(plot.label, "Otu", "OTU"),
         plot.label = fct_reorder(plot.label, tip, .desc = T)) %>%
  ggplot(aes(x = plot.label, y = tip)) +
  geom_errorbar(aes(ymin = tip + var, ymax = tip - var), width = 0.5, linewidth = 1) +
  geom_col(aes(fill = tip), color = "black", linewidth = 1, show.legend = F) +
  geom_textbox(data = label.df, aes(label = label, y = y, x = x), width = unit(3, "inch"), size = 6.5, box.color = NA) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_gradientn(colors = c("White", brewer.pal(8, "Blues")[c(1,5,8)])) +
  theme_classic() +
  labs(y = "Total Indicator Power<br><span style = 'font-size:12pt;'>(Average &plusmn; variance)</span>",
       x = "*Basidiobolous* OTUs") +
  theme(legend.position = "none",
        axis.title.y = element_markdown(size = 18, margin = margin(r = 10)),
        plot.margin = margin(20,85,20,20),
        axis.title.x = element_markdown(size = 18, margin = margin(0)),
        axis.text.x = element_markdown(angle = 45,
                                       hjust = 1,
                                       size = 15,
                                       margin = margin(b = 10, r = 10, t = 5)),
        axis.text.y = element_text(size = 14)) -> basidfout.tip
basidfout.tip

# Generate multipanel figure for perspectives paper ####

plot_grid(
  cor.hm,
  plot_grid(
    bacotu.indpower,
    basidfout.tip,
    nrow = 2,
    labels = c("B", "C"),
    rel_heights = c(1, 0.8),
    label_size = 30,
    label_x = 0.01
  ),
  labels = c("A", ""),
  rel_widths = c(0.8, 1),
  label_size = 30,
  label_x = 0.02
) -> per.multi

# ggsave(filename = "ind.multipanel.png",
#        plot = per.multi,
#        device = "png",
#        width = 6000,
#        height = 4000,
#        units = "px"
#        )
