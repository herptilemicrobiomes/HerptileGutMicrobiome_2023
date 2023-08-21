library(tidyverse)
library(factoextra) 
library(mmtable2)

setwd("G:/My Drive/Basidobolus/ASV analysis")

#import and format data
df <- read.csv(file = "ASVcount_and_Tax_ITSBasi.csv", row.names = 1) %>% 
  select(-starts_with("tax.")) %>% 
  rename_all(~sapply(str_split(string = ., pattern = "_"), `[`, 1)) %>% 
  t() %>% 
  as.data.frame()
basido.asvs <- read.csv(file = "ASVcount_and_Tax_ITSBasi.csv", row.names = 1) %>%
  filter(tax.Genus == "g__Basidiobolus") %>%
  rownames_to_column(var = "asv") %>% 
  select(asv) %>% pull()

#count number of asvs per isolate
df.n_asvs <- df %>% 
  mutate(across(everything(), ~ifelse(. > 0, 1, 0))) %>%
  mutate(asv.count = rowSums(.),
         basido.asv.count = rowSums(.[,basido.asvs])) %>%
  rownames_to_column(var = "isolate") %>% 
  select(-starts_with("ASV_")) %>% 
  arrange(asv.count)
# visualize this result
df.n_asvs %>% 
  mutate(isolate = fct_reorder(isolate, basido.asv.count, .desc = T)) %>%
  ggplot(aes(x = isolate, y = basido.asv.count, fill = isolate)) + 
  geom_col(color = "black")  + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20), breaks = seq(0,18,3)) + 
  labs(y = "ASVs of Basidiobolus", 
       fill = "Fungal Isolate") + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(0.8, 0.8)
  )

#scale and normalize ASV counts
means <- apply(df,2,mean)
sds <- apply(df,2,sd)
nor <- scale(df,center=means,scale=sds)
#search for optimal number of clusters
factoextra::fviz_nbclust(nor, FUNcluster = kmeans, method = "silhouette", k.max = (nrow(nor)-1))
#perform k-means clusterings
mykmeans <- kmeans(nor, centers = 4, nstart = 25)
#plot kmeans clustering output
fviz_cluster(mykmeans, data = df, labelsize = 15) + 
  scale_x_continuous(limits = c(-9,9)) +
  labs(title = "K-means Clustering") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.margin = margin(10,10,10,10))

#extract k means results 
kmeans.results <- df %>%
  mutate(cluster = mykmeans$cluster) %>% 
  rownames_to_column(var = "isolate") %>% 
  select(isolate, cluster)  
  
#Copy data to clipboard to generate a table 
right_join(kmeans.results, df.n_asvs, by = "isolate") %>% 
  arrange(cluster, basido.asv.count) %>% 
  write.table(., "clipboard", sep="\t", row.names=FALSE)
