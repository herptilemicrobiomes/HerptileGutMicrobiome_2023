#load packages 
library(lubridate)
library(decontam)
library(ggpubr)
library(tidyverse)
library(gridExtra)

#source custom functions/plotting theme
setwd("G:/My Drive/custom_r_functions")
source("gg_color_hue.R")
source("romer_theme.R")

#set working directory
setwd("G:/My Drive/Cryo_alex")

#get system date and save 
date <- format(Sys.Date(), "%y%m%d")

#specify data set file paths and name data sets 
datasetnames <- c("Wood Frog - 16s", "Wood Frog - ITS",
                  "Cryo Collection - 16s", "Cryo Collection - ITS")

datasetfiles <- c('220901_ASVcounts_and_Tax_16S_WF.csv',
                  "220901_ASVcounts_and_Tax_ITSShcmidtWF.csv",
                  "220901_ASVcount_and_Tax_16s_Cryo.csv",
                  "220901_ASVcount_and_Tax_ITS_Cryo.csv")

#load biosample / sample concentration codec
codec <- read.csv("220901_biosamplecodec.csv")

#generate list to save loop output 
decontam.out <- list() 

#loop that iterates decontam pipeline across each data set
#with multiple threshold values to determine optimal value  
for(i in 1:length(datasetfiles)){
  
  #import and name data set in loop iteration  
  df <- read.csv(datasetfiles[i])
  datasetname <- datasetnames[i]
  
  #generate taxonomy table
  df.tax <- df %>%   
    select(1, contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))) %>%
    column_to_rownames(var = 'X') %>% 
    rename_all(~str_to_lower(gsub(x = ., pattern = "tax[.]", replacement = ""))) %>% 
    rownames_to_column("ASV")
  
  #Run decontamination procedures 
  if(grepl("WF", datasetfiles[i])==T){
    
    #generate ASV abundance table  
    df.asv <- df %>%
      select(contains("UHM"), X) %>%
      column_to_rownames(var = 'X') %>%
      rename_all(~sapply(str_split(string = ., pattern = "_"), `[`, 1)) %>%
      t() %>% 
      as.data.frame() %>% 
      mutate(sample = ifelse(grepl("NTC", rownames(.)), "NTC", "Sample"),
             .before = 1) %>% 
      rownames_to_column(var = "sample.id") 
    
    #load concentration data based on marker data 
    con_values <- 
      if(grepl("16s", datasetname) == T){
        codec %>%
          filter(label.16s %in% df.asv$sample.id) %>% 
          select(label.16s, con_16s)
      } else 
        if(grepl("ITS", datasetname) == T){
          codec %>%
            filter(label.its %in% df.asv$sample.id) %>% 
            select(label.its, con_its)
        }
    
    #Add concentration values to ASV table
    names(con_values) <- c("sample.id", "sample.con")
    df.asv <- df.asv %>% 
      left_join(., con_values, by = "sample.id") %>% 
      relocate(sample.con, .after = sample.id)
    
    #specify sequence of threshold probability values to iterate through 
    threshold <- round(seq(0.05,0.95, length.out = 10), 3)
    
    #generate empty lists for loop output items 
    sapply(c("decontam.results.list", 'decontam.summary', "plots.sampletype", "plots.persample"),
           function(k) assign(k, list(), envir = .GlobalEnv))
    
    #specify dataset being process during loop
    print(datasetname)
    
    for(j in 1:length(threshold)){
      #identify contaminating taxa via decontam 
      decontam.results <- df.asv %>% 
        select(starts_with("ASV_")) %>% 
        as.matrix() %>% 
        isContaminant(seqtab = ., conc = df.asv$sample.con, method="frequency", threshold = threshold[j]) %>% 
        rownames_to_column("ASV") %>% 
        left_join(., df.tax, by = "ASV") %>% 
        filter(contaminant == T)
      #report threshold value to track loop 
      print(paste("threshold:", threshold[j], sep = ""))
      #save results into list
      decontam.results.list[[j]] <- decontam.results
      #name each dataframe according to the data set and threshold value 
      names(decontam.results.list)[j] <- paste("Threshold", threshold[j], sep=": ")
      
      #extract list of contaminant ASVs 
      contaminants <- filter(decontam.results, contaminant == T)$ASV
      
      #generate long dataframe and classify taxa according to decontam results 
      df.asv.long <- df.asv %>%
        rename(sample.type = sample) %>%
        mutate(sample = rownames(.)) %>% 
        pivot_longer(starts_with("ASV_") | starts_with("Otu"), names_to = "asv", values_to = "reads") %>% 
        mutate(asv.type = ifelse(asv %in% contaminants, "contamination", "sample taxa")) %>% 
        group_by(sample.type) %>% 
        mutate(reads = as.numeric(reads),
               total.reads = sum(reads))  %>% 
        ungroup()
      
      #plot decontam results according to sample type 
      df.plot.sampletype <- df.asv.long %>% 
        group_by(sample.type, asv.type) %>% 
        mutate(asv.type.reads = sum(reads),
               asv.proportion = asv.type.reads/total.reads) %>% 
        sample_n(1) %>% 
        mutate(asv.type = as.factor(str_to_title(asv.type)))
      plot.sampletype <- df.plot.sampletype %>% 
        ggplot(aes(x = sample.type, y = asv.proportion, fill = asv.type)) + 
        geom_col() + 
        scale_y_continuous(expand = c(0,0), labels = scales::percent) +
        scale_fill_manual(values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")) +
        labs(title = datasetname,
             subtitle = paste("Threshold =", threshold[j], sep = " ")) + 
        theme(plot.title = element_text(hjust = 0),
              plot.subtitle = element_text(size = 14))
      #extract summary values for decontamination
      decontam.summary[[j]] <- df.plot.sampletype %>% 
        mutate(threshold = threshold[j]) %>% 
        select(sample.type, threshold, asv.type, asv.proportion) %>%
        as.data.frame()
      #save plots in a named list
      plots.sampletype[[j]] <- plot.sampletype
      names(plots.sampletype)[j] <- paste("Threshold", threshold[j], sep=": ")
      #plot decontam results on a per-sample basis 
      df.persample <- df.asv.long %>%
        group_by(sample) %>% 
        mutate(sample.reads = sum(reads)) %>% 
        ungroup() %>%
        mutate(sample = as.factor(sample),
               sample.reads = as.numeric(sample.reads),
               sample = fct_reorder(sample, sample.reads, .desc = T)) %>%
        group_by(sample, asv.type) %>%
        mutate(read.counts = sum(reads)) %>% 
        sample_n(1) %>% 
        select(-asv) 
      plot.persample <- df.persample %>%
        ggplot(aes(x = sample, y = read.counts, fill = asv.type)) + 
        geom_col(width = 1) + 
        scale_y_continuous(expand = c(0,0)) + 
        coord_cartesian(ylim = c(0,as.numeric(quantile(df.persample$read.counts, 0.995)))) + 
        labs(title = datasetname,
             subtitle = paste("Threshold =", threshold[j], sep = " ")) + 
        scale_fill_manual(values = c("contamination" = "#F8766D", "sample taxa" = "#00BFC4")) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              plot.title = element_text(hjust = 0),
              plot.subtitle = element_text(size = 14))
      #save plots in a named list
      plots.persample[[j]] <- plot.persample
      names(plots.persample)[j] <- paste("Threshold", threshold[j], sep=": ")
    }
    #generate plot to select optimal threshold value
    decontam.summary.plot <- decontam.summary %>% 
      bind_rows() %>%
      pivot_wider(names_from = "asv.type", values_from = "asv.proportion") %>% 
      group_by(sample.type) %>% 
      mutate(delta.contam = ifelse(threshold == min(threshold),
                                   Contamination,
                                   Contamination - lag(Contamination)),
             threshold = as.factor(round(threshold, 3))) %>% 
      ungroup() %>% 
      ggplot(aes(x = threshold, y = delta.contam, fill = sample.type)) +
      geom_col(position = position_dodge(), color = "black") + 
      scale_y_continuous(labels = scales::percent) +
      labs(y = "Increase in Removed Taxa",
             x = "Threshold Value",
           title = datasetnames[i]) +
      theme(legend.title = element_blank())
    #combine results into nested list
    i.list <- list(decontam.results.list, plots.sampletype, plots.persample, decontam.summary.plot)
    #assign meaningful names to list elements
    names(i.list) <- c("decontam.results", "plots.sampletype", "plots.persample", "decontam.summary.plot")
    #add dataset results to loop output
    decontam.out <- append(decontam.out, list(i.list))
    names(decontam.out)[[i]] <- datasetname #name directory with dataset name 
    
  } else if(grepl("Cryo", datasetfiles[i]) == T){  
      
    #generate ASV abundance table 
    df.asv <- df %>%
      select(contains(c("STP", "NTC")), X) %>%
      column_to_rownames(var = 'X') %>%
      rename_all(~sapply(str_split(string = ., pattern = "_"), `[`, 1)) %>%
      t() %>% 
      as.data.frame() %>% 
      mutate(sample = ifelse(grepl("NTC", rownames(.)), "NTC", "Sample"),
             .before = 1) %>% 
      rownames_to_column(var = "sample.id") 
    
    #load concentration data based on marker data 
    con_values <- 
      if(grepl("16s", datasetname) == T){
        codec %>%
          filter(label.16s %in% df.asv$sample.id) %>% 
          select(label.16s, con_16s)
      } else 
        if(grepl("ITS", datasetname) == T){
          codec %>%
            filter(label.its %in% df.asv$sample.id) %>% 
            select(label.its, con_its)
        }
      
      #Add concentration values to ASV table
      names(con_values) <- c("sample.id", "sample.con")
      df.asv <- df.asv %>% 
        left_join(., con_values, by = "sample.id") %>% 
        relocate(sample.con, .after = sample.id) %>% 
        #annotate wether samples were processed in-house or not 
        mutate(processed = ifelse(grepl("B.", sample.id), "external", "inhouse"),
               .after = sample.con) 
      
      #extract inhouse samples (w/ concentraction values) for decontamination analysis 
      df.decontam <- df.asv %>% 
        filter(processed == "inhouse")
      
      #specify sequence of threshold probability values to iterate through 
      threshold <- if(grepl("16s", datasetname) == T){
        round(seq(0.05,0.95, length.out = 10), 3)
      } else
        if(grepl("ITS", datasetname) == T){
          round(seq(0.35,0.65, length.out = 10), 3)
        }
      
      #generate empty lists for loop output items 
      sapply(c("decontam.results.list", 'decontam.summary', "plots.sampletype", "plots.persample"),
             function(k) assign(k, list(), envir = .GlobalEnv))
      
      #specify dataset being process during loop
      print(datasetname)
      
      #loop cutoff decontamination analysis across specified values  
      for(j in 1:length(threshold)){
        
        #identify contaminating taxa via decontam 
        decontam.results <- df.decontam %>% 
          select(starts_with("ASV_")) %>% 
          as.matrix() %>% 
          isContaminant(seqtab = ., conc = df.decontam$sample.con, method="frequency", threshold = threshold[j]) %>% 
          rownames_to_column("ASV") %>% 
          left_join(., df.tax, by = "ASV") %>% 
          filter(contaminant == T)
        #report threshold value to track loop 
        print(paste("threshold:", threshold[j], sep = ""))
        #save results into list
        decontam.results.list[[j]] <- decontam.results
        #name each dataframe according to the data set and threshold value 
        names(decontam.results.list)[j] <- paste("Threshold", threshold[j], sep=": ")
        
        #extract list of contaminant ASVs 
        contaminants <- filter(decontam.results, contaminant == T)$ASV
        
        #generate long dataframe and classify taxa according to decontam results 
        df.asv.long <- df.asv %>%
          rename(sample.type = sample) %>%
          mutate(sample = rownames(.)) %>% 
          pivot_longer(starts_with("ASV_") | starts_with("Otu"), names_to = "asv", values_to = "reads") %>% 
          mutate(asv.type = ifelse(asv %in% contaminants, "contamination", "sample taxa")) %>% 
          group_by(sample.type) %>% 
          mutate(reads = as.numeric(reads),
                 total.reads = sum(reads))  %>% 
          ungroup()
        
        #plot decontam results according to sample type 
        df.plot.sampletype <- df.asv.long %>% 
          group_by(sample.type, asv.type) %>% 
          mutate(asv.type.reads = sum(reads),
                 asv.proportion = asv.type.reads/total.reads) %>% 
          sample_n(1) %>% 
          mutate(asv.type = as.factor(str_to_title(asv.type)))
        plot.sampletype <- df.plot.sampletype %>% 
          ggplot(aes(x = sample.type, y = asv.proportion, fill = asv.type)) + 
          geom_col() + 
          scale_y_continuous(expand = c(0,0), labels = scales::percent) +
          scale_fill_manual(values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")) +
          labs(title = datasetname,
               subtitle = paste("Threshold =", threshold[j], sep = " ")) + 
          theme(plot.title = element_text(hjust = 0),
                plot.subtitle = element_text(size = 14))
        #extract summary values for decontamination
        decontam.summary[[j]] <- df.plot.sampletype %>% 
          mutate(threshold = threshold[j]) %>% 
          select(sample.type, threshold, asv.type, asv.proportion) %>%
          as.data.frame()
        #save plots in a named list
        plots.sampletype[[j]] <- plot.sampletype
        names(plots.sampletype)[j] <- paste("Threshold", threshold[j], sep=": ")
        #plot decontam results on a per-sample basis 
        df.persample <- df.asv.long %>%
          group_by(sample) %>% 
          mutate(sample.reads = sum(reads)) %>% 
          ungroup() %>%
          mutate(sample = as.factor(sample),
                 sample.reads = as.numeric(sample.reads),
                 sample = fct_reorder(sample, sample.reads, .desc = T)) %>%
          group_by(sample, asv.type) %>%
          mutate(read.counts = sum(reads)) %>% 
          sample_n(1) %>% 
          select(-asv) 
        plot.persample <- df.persample %>%
          ggplot(aes(x = sample, y = read.counts, fill = asv.type)) + 
          geom_col(width = 1) + 
          scale_y_continuous(expand = c(0,0)) + 
          coord_cartesian(ylim = c(0,as.numeric(quantile(df.persample$read.counts, 0.995)))) + 
          labs(title = datasetname,
               subtitle = paste("Threshold =", threshold[j], sep = " ")) + 
          scale_fill_manual(values = c("contamination" = "#F8766D", "sample taxa" = "#00BFC4")) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                plot.title = element_text(hjust = 0),
                plot.subtitle = element_text(size = 14))
        #save plots in a named list
        plots.persample[[j]] <- plot.persample
        names(plots.persample)[j] <- paste("Threshold", threshold[j], sep=": ")
      }
      #generate plot to select optimal threshold value
      decontam.summary.plot <- decontam.summary %>% 
        bind_rows() %>%
        pivot_wider(names_from = "asv.type", values_from = "asv.proportion") %>% 
        group_by(sample.type) %>% 
        mutate(delta.contam = ifelse(threshold == min(threshold),
                                     Contamination,
                                     Contamination - lag(Contamination)),
               threshold = as.factor(round(threshold, 3))) %>% 
        ungroup() %>% 
        ggplot(aes(x = threshold, y = delta.contam, fill = sample.type)) +
        geom_col(position = position_dodge(), color = "black") + 
        scale_y_continuous(labels = scales::percent) +
        labs(y = "Increase in Removed Taxa",
             x = "Threshold Value",
             title = datasetnames[i]) +
        theme(legend.title = element_blank())
      #combine results into nested list
      i.list <- list(decontam.results.list, plots.sampletype, plots.persample, decontam.summary.plot)
      #assign meaningful names to list elements
      names(i.list) <- c("decontam.results", "plots.sampletype", "plots.persample", "decontam.summary.plot")
      #add dataset results to loop output
      decontam.out <- append(decontam.out, list(i.list))
      names(decontam.out)[[i]] <- datasetname #name directory with dataset name
  }
}
      
#keep only final result list (and data set names/paths)
gdata::keep(decontam.out, datasetnames, datasetfiles, date, sure = T)  

#examine visualizations to select optimal value (last value where NTC > Sample)
for(i in 1:length(decontam.out)){
  print(decontam.out[[i]][["decontam.summary.plot"]])
}

#set threshold value for each dataset 
threshold.asg <- data.frame(
  dataset = datasetnames,
  threshold = c(0.45, 0.45, 0.15, 0.483)
)
decontaminated.datasets <- list()

#decontaminate each dataset & export visualizations
for(i in 1:length(decontam.out)){
  df.name <- (datasetnames[i])
  df <- read.csv(datasetfiles[i])
  #generate ASV abundance table 
  df.asv <- df %>%
    select(contains(c("STP", "NTC", "UHM")), X) %>%
    column_to_rownames(var = 'X') %>%
    rename_all(~sapply(str_split(string = ., pattern = "_"), `[`, 1)) %>%
    t() %>% 
    as.data.frame() %>% 
    mutate(sample = ifelse(grepl("NTC", rownames(.)), "NTC", "Sample"),
           .before = 1) %>% 
    rownames_to_column(var = "sample.id") 
  #designation name of dataframe that needs to be pulled from list 
  i.threshold <- threshold.asg %>% 
    filter(dataset == df.name) %>% 
    pull(threshold) %>% 
    paste("Threshold:", ., sep = " ")
  #pull relevant result list 
  i.contaminants <- decontam.out[[i]][["decontam.results"]][[i.threshold]] %>% 
    pull(ASV)
  #remove flagged contaminants 
  df.asv.contam_out <- df.asv %>% 
    select(-any_of(i.contaminants)) %>% 
    select_if(~ !is.numeric(.) || sum(.) > 5) #remove ASVs which were observed 5 times or less 
  #save datasets into list
  decontaminated.datasets[[i]] <- df.asv.contam_out
  names(decontaminated.datasets)[i] <- datasetnames[i]
  #arrange visualizations
  l1 <- decontam.out[[i]][['decontam.summary.plot']]
  l2 <- decontam.out[[i]][["plots.sampletype"]][[i.threshold]]
  l3 <- decontam.out[[i]][["plots.persample"]][[i.threshold]]
  L <- list(l1, l2, l3)
  gridplot <- marrangeGrob(grobs = L, nrow = 2, ncol = 2,
                           layout_matrix = rbind(c(1,1),
                                                 c(2,3)))
  
  #get dataset name
  dataset <- names(decontam.out)[i] %>%
    gsub(" ", "", x = .)
  
  #export plots as single pdf with meaningful name 
  # ggsave(
  #   filename = paste(dataset, "_decontamination.plots_final", ".pdf", sep = ""),
  #   plot = gridplot,
  #   width = 15, height = 9
  # )
}

# conduct rarefaction on each dataset individually and save them to working directory #####
for(i in 1:length(decontaminated.datasets)){
  
  filenames <- c("WF16s", "WFITS", "Cryo16s", "CryoITS")
  dataset <- names(decontaminated.datasets[i])
  dat <- decontaminated.datasets[[i]]
  
  dat.tax <- read.csv(datasetfiles[i]) %>%   
    select(1, contains(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))) %>%
    column_to_rownames(var = 'X') %>% 
    rename_all(~str_to_lower(gsub(x = ., pattern = "tax[.]", replacement = ""))) %>% 
    rownames_to_column("ASV")
  
  dat.analyze <- dat %>%
    select(-sample) %>% 
    column_to_rownames("sample.id") %>% 
    mutate(depth = rowSums(.)) %>%
    rownames_to_column("sample.id") %>%
    select(-starts_with("ASV")) %>%
    arrange(depth)
  rare.cutoff <- quantile(dat.analyze$depth,
                          if(dataset == "Wood Frog - 16s"){
                              .03 }
                          else 
                            if(dataset == "Wood Frog - ITS") {
                              .03 }
                          else 
                            if(dataset == "Cryo Collection - 16s") {
                              .175 }
                          else 
                            if(dataset == "Cryo Collection - ITS") {
                              .175 }
  )
  rare.plot <-   dat.analyze %>% 
    ggplot(aes(x = 1:nrow(.), y = depth)) +
    geom_line() +
    labs(y = "n Reads",
         x = "Samples arranged by Sequencing Depth",
         title = dataset) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_hline(yintercept = rare.cutoff) 
  print(rare.plot)
  
  sample.id.rarified <-
    dat.analyze %>% 
    filter(depth > rare.cutoff) %>% 
    pull(sample.id)
  
  rarified.dat <- filter(dat, sample.id %in% sample.id.rarified) 
  
  asvs.rarified <- 
    rarified.dat %>% 
    select(starts_with("ASV")) %>%
    names()
    
  rarified.tax <- filter(dat.tax, ASV %in% asvs.rarified)
  
  write.csv(rarified.dat, file = paste(date, filenames[i], "rarified.dataset.csv", sep = "_"), row.names = F)
  write.csv(rarified.tax, file = paste(date, filenames[i], "taxonomyfile.csv", sep = "_"), row.names = F)
  
}
