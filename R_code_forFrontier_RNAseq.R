#functions
get.RNA_data_combine <- function(list_dataset, data_type) {
  
  RNA_data_combine  <- merge(list_dataset[[1]][[data_type]], list_dataset[[2]][[data_type]], by = "row.names")
  RNA_data_combine  <- get.rownames(RNA_data_combine)
  
  if (length(list_dataset) > 2) {
    for (i in c(3:length(list_dataset))) {
      RNA_data_combine  <- merge(RNA_data_combine, list_dataset[[i]][[data_type]], by = "row.names")
      RNA_data_combine  <- get.rownames(RNA_data_combine)
    }
  }
  
  return(RNA_data_combine)
}



get.rownames <- function(data) {
  rownames(data) <- data[, 1]
  data  <- data[, -1]
  return(data)
}


get.norm_data <- function(RNA_data_count) {
  metadata <- generate_metadata_for_sample(colnames(RNA_data_count))
  library("DESeq2")
  dds <- DESeqDataSetFromMatrix(countData = round(RNA_data_count), 
                                colData = metadata, design = ~ Condition)
  dds2 <- estimateSizeFactors(dds)
  norm_data <- counts(dds2,normalized=TRUE)
  return(norm_data)
}

generate_metadata_for_sample <- function(list_data) {
  output <- matrix(data=NA, ncol=4, nrow= length(list_data))
  colnames(output) <- c("Compound", "Dose", "Time", "Condition")
  rownames(output) <- list_data
  output[, "Compound"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 1)
  output[, "Dose"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 2)
  output[, "Time"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 3)
  output[, "Condition"] <- paste0(output[, "Compound"], "_", output[, "Dose"], "_", output[, "Time"])
  return(output)
}


Select_Part_Names <- function(list, separator, selected_position) {
  # eg. for separator = "|" --> "[|]"
  output<-c()
  for (i in 1: length(list)) output[i] <- strsplit(list[i], paste0("[", separator, "]"))[[1]][selected_position]
  return(output)
}


get.log <- function(data, log_base) {
  data[data == 0]  <- 1
  output           <- log(data, log_base)
  return(output)
}


Select_if_in_specific_list <- function (data, selected_list) {
  output <- data[rownames(data) %in% selected_list,]
  return(output)
}

get.geneExpression <- function(Selected_genes, norm_data) {
  library("tidyverse")			
  library("gridExtra")
  library("magrittr")
  library("ggpubr")
  library("plyr")	
  
  metadata <- data.frame(generate_metadata_for_sample(colnames(norm_data)), stringsAsFactors=FALSE)
  metadata$Dose <- paste0(metadata$Compound, "_", metadata$Dose)
  
  # plot
  plist <- list()
  for(i in 1:nrow(Selected_genes)) {
    Selected_data <- Select_if_in_specific_list(norm_data, Selected_genes$ENSEMBL[i])
    
    Log_invitro_tem  <- as.data.frame(cbind(metadata, Selected_data)) %>%
      group_by(Condition) %>%
      summarise_at(vars('Selected_data'), list(mean = mean))
    
    plot_dat <- data.frame(generate_metadata_for_sample(Log_invitro_tem$Condition), stringsAsFactors=FALSE) %>% 
      full_join(Log_invitro_tem, by = c("Condition"))
    
    plist[[i]] <- ggplot(plot_dat, aes(x = Time, y = mean, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + coord_cartesian(ylim =c(0, NA)) +
      xlab("Time (hours)") + ylab("gene expression") +
      ggtitle(Selected_genes$SYMBOL[i]) +  theme_bw()
  }
  ml <- marrangeGrob(plist, nrow=3, ncol=3)
  pdf("GeneExpression.pdf", width=25,height=12)
  print(ml)
  dev.off()
}
#code -------------------------------------
load("./data/EPIdata_Cardiac_NN_20191112.RData")
load("./data/Con_DF2data_Cardiac_NN_20191112.RData")

RNA_data_count  <- get.RNA_data_combine(list(Con_DF2, EPI), "expected_count")
RNA_data_count  <- RNA_data_count[, -grep("Con_DF2_000_", colnames(RNA_data_count))]
norm_data       <- get.norm_data(RNA_data_count)

# selected genes
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
selected_genes <- c("ANKRD11", "DENND3", "EHMT1", "HDAC4",
                    "MAD1L1", "NCOR2", "PFKP", "SDHA", "TCF25")
annot <- AnnotationDbi::select(org.Hs.eg.db, keys = selected_genes, 
                               column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                               keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot, norm_data)
