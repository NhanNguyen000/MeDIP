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

get.geneExpression <- function(Selected_genes, norm_data, filename) {
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
    if (is.na(Selected_genes$ENSEMBL[i])) {
      print(paste0("no Ensemble ID for this gene: ", Selected_genes$SYMBOL[i]))
    } else {
      Selected_data <- Select_if_in_specific_list(norm_data, Selected_genes$ENSEMBL[i])
      if (length(Selected_data) == 0) {
        print(paste0("no gene expression for this gene: ", 
                     Selected_genes$SYMBOL[i], "-", Selected_genes$ENSEMBL[i]))
      } else {
        dat_tem  <- as.data.frame(cbind(metadata, Selected_data)) %>%
          group_by(Condition) %>%
          summarise_at(vars('Selected_data'), list(mean = mean))
        
        plot_dat <- data.frame(generate_metadata_for_sample(dat_tem$Condition), stringsAsFactors=FALSE) %>% 
          full_join(dat_tem, by = c("Condition"))
        
        plist[[i]] <- ggplot(plot_dat, aes(x = Time, y = mean, colour = Dose, group = Dose)) + 
          geom_line(size = 1) + geom_point(size = 3) + coord_cartesian(ylim =c(0, NA)) +
          xlab("Time (hours)") + ylab("gene expression") +
          ggtitle(Selected_genes$SYMBOL[i]) + theme_bw() +
          theme(plot.title = element_text(size = 16, face = "bold.italic"),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 14),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 16))
      }
    }
  }
  plist <- plist[lengths(plist) != 0]
  ml <- marrangeGrob(grobs=plist, nrow=3, ncol=3)
  pdf(paste0(filename, ".pdf"), width=25,height=12)
  print(ml)
  dev.off()
}

get.geneExpression_multiID <- function(Selected_genes, norm_data, filename) {
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
    if (is.na(Selected_genes$ENSEMBL[i])) {
      print(paste0("no Ensemble ID for this gene: ", Selected_genes$SYMBOL[i]))
    } else {
      Selected_data <- Select_if_in_specific_list(norm_data, Selected_genes$ENSEMBL[i])
      if (length(Selected_data) == 0) {
        print(paste0("no gene expression for this gene: ", 
                     Selected_genes$SYMBOL[i], "-", Selected_genes$ENSEMBL[i]))
      } else {
        dat_tem  <- as.data.frame(cbind(metadata, Selected_data)) %>%
          group_by(Condition) %>%
          summarise_at(vars('Selected_data'), list(mean = mean))
        
        plot_dat <- data.frame(generate_metadata_for_sample(dat_tem$Condition), stringsAsFactors=FALSE) %>% 
          full_join(dat_tem, by = c("Condition"))
        
        plist[[i]] <- ggplot(plot_dat, aes(x = Time, y = mean, colour = Dose, group = Dose)) + 
          geom_line(size = 1) + geom_point(size = 3) + coord_cartesian(ylim =c(0, NA)) +
          xlab("Time (hours)") + ylab("gene expression") +
          ggtitle(paste0(Selected_genes$SYMBOL[i], "-", Selected_genes$ENSEMBL[i])) + theme_bw() +
          theme(plot.title = element_text(size = 16, face = "bold.italic"),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 14),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 16))
      }
    }
  }
  plist <- plist[lengths(plist) != 0]
  ml <- marrangeGrob(grobs=plist, nrow=3, ncol=3)
  pdf(paste0(filename, ".pdf"), width=25,height=12)
  print(ml)
  #print(plist)
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
# overlaped genes between EPI-treated condition
Overlap_genes <- c("MAD1L1", "PRDM15", "NCOR2", "SUN1", 
                    "SPG7", "ANKRD11", "DENND3", "ATP11A")
annot <- AnnotationDbi::select(org.Hs.eg.db, keys = Overlap_genes, 
                               column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                               keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot, norm_data, filename = "GeneExpression_overlap_20220216")

# selected genes in EPI therapeutic-treated condition
hypeMethy_EPI_The <- c("PIGG","PALM","ADAP1","LAMA5",
                    "TSC2","DNM2","MCF2L","MAD1L1",
                    "TCF25","NCOR2","GET4","PRKCZ",
                    "DPP9","KIF1A","SNHG14")
hypeMethy_EPI_The <- hypeMethy_EPI_The[!hypeMethy_EPI_The %in% Overlap_genes]
annot_hypeMethy_EPI_The <- AnnotationDbi::select(org.Hs.eg.db, keys = hypeMethy_EPI_The, 
                                     column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                                     keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot_hypeMethy_EPI_The, norm_data, 
                   filename = "GeneExpression_hypeMethy_EPI_The_20220216")
#hypomethlated genes with multiple Ensemble ID: CTTN, SEPTIN9
hypoMethy_EPI_The <- c("SMARCA4","NPHP4","DNMT1","TNK2","HDAC4",
                       "PRDM15","ANKRD11","SPG7","DENND3","HDLBP",
                       "RNF213","PKN1","ZC3H18","CTTN","NADSYN1","CHFR",
                       "SUN1","CCDC57","RGS12","SEPTIN9","ATP11A","SPTAN1")
hypoMethy_EPI_The <- hypoMethy_EPI_The[!hypoMethy_EPI_The %in% Overlap_genes]
annot_hypoMethy_EPI_The <- AnnotationDbi::select(org.Hs.eg.db, keys = hypoMethy_EPI_The, 
                                                 column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                                                 keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot_hypoMethy_EPI_The, norm_data, 
                   filename = "GeneExpression_hypoMethy_EPI_The_20220216")
# selected genes in EPI toxic-treated condition
#hypomethlated genes with multiple Ensemble ID: POLR2A
hypoMethy_EPI_Tox <- c("EIF3B", "BRD9", "ATP11A", "SUN1", "SDHA", 
                       "NCOR2", "MAD1L1", "PRDM15", "LINC02188", 
                       "CCDC187", "ANKLE2", "AGPAT3", "EHMT1", 
                       "PFKP", "ANKRD11", "DENND3", "POLR2A", "PPP6R2")
hypoMethy_EPI_Tox <- hypoMethy_EPI_Tox[!hypoMethy_EPI_Tox %in% Overlap_genes]
annot_hypoMethy_EPI_Tox <- AnnotationDbi::select(org.Hs.eg.db, keys = hypoMethy_EPI_Tox, 
                                                 column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                                                 keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot_hypoMethy_EPI_Tox, norm_data, 
                   filename = "GeneExpression_hypoMethy_EPI_Tox_20220216")
# Make figure
selected_genes <- c("SMARCA4", "PKN1", "RGS12", "HDAC4",
                     "DPP9", "SDHA", "POLR2A", "AGPAT3")
annot_selected_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = selected_genes, 
                                                 column = c("SYMBOL", "ENSEMBL", "GENENAME"),
                                                 keytype = "SYMBOL", multiVals = "list")
get.geneExpression(annot_selected_genes, norm_data, 
                   filename = "GeneExpression_selected_genes_20220216")

