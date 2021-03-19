get.avg_DMR_genes <- function(QSEA_outcome, DMR_cutoff, list_promoter_regions) {
  library(tidyverse)
  DMR_data <- makeGRangesFromDataFrame(QSEA_outcome)
  dm_annotated = annotate_regions( regions = DMR_data,
                                   annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
  
  Select_info <- c("seqnames", "start", "end", "width", "strand",
                   "annot.gene_id", "annot.symbol", "annot.type" )
  
  DMR_gene_annotated <- na.omit(data.frame(dm_annotated))[, Select_info] %>% distinct() %>% # remove the dublicated rows
    group_by(annot.symbol) %>% mutate(count = n()) %>% filter(count >= DMR_cutoff) %>% # remove gene has low DMR than cutoff 
    mutate(promoter = sum(annot.type %in% list_promoter_regions)) %>% filter(promoter > 0) # remove gene have no DMR in promoter resgions
  
  DMR_gene_annotated <- merge(DMR_gene_annotated, The_002[,c("chr", "window_start", "window_end",
                                                             "CpG_density", "TvN_log2FC", "TvN_pvalue", "TvN_adjPval")], 
                              by.x = c("seqnames", "start", "end"), by.y = c("chr", "window_start", "window_end"))
  
  avg_DMR_genes <- DMR_gene_annotated %>% group_by(annot.symbol) %>% summarize(log2FC_avg = mean(TvN_log2FC),
                                                                               pvale_avg = mean(TvN_pvalue),
                                                                               adjPval_avg = mean(TvN_adjPval))
  return(avg_DMR_genes)
}
pre.volcano_plot <- function(avg_DMR_genes, Pvalue, Log2FC_cutoff) {
  avg_DMR_genes$Direction <- "Mild"
  avg_DMR_genes$Direction[avg_DMR_genes$log2FC_avg > Log2FC_cutoff] <- "Up"
  avg_DMR_genes$Direction[avg_DMR_genes$log2FC_avg < -Log2FC_cutoff] <- "Down"
  avg_DMR_genes$label <- ifelse(avg_DMR_genes$Direction == "Mild", NA, avg_DMR_genes$annot.symbol)
  return(avg_DMR_genes)
}
get.volcano_plot<- function(avg_DMR_genes, Pvalue, Log2FC_cutoff) {
  library(ggrepel)
  mycolors <- c("blue", "black", "red")
  names(mycolors) <- c("Down", "Mild", "Up")
  
  ggplot(data = avg_DMR_genes, mapping = aes(x=log2FC_avg, y=-log10(adjPval_avg), col=Direction, label = label)) +
    geom_point() + theme_minimal() + geom_text_repel() +
    geom_vline(xintercept = c(-Log2FC_cutoff, Log2FC_cutoff), col = "red") +
    geom_hline(yintercept = -log10(Pvalue), col = "red") +
    scale_color_manual(values = mycolors)
}
