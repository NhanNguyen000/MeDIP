get.annotated_genes <- function(QSEA_outcome, annotations) {
  library(tidyverse)
  
  DMR_data <- makeGRangesFromDataFrame(QSEA_outcome)
  dm_annotated <- annotate_regions(regions = DMR_data,
                                   annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
  Select_info <- c("seqnames", "start", "end", "width", "strand",
                   "annot.gene_id", "annot.symbol", "annot.type")
  list_promoter_regions <- c("hg38_genes_1to5kb", "hg38_genes_promoters")
  
  annotated_genes <- na.omit(data.frame(dm_annotated))[, Select_info] %>% 
    filter( is.na(annot.symbol) == FALSE) %>% distinct() %>% # remove unannotated regions (to gene)) and the dublicated rows
    group_by(annot.symbol) %>% mutate(DMR_count = n()) %>% # count the DMRs per gene
    mutate(DMR_in_promoter_count = sum(annot.type %in% list_promoter_regions)) %>% # count the DMRs in the promoter region
    distinct(annot.symbol, .keep_all=T) 
  return(annotated_genes)
}

get.avg_DMRs_genes <- function(selected_genes, QSEA_outcome) {
  library(tidyverse)
  
  DMR_gene_annotated <- merge(selected_genes, QSEA_outcome[,c("chr", "window_start", "window_end",
                                                                  "CpG_density", "TvN_log2FC", "TvN_pvalue", "TvN_adjPval")], 
                              by.x = c("seqnames", "start", "end"), by.y = c("chr", "window_start", "window_end"))
  
  avg_DMR_genes <- DMR_gene_annotated %>% group_by(annot.symbol) %>% summarize(log2FC_avg = mean(TvN_log2FC),
                                                                               pvale_avg = mean(TvN_pvalue),
                                                                               adjPval_avg = mean(TvN_adjPval))
  return(avg_DMR_genes)
}

pre.volcano_plot <- function(dat, selected_gene, Log2FC_cutoff) {
  dat$Direction <- "Non_selected"
  dat$Direction[dat$annot.symbol %in% selected_gene] <- "selected"
  dat$Direction[dat$Direction == "selected" & dat$log2FC_avg > Log2FC_cutoff] <- "Up"
  dat$Direction[dat$Direction == "selected" &dat$log2FC_avg < -Log2FC_cutoff] <- "Down"
  dat$label <- ifelse(dat$Direction == "Non_selected", NA, dat$annot.symbol)
  return(dat)
}

get.volcano_plot<- function(dat, Pvalue =0.01, Log2FC_cutoff) {
  library(ggrepel)
  ggplot(data = dat, aes(x=log2FC_avg, y=-log10(adjPval_avg), fill= Direction)) +
    geom_point(data = dat[dat$Direction == "Non_selected", ],
               colour="grey", alpha = 0.5, size=1, shape =16) +
    geom_point(data = dat[dat$Direction == "Down", ], 
               size=2, shape = 21, fill = "steelblue", colour = "black") +
    geom_point(data = dat[dat$Direction == "Up",], 
               size=2, shape = 21, fill = "firebrick", colour = "black") +
    geom_text_repel(data = dat[!is.na(dat$label),], aes(label = label),
                    size = 3, max.overlaps = Inf) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10))) + 
    geom_vline(xintercept = c(-Log2FC_cutoff, Log2FC_cutoff), linetype="dashed") +
    geom_hline(yintercept = -log10(Pvalue), linetype="dashed")
}
## check if need these code ----------------------------------

get.sum_region <- function(result) {
  gene_region <-strsplit(paste(result$id, collapse = ", "), "[,]")[[1]]
  regions <- matrix(unlist(strsplit(gene_region, "[:]")), ncol=2, byrow = T)[,1]
  regions <- sort(unique(gsub("[[:blank:]]", "", regions)))
  
  output <- matrix(NA, ncol = length(regions), nrow = nrow(result))
  colnames(output) <- regions
  for(i in 1:nrow(result)) {
    for(j in 1:ncol(output)) {
      res_tem <- grep(colnames(output)[j], result$id[i])
      if (length(res_tem) ==0) output[i,j] <- 0 else output[i,j] <- res_tem
    }
  }
  return(colSums(output))
}

get.QC_plots <- function(qseaSet_blind) {
  print("Enrichment profile")
  print(getOffset(qseaSet_blind, scale="fraction"))
  plotEPmatrix(qseaSet_blind) # enrichment matrix
  plotCNV(qseaSet_blind) # a Heatmap-like Overview of the CNVs
  pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
  plotPCA(pca_cgi)
}