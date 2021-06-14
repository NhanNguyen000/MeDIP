## Load packages and function -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("qsea")
#BiocManager::install("BSgenome")
library("BSgenome")
#available.genomes()

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg38)

source("R_functions.R")

# Make annotation --------------------------------------------------------
# Annotation - option 2.1: using annotatr package
#explaining the annotation: https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
#BiocManager::install("annotatr")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(annotatr)
library("org.Hs.eg.db")
annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)

genome(annotations) <- "BSgenome.Hsapiens.UCSC.hg38"

id<- annotations[,1]
tx_id <- annotations[,2]
gene_id <- annotations[,3]
symbol <- annotations[,4]
type <- annotations[,5]
ROIs_2 <- list(id, tx_id, gene_id, symbol, type)
names(ROIs_2) <- c("id", "tx_id", "gene_id", "symbol", "type")

regions <- c("genes_promoter", "genes_1to5kb", "genes_5UTR", "genes_exon", 
             "genes_intron", "genes_3UTR", "genes_intergenic",
             "cpg_island", "cpg_shore", "cpg_shelve", "cpg_inter")
save(ROIs_2, regions, file = "ROIs_2_2021Jan19.RData")
rm(annots, annotations, id, tx_id, gene_id, symbol, type, ROIs_2, regions)

## QSEA code (run in serve) -------------------------------
## Load file for analysis
# Quality control: enrichment profile
#getOffset(qseaSet_blind, scale="fraction")
#png("EPi_matrix.png")
#plotEPmatrix(qseaSet_blind)
#dev.off()

#Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
#png("test2.png")
#plotCNV(qseaSet_blind)
#dev.off()

## PCA of samples without specification
#pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
#png("pca_test2.png")
#col<- rep(c("red", "green"), 3)
#plotPCA(pca_cgi, bgColor=col)
#plotPCA(pca_cgi)
#dev.off()
## run code for Rif MeDIP data ------------------------------------------------
library(qsea)
qsea_outcome <- list.files("./data/")[grepl("Rif" ,list.files("./data/"))]
#QC_plots <- get.QC_plots(qseaSet_blind)
pdf("test.pdf", onefile = T)
for(qsea_result in qsea_outcome) {
  load(qsea_result)
  print(qsea_result)
  get.QC_plots(qseaSet_blind)
}
dev.off()

load("./data/qsea_outcome_Rif_The_allSamples.RData")
load("./data/qsea_outcome_Rif_Tox_allSamples.RData")
pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
time_series <- c("2", "8", "24", "72", "168", "240", "336")
pca_cgi@sample_names <- rep((rep(time_series, each=3)), 2)
plotPCA(pca_cgi, bg=rep(c("red", "green"), each=21))

# using annotatr package: -------------------------
#BiocManager::install("annotatr")
#BiocManager::install("org.Hs.eg.db")
library(annotatr)
library(org.Hs.eg.db)

annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
genome(annotations) <- "BSgenome.Hsapiens.UCSC.hg38"

DMR_cutoff <- 4 # genes have at leat this cutoff DMR (methylation detections) are selected
list_promoter_regions <- c("hg38_genes_1to5kb", "hg38_genes_promoters")
Pvalue <- 0.01; Log2FC_cutoff <- 1.5

load("./data/ROIs_2_2021Jan19.RData")
condition <- c("The072", "The168", "Tox072", "Tox168")
avg_DMR_genes <- list()
plist <- list()
for(i in condition) {
  load(paste0("./data/qsea_outcome_Rif_", i, ".RData"))
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  QSEA_outcome <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                            keep=sig, annotation=c(ROIs_2), norm_method="beta")
  avg_DMR_genes[[i]] <- get.avg_DMR_genes(QSEA_outcome, DMR_cutoff, list_promoter_regions)
  plist[[i]] <- get.volcano_plot(pre.volcano_plot(avg_DMR_genes[[i]], Pvalue, Log2FC_cutoff),Pvalue, Log2FC_cutoff)
}
pdf("Rif_vocalno.pdf", onefile = T)
print(plist)
dev.off()
save(avg_DMR_genes, file = "Rif_avg_DMR_gene.RData")

## run code for EPI MeDIP data ------------------------------------------------
library(qsea)
qsea_outcome <- list.files("./data/")[grepl("EPI" ,list.files("./data/"))]
#QC_plots <- get.QC_plots(qseaSet_blind)
pdf("test.pdf", onefile = T)
for(qsea_result in qsea_outcome) {
  load(qsea_result)
  print(qsea_result)
  get.QC_plots(qseaSet_blind)
}
dev.off()

load("./data/qsea_outcome_EPI_The_allSamples.RData")
load("./data/qsea_outcome_EPI_Tox_allSamples.RData")
pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
time_series <- c("2", "8", "24", "72", "168", "240", "336")
pca_cgi@sample_names <- rep((rep(time_series, each=3)), 2)
plotPCA(pca_cgi, bg=rep(c("red", "green"), each=21))

## Get the DMRs of EPI with annotation -----------------------------------------------
#Check these papers:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0241515
#https://academic.oup.com/eep/article/3/3/dvx016/4098081?login=true

load("./data/ROIs_2_2021Jan19.RData")
library(GenomicRanges)

load("./data/qsea_outcome_EPI_The002.RData")
sig <- isSignificant(qseaGLM, fdr_th=0.01)
The_002 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=c(ROIs_2), norm_method="beta")

# using annotatr package & make vocanol plot - DMRs of EPI: -------------------------
#BiocManager::install("annotatr")
#BiocManager::install("org.Hs.eg.db")
library(annotatr)
library(org.Hs.eg.db)

annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
genome(annotations) <- "BSgenome.Hsapiens.UCSC.hg38"

DMR_cutoff <- 4 # genes have at leat this cutoff DMR (methylation detections) are selected
list_promoter_regions <- c("hg38_genes_1to5kb", "hg38_genes_promoters")
Pvalue <- 0.01; Log2FC_cutoff <- 1.5

load("./data/ROIs_2_2021Jan19.RData")
time <- c("002", "008", "024", "072", "168", "240", "336")
avg_DMR_genes_The <- list()
plist <- list()
for(i in time) {
  load(paste0("./data/qsea_outcome_EPI_The", i, ".RData"))
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  QSEA_outcome <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                       keep=sig, annotation=c(ROIs_2), norm_method="beta")
  avg_DMR_genes_The[[i]] <- get.avg_DMR_genes(QSEA_outcome, DMR_cutoff, list_promoter_regions)
  plist[[i]] <- get.volcano_plot(pre.volcano_plot(avg_DMR_genes_The[[i]], Pvalue, Log2FC_cutoff),Pvalue, Log2FC_cutoff)
}
pdf("The_vocalno.pdf", onefile = T)
print(plist)
dev.off()

avg_DMR_genes_Tox <- list()
plist <- list()
for(i in time) {
  load(paste0("./data/qsea_outcome_EPI_Tox", i, ".RData"))
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  QSEA_outcome <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                            keep=sig, annotation=c(ROIs_2), norm_method="beta")
  avg_DMR_genes_Tox[[i]] <- get.avg_DMR_genes(QSEA_outcome, DMR_cutoff, list_promoter_regions)
  plist[[i]] <- get.volcano_plot(pre.volcano_plot(avg_DMR_genes_Tox[[i]], Pvalue, Log2FC_cutoff),Pvalue, Log2FC_cutoff)
}
pdf("Tox_vocalno.pdf", onefile = T)
print(plist)
dev.off()

## general analysis ----------------------------------------------
## explain how to pick the cut off of 3 or 4 DMRs:----------------------------------
library(tidyverse)
DMR_data <- makeGRangesFromDataFrame(QSEA_outcome)
dm_annotated = annotate_regions(regions = DMR_data,
                                 annotations = annotations, ignore.strand = TRUE, quiet = FALSE)

Select_info <- c("seqnames", "start", "end", "width", "strand",
                 "annot.gene_id", "annot.symbol", "annot.type" )


DMR_gene_annotated <- na.omit(data.frame(dm_annotated))[, Select_info] %>% distinct() %>% # remove the dublicated rows
  group_by(annot.symbol) %>% mutate(count = n()) %>% # count the DMRs per gene
  select(c("annot.symbol", "count")) %>% distinct()
  
ggplot(DMR_gene_annotated) + geom_bar(mapping = aes(x=count, fill = count)) + 
  scale_x_continuous(name="Number of DMRs", breaks=seq(0, 10, 1)) +
  geom_vline(xintercept = 3.5, colour = "red")

# pathway analysis:----------------------------------------------------------
library(pathfindR)
time <- c("002", "008", "024", "072", "168", "240", "336")
output_The <- list()
clustered_The <- list()
pdf("pathway_The.pdf", onefile = T)
for (i in time){
  input <- as.data.frame(avg_DMR_genes_The[[i]][,c(1:3)])
  output_The[[i]] <- run_pathfindR(input)
  clustered_The[[i]] <- cluster_enriched_terms(output_df)
  term_gene_heatmap(result_df = output_The[[i]], genes_df = input)
}
dev.off()

write.table(avg_DMR_genes$annot.symbol, "test.txt", col.names = F, row.names = F)

# need to check the MeDIP in cpg island?
# check with the gene expression in RNAseq: MeDIP log2FC up/down --> gene expression log2 up/down
# from: https://www.researchgate.net/post/How_to_find_cpg_islands_in_promoter_region_of_given_gene
# get the gene with cpg island

## get the gene information ------------------------------------------------
library(org.Hs.eg.db)
Ensemble_database <-read.csv("D:/TGX/GitHub/lncRNA_EPI/data/Ensemble_mart_export_NN_20190815.txt")

Gene_Enseml_Ids <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(DMR_gene_annotated$annot.symbol), 
                                         columns = c("SYMBOL",  "ENSEMBL"), keytype = "SYMBOL")
Gene_Enseml_Ids$SYMBOL[which(is.na(Gene_Enseml_Ids$ENSEMBL))]
Gene_info <- Ensemble_database[which(Ensemble_database$Gene.stable.ID %in% Gene_Enseml_Ids$ENSEMBL), ]

lncRNAs <- unique(Gene_info$Gene.stable.ID[Gene_info$Gene.type == "lncRNA"])
protein_coding_genes <- unique(length(Gene_info$Gene.stable.ID[which(Gene_info$Gene.type=="protein_coding")]))


# check the gene region: ------------------------------------------------------------
qsea_outcome_region <- list()
qsea_outcome <- list.files("./data/")[grep("EPI",list.files("./data/"))]
qsea_outcome <- qsea_outcome[c(2:8, 10:16)]
for (i in qsea_outcome) {
  load(paste0("./data/", i))
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                      keep=sig, annotation=c(ROIs_2), norm_method="beta")
  qsea_outcome_region[[i]] <- c(get.sum_region(result), "total_hit" = nrow(result))
}

DMR_all_conditions <- bind_rows(qsea_outcome_region)
DMR_all_conditions$condition <- substring(names(qsea_outcome_region), 14, 23)
#save(qsea_outcome_region,DMR_all_conditions ,  file = "qsea_EPIoutcome_region_2021Feb18.RData")
#load("./data/qsea_EPIoutcome_region_2021Feb18.RData")

data_test<-gather(DMR_all_conditions, key = "test", value = "value", colnames(DMR_all_conditions)[1:11])
data_test$dose <- substring(data_test$condition, 5, 7)
data_test$time <- substring(data_test$condition, 8, 11)

library(tidyverse)
library(viridis)
ggplot(data = data_test, mapping = aes(x=time, y=value, fill = test)) + 
  geom_bar(stat = "identity", position = "dodge") +facet_grid(dose~.) +
  scale_color_viridis(discrete=TRUE)

library(dplyr)
data_test2 <- data_test %>% distinct(condition, dose, time, total_hit)
ggplot(data = data_test2, mapping = aes(x=time, y=total_hit)) + 
  geom_bar(stat= "identity") +facet_grid(dose~.) +
  scale_color_viridis(discrete=TRUE)
