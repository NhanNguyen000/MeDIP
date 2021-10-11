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

# QSEA pipeline (R code in server) --------------------------------------------

# Futher analysis in Rstudio ---------------------------------------------------
library(GenomicRanges)
library(qsea)
library(tidyverse)
library(annotatr)
source("R_functions.R")
load("./data/ROIs_2_2021Jan19.RData")
load("./data/annotations_2021Sep24.RData")

# Extract and annotate DMRs for all samples --------------------------------------
load("./data/qsea_outcome_EPI_all_samples_20211001.RData")

sig <- isSignificant(qseaGLM, fdr_th=0.01)
QSEA_outcome <- makeTable(qseaSet, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet), 
                          keep=sig, annotation=c(ROIs_2), norm_method="beta")
annotated_genes <-get.annotated_genes(QSEA_outcome, annotations)

# select 5% gene with high DMRs count in gene regions
gene_DMRs <-  annotated_genes[which(annotated_genes$DMR_count>quantile(annotated_genes$DMR_count, .95)),]
ggplot(annotated_genes, aes(x=DMR_count)) + geom_histogram(binwidth = 1, color="gray") + 
  geom_histogram(data=gene_DMRs, aes(x=DMR_count), binwidth=1, color="transparent", fill="red")
dim(annotated_genes)
View(head(annotated_genes))
dim(gene_DMRs)
View(head(gene_DMRs))
# select 5% gene with high DMRs count in promoter regions
genes_DMRs_promoter <- gene_DMRs[which(gene_DMRs$DMR_in_promoter_count>quantile(gene_DMRs$DMR_in_promoter_count, .95)),] 
ggplot(gene_DMRs, aes(x=DMR_in_promoter_count)) + geom_histogram(binwidth = 1) +
  geom_histogram(data=genes_DMRs_promoter, aes(x=DMR_in_promoter_count), binwidth = 1, color="transparent", fill="red")
dim(genes_DMRs_promoter)
View(head(genes_DMRs_promoter))

# calculate average log2FC and p-values per gene & volvano plot
avg_genes_DMRs_protomer <- get.avg_DMRs_genes(genes_DMRs_promoter, QSEA_outcome)
selected_genes <- avg_genes_DMRs_protomer[which(abs(avg_genes_DMRs_protomer$log2FC_avg)>0.5),]
dim(selected_genes)

avg_genes <- get.avg_DMRs_genes(annotated_genes, QSEA_outcome)
volcano_gene <- pre.volcano_plot(avg_genes, selected_genes$annot.symbol, Log2FC_cutoff =0.5)
get.volcano_plot(volcano_gene, Log2FC_cutoff = 0.5)
write.csv(selected_genes, file = "selected_genes_20211006.csv")

# make the pdf file
pdf("EPI_allSamples.pdf", onefile=TRUE)
ggplot(annotated_genes, aes(x=DMR_count)) + geom_histogram(binwidth = 1, color="gray") + 
  geom_histogram(data=gene_DMRs, aes(x=DMR_count), binwidth=1, color="transparent", fill="red")
ggplot(gene_DMRs, aes(x=DMR_in_promoter_count)) + geom_histogram(binwidth = 1) +
  geom_histogram(data=genes_DMRs_promoter, aes(x=DMR_in_promoter_count), binwidth = 1, color="transparent", fill="red")
get.volcano_plot(volcano_gene, Log2FC_cutoff = 0.5)
dev.off()

# clean the terminal
rm(sig, QSEA_outcome, annotated_genes, gene_DMRs, genes_DMRs_promoter, avg_genes_DMRs_protomer, selected_genes, avg_genes, volcano_gene)

# Extract and annotate DMRs per time points --------------------------------------
load("./data/qsea_outcome_EPI_per_timepoint_20211001.RData")
names(outcome)

EPI_time <- list()
for (treat in names(outcome)) {
  sig <- isSignificant(outcome[[treat]]$qseaGLM, fdr_th = 0.01)
  QSEA_outcome <- makeTable(outcome[[treat]]$qseaSet, glm = outcome[[treat]]$qseaGLM,
                            groupMeans=getSampleGroups(outcome[[treat]]$qseaSet),
                            keep=sig, annotation = c(ROIs_2), norm_methods = "beta")
  annotated_genes <- get.annotated_genes(QSEA_outcome, annotations)
  
  gene_DMRs <-  annotated_genes[which(annotated_genes$DMR_count>quantile(annotated_genes$DMR_count, .95)),]
  genes_DMRs_promoter <- gene_DMRs[which(gene_DMRs$DMR_in_promoter_count>quantile(gene_DMRs$DMR_in_promoter_count, .95)),] 
  avg_genes_DMRs_protomer <- get.avg_DMRs_genes(genes_DMRs_promoter, QSEA_outcome)
  selected_genes <- avg_genes_DMRs_protomer[which(abs(avg_genes_DMRs_protomer$log2FC_avg)>0.5),]
  avg_genes <- get.avg_DMRs_genes(annotated_genes, QSEA_outcome)
  volcano_gene <- pre.volcano_plot(avg_genes, selected_genes$annot.symbol, Log2FC_cutoff =0.5)
  
  pdf(paste0("EPI_", treat, ".pdf"), onefile=TRUE)
  ggplot(annotated_genes, aes(x=DMR_count)) + geom_histogram(binwidth = 1, color="gray") + 
    geom_histogram(data=gene_DMRs, aes(x=DMR_count), binwidth=1, color="transparent", fill="red")
  ggplot(gene_DMRs, aes(x=DMR_in_promoter_count)) + geom_histogram(binwidth = 1) +
    geom_histogram(data=genes_DMRs_promoter, aes(x=DMR_in_promoter_count), binwidth = 1, color="transparent", fill="red")
  get.volcano_plot(volcano_gene, Log2FC_cutoff = 0.5)
  dev.off()
  
  EPI_time[[treat]] <- list("sig" = sig, "QSEA_outcome" = QSEA_outcome, "annotated_genes" = annotated_genes,
                            "gene_DMRs" = gene_DMRs, "genes_DMRs_promoter" = genes_DMRs_promoter, "avg_genes_DMRs_protomer" = avg_genes_DMRs_protomer, 
                            "selected_genes" = selected_genes, "avg_genes" = avg_genes, "volcano_gene" = volcano_gene)  
  rm(sig, QSEA_outcome, annotated_genes, gene_DMRs, genes_DMRs_promoter,
     avg_genes_DMRs_protomer, selected_genes, avg_genes, volcano_gene)
}

names(EPI_time)
names(EPI_time$'002')
rm(outcome)
# Extract and annotate DMRs per does and time points --------------------------------------
load("./data/qsea_outcome_EPI_per_dose_timepoint_20211001.RData")
names(outcome)

names(outcome)

EPI_dose_time <- list()
for (treat in names(outcome)) {
  sig <- isSignificant(outcome[[treat]]$qseaGLM, fdr_th = 0.01)
  QSEA_outcome <- makeTable(outcome[[treat]]$qseaSet, glm = outcome[[treat]]$qseaGLM,
                            groupMeans=getSampleGroups(outcome[[treat]]$qseaSet),
                            keep=sig, annotation = c(ROIs_2), norm_methods = "beta")
  annotated_genes <- get.annotated_genes(QSEA_outcome, annotations)
  
  gene_DMRs <-  annotated_genes[which(annotated_genes$DMR_count>quantile(annotated_genes$DMR_count, .95)),]
  genes_DMRs_promoter <- gene_DMRs[which(gene_DMRs$DMR_in_promoter_count>quantile(gene_DMRs$DMR_in_promoter_count, .95)),] 
  avg_genes_DMRs_protomer <- get.avg_DMRs_genes(genes_DMRs_promoter, QSEA_outcome)
  selected_genes <- avg_genes_DMRs_protomer[which(abs(avg_genes_DMRs_protomer$log2FC_avg)>0.5),]
  avg_genes <- get.avg_DMRs_genes(annotated_genes, QSEA_outcome)
  volcano_gene <- pre.volcano_plot(avg_genes, selected_genes$annot.symbol, Log2FC_cutoff =0.5)
  
  pdf(paste0("EPI_", treat, ".pdf"), onefile=TRUE)
  ggplot(annotated_genes, aes(x=DMR_count)) + geom_histogram(binwidth = 1, color="gray") + 
    geom_histogram(data=gene_DMRs, aes(x=DMR_count), binwidth=1, color="transparent", fill="red")
  ggplot(gene_DMRs, aes(x=DMR_in_promoter_count)) + geom_histogram(binwidth = 1) +
    geom_histogram(data=genes_DMRs_promoter, aes(x=DMR_in_promoter_count), binwidth = 1, color="transparent", fill="red")
  get.volcano_plot(volcano_gene, Log2FC_cutoff = 0.5)
  dev.off()
  
  EPI_dose_time[[treat]] <- list("sig" = sig, "QSEA_outcome" = QSEA_outcome, "annotated_genes" = annotated_genes,
                                 "gene_DMRs" = gene_DMRs, "genes_DMRs_promoter" = genes_DMRs_promoter,
                                 "avg_genes_DMRs_protomer" = avg_genes_DMRs_protomer, 
                                 "selected_genes" = selected_genes, "avg_genes" = avg_genes, "volcano_gene" = volcano_gene)
  rm(sig, QSEA_outcome, annotated_genes, gene_DMRs,
     genes_DMRs_promoter, avg_genes_DMRs_protomer, selected_genes, avg_genes, volcano_gene)
}

names(EPI_dose_time$"002_EPI_The")
Time <- substring(names(EPI_dose_time), 1,3)
Dose <- substring(names(EPI_dose_time), 9,11)

DMRs_num <- c()
gene_num <- c()
selected_num <- c()
for(i in names(EPI_dose_time)) {
  DMRs_num <- c(DMRs_num, nrow(EPI_dose_time[[i]]$QSEA_outcome))
  gene_num <- c(gene_num, nrow(EPI_dose_time[[i]]$annotated_genes))
  selected_num <- c(selected_num, nrow(EPI_dose_time[[i]]$selected_genes))
}

dat_plot <- as.data.frame(cbind(Dose, Time, DMRs_num, gene_num, selected_num))
dat_plot

pdf("EPI_dose_timepoint_graphs.pdf", onefile=TRUE)
ggplot(dat_plot, aes(x=Time, y=as.numeric(DMRs_num), colour = Dose)) + 
  geom_point(aes(shape=Dose), size = 5) + geom_line(aes(group = Dose)) + 
  scale_y_continuous(limits=c(0, 15000)) +
  theme(text = element_text(size=15))

ggplot(dat_plot, aes(x=Time, y=as.numeric(gene_num), colour = Dose)) + 
  geom_point(aes(shape=Dose), size = 5) + geom_line(aes(group = Dose)) + 
  scale_y_continuous(limits=c(0, 7000)) +
  theme(text = element_text(size=15))
dev.off()

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


library(annotatr)
annotations = build_annotations(genome = 'hg38', 
                                annotations = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic'))
genome(annotations) <- "BSgenome.Hsapiens.UCSC.hg38"
save(annotations, file = "annotations_2021Sep24.RData")
rm(annotations)

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