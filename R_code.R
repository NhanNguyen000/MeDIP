## QSEA code -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("qsea")
#BiocManager::install("BSgenome")
library("BSgenome")
available.genomes()

# try to run Con_DF data:
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg38)

## QSEA code -------------------------------

folder_Con <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/ConFlucDMSO/"
folder_EPI <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/EPI/"
#file_EPI <- 6794:6796
#file_Con <- 10:12

qsea.get_DMR <- function(samples, output_name) {
  qseaSet=createQseaSet(sampleTable=samples, 
                        BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 1:22),
                        window_size=500)
  qseaSet
  
  qseaSet <- addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)  
  qseaSet <- addCNV(qseaSet, file_name="file_name",window_size=2e6, 
                    paired=TRUE, parallel=FALSE, MeDIP=TRUE)
  qseaSet <- addLibraryFactors(qseaSet) 
  qseaSet <- addPatternDensity(qseaSet, "CG", name="CpG") # have a warning message about the masks selected: AGAPS, AMB
  qseaSet <- addOffset(qseaSet, enrichmentPattern = "CpG")
  
  wd <- which(getRegions(qseaSet)$CpG_density >1 & getRegions(qseaSet)$CpG_density <15)
  signal <- (15-getRegions(qseaSet)$CpG_density[wd]*0.55/15+0.25)
  qseaSet_blind <- addEnrichmentParameters(qseaSet, enrichmentPattern="CpG",
                                           windowIdx=wd, signal=signal)
  
  # Differential methylation analysis
  design<-model.matrix(~group, getSampleTable(qseaSet_blind))
  qseaGLM<-fitNBglm(qseaSet_blind, design, norm_method="beta")
  qseaGLM<-addContrast(qseaSet_blind, qseaGLM, coef=2, name="TvN")
  
  save(qseaSet_blind, qseaGLM, file= paste0("qsea_outcome_", output_name, ".RData"))
}


data_number <- as.data.frame(matrix(NA, nrow = 7, ncol = 3))
rownames(data_number) <-c("002", "008", "024", "072", "168", "240", "336")
colnames(data_number) <- c("Con", "EPI_The", "EPI_Tox")

data_number$Con <- seq(10, 30, 3)
data_number$EPI_The <- seq(6794, 6814, 3)
data_number$EPI_Tox <- seq(6815, 6835, 3)

for(time in row.names(data_number)) {
  file_Con <- data_number[time, "Con"] : (data_number[time, "Con"]+2)
  
  for(EPI_Dose in c("EPI_The", "EPI_Tox")) {
    file_EPI <- data_number[time, EPI_Dose] : (data_number[time, EPI_Dose]+2)
    
    samples<-data.frame(sample_name=c(paste0("EPI_L", file_EPI), paste0("ConDMSO_S", file_Con)),
                        file_name=c(paste0(folder_EPI, "EPI_L", file_EPI, "_pe.sorted.bam"), 
                                    paste0(folder_Con, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                        group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)
    
    qsea.get_DMR(samples, output_name = paste0(EPI_Dose, time))
  }
}

ConDMSO <- seq(10, 30)
EPI_The <- seq(6794, 6814)
EPI_Dox <- seq(6815, 6835)

The_samples <- data.frame(sample_name = c(paste0("EPI_The_L", EPI_The), paste0("ConDMSO_S", ConDMSO)),
                         file_name = c(paste0(folder_EPI, "EPI_L", EPI_The, "_pe.sorted.bam"), 
                                       paste0(folder_Con, "Cardiac_FlucDMSO_S", ConDMSO, "_pe.sorted.bam")),
                         group = c(rep("EPI_The", length(EPI_The)), rep("Control", length(ConDMSO))), stringsAsFactors = FALSE)
qsea.get_DMR(The_samples, output_name = "EPI_The_allSamples")

Tox_samples <- data.frame(sample_name = c(paste0("RPI_Tox_L", RPI_Tox), paste0("ConDMSO_S", ConDMSO)),
                          file_name = c(paste0(folder_EPI, "EPI_L", EPI_Tox, "_pe.sorted.bam"), 
                                        paste0(folder_Con, "Cardiac_FlucDMSO_S", ConDMSO, "_Pe.sorted.bam")),
                          group = c(rep("EPI_Tox", length(EPI_Tox)), rep("Control", length(ConDMSO))), stringsAsFactors = FALSE)
qsea.get_DMR(Tox_samples, output_name = "EPI_Tox_allSamples")

## Load file for analysis ------------------------------------------------
# Quality control: enrichment profile
getOffset(qseaSet_blind, scale="fraction")
png("EPi_matrix.png")
plotEPmatrix(qseaSet_blind)
dev.off()

#Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
png("test2.png")
plotCNV(qseaSet_blind)
dev.off()

## PCA of samples without specification
pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
png("pca_test2.png")
col<- rep(c("red", "green"), 3)
plotPCA(pca_cgi, bgColor=col)
dev.off()


# Annotation - option 1:
library(GenomicRanges)
sig <- isSignificant(qseaGLM, fdr_th=0.01) # No region was selected in EPI_The_002 vs EPI_The_008,but have sig region  EPI_Thee_002 vs control

library("rtracklayer")
gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
myGRanges <- as(gtfRangeData, "GRanges")

result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=list(the_ranges=myGRanges), norm_method="beta") # MAKE LIST 


# Annotation - option 2:

explaining the annotation: https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
library(annotatr)
annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)
genome(annotations) <- "BSgenome.Hsapiens.UCSC.hg38"

id<- annotations[,1]
tx_id <- annotations[,2]
gene_id <- annotations[,3]
symbol <- annotations[,4]
type <- annotations[,5]
ROIs_v2 <- list(id, tx_id, gene_id, symbol, type)
names(ROIs_v2) <- c("id", "tx_id", "gene_id", "symbol", "type")
#a<- str_sub(unique(annotations$type), 6) # problem with "s"at the end
a<- c("genes_promoter", "genes_1to5kb", "genes_5UTR", "genes_exon", 
      "genes_intron,", "genes_3UTR", "genes_intergenic",
      "cpg_island", "cpg_shore", "cpg_shelve", "cpg_inter")

## using txdb
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb
methods(class=class(txdb))

seqlevels(txdb)
columns(txdb)
keytypes(txdb)

#transcript
genetic_features <- c("tx_id", "tx_name", "gene_id")
transcript_reg <- transcripts(txdb, columns = genetic_features)
prom_reg <- promoters(txdb, upstream=2000, downstream=2000, columns = genetic_features) #2000 nu up TSS, and 2000 nu down TSS
exon_reg <- exons(txdb, columns = c("exon_id", genetic_features))
cds_reg <- cds(txdb, columns = c("cds_id", genetic_features))
gene_reg <- genes(txdb)

rm(txdb, genetic_features)

#gene_reg <- genes(txdb, single.strand.genes.only=FALSE)
#gene_reg2 <- genes(txdb, columns = genetic_features, single.strand.genes.only=FALSE)
#threeUTR <- threeUTRsByTranscript(txdb)
#fiveUTR <- fiveUTRsByTranscript(txdb)

#fix the genome name
genome(transcript_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(prom_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(exon_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(cds_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(gene_reg) <- "BSgenome.Hsapiens.UCSC.hg38"

ROIs <- list(transcript_reg, prom_reg, exon_reg, cds_reg, gene_reg)
names(ROIs) <- c("transcript", "promoter", "exon", "coding_region", "gene_region")
save(ROIs, file = "ROIs_2021Jan19.RData")
rm(ROIs, transcript_reg, prom_reg, exon_reg, cds_reg, gene_reg)


## Get the DMRs -----------------------------------------------
Check these papers:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0241515
https://academic.oup.com/eep/article/3/3/dvx016/4098081?login=true

load("qsea_outcome_EPI_The002.RData")
load("ROIS_2021Jan19.RData")
#library(GenomicRanges)
sig <- isSignificant(qseaGLM, fdr_th=0.01)
result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                    keep=sig, annotation=ROIs, norm_method="beta")
knitr::kable(head(result))

result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=c(ROIs, ROIs_v2), norm_method="beta")
result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                    keep=sig, annotation=ROIs_v2, norm_method="beta")
which(result$id=="") # --> annoation using "annotatr" package provide annotation for all region 

k <- matrix(NA, ncol = length(a), nrow=nrow(result))
colnames(k) <- a
result_v2 <- cbind(result, k)

for (i in 1: nrow(result_v2)) {
  for (region in a) {
    if (length(grep(unlist(strsplit(region, "_"))[2], result_v2$id[i])) >0) result_v2[i, region] <- 1
  }
}


gsub("[[:blank:]]","", unlist(strsplit(result_v2$id[1], ",")))

# conver gene name
get.gene_symbol <- function(res) {
  library(org.Hs.eg.db)
  GeneSymbols <- select(org.Hs.eg.db,
                        keys = res$gene_region,
                        columns = c("SYMBOL","ENTREZID"),
                        keytype = "ENTREZID")
  multi_gene_ids <- grep(",", res$gene_region)
  if (length(multi_gene_ids) >0) {
    for (i in multi_gene_ids) {
      gene_ids <- strsplit(res$gene_region[i], "[,]")[[1]]
      convert_table <- select(org.Hs.eg.db, keys = gsub("[[:blank:]]", "", gene_ids),
                              columns = c("SYMBOL", "ENTREZID"), keytype = "ENTREZID")
      GeneSymbols[i, "SYMBOL"] <- paste(convert_table$SYMBOL, collapse = ", ")
    }
  }
  return(GeneSymbols)
}

res_gene_symbol <- get.gene_symbol(result)
identical(res_gene_symbol$ENTREZID, result$gene_region)
res_gene_symbol <- cbind(result, res_gene_symbol$SYMBOL)

# do the volcano plot
tmp<-res_gene_symbol[,c("TvN_log2FC", "TvN_adjPval", "res_gene_symbol$SYMBOL")]
colnames(tmp) <- c("log2FC", "adjPval", "gene_symbol")
# which(tmp$adjPval>0.05) # no genes
tmp <- tmp[complete.cases(tmp),]
# calculate mean for duplicated gene symbols
library(tidyverse)
filted_tmp<- tmp %>% group_by(gene_symbol) %>% summarize(log2FC_mean = mean(log2FC), 
                                                adjPval_mean = mean(adjPval))

Pvalue <- 0.001
Log2FC_cutoff <- 1.5
filted_tmp$diffexpressed <- "MILD"
filted_tmp$diffexpressed[filted_tmp$log2FC_mean > Log2FC_cutoff & filted_tmp$adjPval_mean <Pvalue] <- "UP"
filted_tmp$diffexpressed[filted_tmp$log2FC_mean < -Log2FC_cutoff & filted_tmp$adjPval_mean <Pvalue] <- "DOWN"
filted_tmp$label[filted_tmp$diffexpressed != "MILD"] <- filted_tmp$gene_symbol[filted_tmp$diffexpressed != "MILD"]


library(ggrepel)
mycolors <- c("blue", "black", "red")
names(mycolors) <- c("DOWN", "MILD", "UP")

ggplot(data = filted_tmp, mapping = aes(x=log2FC_mean, y=-log10(adjPval_mean), col=diffexpressed)) +
  geom_point() + theme_minimal() + 
  geom_vline(xintercept = c(-Log2FC_cutoff, Log2FC_cutoff), col = "red") +
  geom_hline(yintercept = -log10(Pvalue), col = "red") +
  scale_color_manual(values = mycolors)


# with gene_symbol
ggplot(data = filted_tmp, mapping = aes(x=log2FC_mean, y=-log10(adjPval_mean), col=diffexpressed, label = label)) +
  geom_point() + theme_minimal() + geom_text_repel() +
  geom_vline(xintercept = c(-Log2FC_cutoff, Log2FC_cutoff), col = "red") +
  geom_hline(yintercept = -log10(Pvalue), col = "red") +
  scale_color_manual(values = mycolors)

# check the lncRNA -

# covragane range
exon_gene <- exonsBy(txdb, by="gene")
genes_cvg <- coverageByTranscript(genes_txdb, exon_gene)
txBygene<- transcriptsBy(txdb, by="gene")
cds_gene <- cdsBy(txdb, "gene")


# CpG island
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")

library("AnnotationHub")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("annotatr")
BiocManager::install("org.Hs.eg.db")
library("annotatr")
library("org.Hs.eg.db")
annots = c('hg38_cpgs', 'hg38_basicgenes')
annotations = build_annotations(genome = 'hg38', annotations = annots)
dm_annotated = annotate_regions(regions = genes_txdb, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
dim(df_dm_annotated)
View(df_dm_annotated)
unique(df_dm_annotated$annot.type) 
# --> promoter, 1to5kb, 5UTR, exon, intron, 3utrs cpg island, cpg-shore, cpg-shelve, cpg inter

df_genes_txdb = data.frame(genes_txdb)
dim(df_genes_txdb)

# TFBS
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("liftOver")

library(rtracklayer)
library(liftOver)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

seqlevelsStyle(genes_txdb) = "UCSC"  # necessary
cur19 = liftOver(genes_txdb, ch)
class(cur19)

cur19 = unlist(cur19)
genome(cur19) = "hg19"
cur19 = new("gwaswloc", cur19)
cur19 # enhancer from hg19 are lifted to hg38

length(genes_txdb)-length(cur19) # loss/gain locus
setdiff(mcols(cur)$SNPS, mcols(cur19)$SNPS)

## Compute the transcript coverage with coverageByTranscript():
tx_cvg <- coverageByTranscript(txdb, transcripts)
tx_cvg


# In .Seqinfo.mergexy(x, y) :
The 2 combined objects have no sequence levels in common. (Use  suppressWarnings() to suppress this warning.)
knitr::kable(head(result2))

write.table(result2, "test_2020Now05.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
result_homer <- cbind(c(1:nrow(result2)), result2[,c(1:3)], rep(0, nrow(result2)))
colnames(result_homer) <- c("PeakID", "Chr", "Start", "End", "Strand")
write.table(result_homer, "test_2020Now05.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

result_homer1 <- cbind(c(1:nrow(result2)), result2[,c(1:3)], rep(1, nrow(result2)))
colnames(result_homer1) <- c("PeakID", "Chr", "Start", "End", "Strand")
write.table(result_homer1, "test1_2020Now05.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)



# stop here
#setwd("D:/TGX/GitHub/MeDIP")
output0 <-read.table("output0_2020Nov05.txt", fill = TRUE, sep = "\t", header=TRUE)
output1 <-read.table("output1_2020Nov05.txt", fill = TRUE, sep = "\t", header=TRUE)

$ check if the strand +/- have the same result
output0 <-read.table("output0_2020Nov05_order.txt", fill = TRUE, sep = "\t", header=TRUE)
output1 <-read.table("output1_2020Nov05_order.txt", fill = TRUE, sep = "\t", header=TRUE)

ncol(output0)
for(i in 1:19) {
  print(identical(output0[,i], output1[,i]))
  
}
--> not all the same #There is some different between adding +/-
length(unique(output0$Entrez.ID))
length(unique(output1$Entrez.ID))
length(unique(output1$Entrez.ID, output1$Entrez.ID))

# using gtf customer file did not help
output0_gtf <-read.table("output0_2020Nov05_gtf.txt", fill = TRUE, sep = "\t", header=TRUE) # --> no annotation
output0_gtf_release <-read.table("output0_2020Nov05_gtf_release.txt", fill = TRUE, sep = "\t", header=TRUE) # --> no annotation??

#look at data
unique(output0$Gene.Type)
unique(output1$Gene.Type)
length(which(output0$Gene.Type == "protein-coding"))

get.gene_type_count <- function(data, col.name) {
  outcome <- matrix(NA, nrow = length(unique(data[[col.name]])), ncol = 1)
  rownames(outcome) <- unique(data[[col.name]])
  for (i in 1:length(unique(data[[col.name]]))) {
    outcome[i,] = length(which(data[[col.name]] == unique(data[[col.name]])[i]))
  }
  return(outcome)
}
get.gene_type_count(output0, "Gene.Type")
get.gene_type_count(output1, "Gene.Type")

a0<-output0[which(output0$Gene.Type == "ncRNA"), ]
a1<-output1[which(output1$Gene.Type == "ncRNA"), ]


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsembl()
View(listAttributes(ensembl))
k0 <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'description', 
                     'start_position', 'end_position', 'strand', 
                     'external_gene_name', 'gene_biotype'),
      filters = 'ensembl_gene_id', 
      values = a0$Nearest.Ensembl, mart = ensembl)
k1 <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'description', 
                           'start_position', 'end_position', 'strand', 
                           'external_gene_name', 'gene_biotype'),
            filters = 'ensembl_gene_id', 
            values = a1$Nearest.Ensembl, mart = ensembl)
get.gene_type_count(k0, 'gene_biotype')
get.gene_type_count(k1, 'gene_biotype')

length(unique(c(k0$ensembl_gene_id, k1$ensembl_gene_id))) # 337 unique

g0<-k0[which(k0$gene_biotype == "lncRNA"), ]
g1<-k1[which(k1$gene_biotype == "lncRNA"), ] # sinilar lncRNAs

write.table(t(g0$ensembl_gene_id), "test.txt", sep = ", ",
            quote = F, row.names = F, col.names = F)
write.table(t(a0$Entrez.ID), "test_Entrez.txt", sep = ", ",
            quote = F, row.names = F, col.names = F)
a<-read.csv("LncTarD_output0_2020Nov06.csv")
unique(a$Regulator)
b<-read.csv("LncTarD_output0_2020Nov06_Entrez.csv")
unique(b$Regulator)
# another way of annotation (2020Oct 09)
BiocManager::install("rtracklayer")


library("rtracklayer")
gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
myGRanges <- as(gtfRangeData, "GRanges")
a<- myGRanges$type
head(as.vector(a))
length(as.vector(a))
unique(as.vector(a))

result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), keep=sig, annotation=myGRanges, norm_method="beta") # could not use, require a name list of Grange object
result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), keep=sig, annotation=list(gene_body=myGRanges), norm_method="beta") # MAKE LIST 
knitr::kable(head(result2))

# annotation:
library(GenomicRanges)
gr1 <- GRanges(seqnames=Rle(c("ch1", "chMT"), c(2, 4)), ranges=IRanges(16:21, 20), strand=rep(c("+", "-", "*"), 2)) # but nothing knew about the CpG island, intron, exon,...

perl /home/nnguyen/homer/configureHomer.pl -list
perl /home/nnguyen/homer/configureHomer.pl -install hg38 # v6.4	human genome and annotation for UCSC hg38
annotatePeaks.pl test_2020Oct05.txt hg38   > outputfile.txt # need to add .gtf file or not? 


sigList=list(gain=isSignificant(qseaGLM, fdr_th=.1,direction="gain"),
             loss=isSignificant(qseaGLM, fdr_th=.1,direction="loss"))
roi_stats=regionStats(qseaSet, subsets=sigList, ROIs=ROIs) # could not run , no ROIs
knitr::kable(roi_stats)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
# annotation for the real data:
load("D:/TGX/GitHub/MeDIP/myGRanges")

# another way of annotation:
library(BSgenome.Hsapiens.UCSC.hg38)
test = BSgenome.Hsapiens.UCSC.hg38
gene_body = genes(test)
