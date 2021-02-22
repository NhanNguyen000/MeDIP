## QSEA code -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("qsea")
#BiocManager::install("BSgenome")
library("BSgenome")
available.genomes()

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg38)

## QSEA code (run in serve) -------------------------------

## Load file for analysis ------------------------------------------------
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

get.QC_plots <- function(qseaSet_blind) {
  print("Enrichment profile")
  print(getOffset(qseaSet_blind, scale="fraction"))
  plotEPmatrix(qseaSet_blind) # enrichment matrix
  plotCNV(qseaSet_blind) # a Heatmap-like Overview of the CNVs
  pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
  plotPCA(pca_cgi)
}

library(qsea)
qsea_outcome <- list.files()[grepl("qsea_outcome_" ,list.files())]
#QC_plots <- get.QC_plots(qseaSet_blind)
pdf("test.pdf", onefile = T)
for(qsea_result in qsea_outcome) {
  load(qsea_result)
  print(qsea_result)
  get.QC_plots(qseaSet_blind)
}
dev.off()

load("qsea_outcome_EPI_The_allSamples.RData" )
load("qsea_outcome_EPI_Tox_allSamples.RData" )
pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
time_series <- c("2", "8", "24", "72", "168", "240", "336")
pca_cgi@sample_names <- rep((rep(time_series, each=3)), 2)
plotPCA(pca_cgi, bg=rep(c("red", "green"), each=21))

# Annotation --------------------------------------------------------
# Annotation - option 1:---------------------------------------------
#library(GenomicRanges)
#sig <- isSignificant(qseaGLM, fdr_th=0.01) # No region was selected in EPI_The_002 vs EPI_The_008,but have sig region  EPI_Thee_002 vs control
#library("rtracklayer")
#gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
#myGRanges <- as(gtfRangeData, "GRanges")
#result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
#                     keep=sig, annotation=list(the_ranges=myGRanges), norm_method="beta") # MAKE LIST 


# Annotation - option 2.1: using annotatr package: -------------------------
#explaining the annotation: https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
library(annotatr)
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

# Annotation - option 2.2: using txdb: -------------------------
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


## Get the DMRs with annotation -----------------------------------------------
#Check these papers:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0241515
#https://academic.oup.com/eep/article/3/3/dvx016/4098081?login=true


load("ROIs_2021Jan19.RData")
load("ROIs_2_2021Jan19.RData")
library(GenomicRanges)

load("qsea_outcome_EPI_The002.RData")
sig <- isSignificant(qseaGLM, fdr_th=0.01)
The_002 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")

load("qsea_outcome_EPI_The008.RData")
sig <- isSignificant(qseaGLM, fdr_th=0.01)
The_008 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                    keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")

load("qsea_outcome_EPI_The024.RData")
sig <- isSignificant(qseaGLM, fdr_th=0.01)
The_024 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")

time_series <- c("002", "008", "024", "072", "168", "240", "336")
EPI_The <- list()
for(time in time_series) {
  load(paste0("qsea_outcome_EPI_The", time, ".RData"))
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  output_tem <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                       keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")
  EPI_The[[time]] <- dplyr::select(output_tem, c(chr, window_start, window_end))
  
}

# using flow diagram:
a1<- dplyr::select(The_002, c(chr, window_start, window_end)) #2393 hit
a2<- dplyr::select(The_008, c(chr, window_start, window_end)) # 1931 hit
a3<- dplyr::select(The_024, c(chr, window_start, window_end)) # 14590 hit
k1<-inner_join(a1, a2) # 476 hit
k2<-inner_join(a2, a3) #841 hit
k12 <- inner_join(k1, k2) # 374 hit

nodes = data.frame("name" = paste0("The_", names(EPI_The)))
links<-matrix(NA, nrow = 1, ncol=3)
for (i in 1:length(EPI_The)) {
  for (j in (i+1):length(EPI_The)) {
    outcome_tem <- c((i-1), (j-1), nrow(inner_join(EPI_The[[i]], EPI_The[[j]])))
    links <- rbind(links, outcome_tem)
  }
}
links <- as.data.frame(links[-1,])
names(links) = c("source", "target", "value")
links$group <- as.factor(ifelse(links$source==0, "002", 
                                ifelse(links$source == 1, "008", 
                                       ifelse(links$source==2, "024", 
                                              ifelse(links$source==3, "072",
                                                     ifelse(links$source==4, "168",
                                                            ifelse(links$source==5, "240", "336")))))))

nodes$group <- as.factor(substring(nodes$name, 5, 9))
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              LinkGroup = 'group', NodeGroup = "group",
              fontSize= 12, nodeWidth = 30)

# check the gene region:
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

qsea_outcome_region <- list()
qsea_outcome <- list.files()[grep("qsea_outcome",list.files())]
qsea_outcome <- qsea_outcome[c(2:8, 10:16)]
for (i in 1:length(qsea_outcome)) {
  load(qsea_outcome[i])
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                      keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")
  qsea_outcome_region[[qsea_outcome[i]]] <- c(get.sum_region(result), "total_hit" = nrow(result))
}

DMR_allEPI_con <- bind_rows(qsea_outcome_region)
DMR_allEPI_con$condition <- substring(names(qsea_outcome_region), 14, 23)
save(qsea_outcome_region,DMR_allEPI_con ,  file = "qsea_outcome_region_2021Feb18.RData")
load("qsea_outcome_region_2021Feb18.RData")

data_test<-gather(DMR_allEPI_con, key = "test", value = "value", colnames(DMR_allEPI_con)[1:11])
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

#result_v2<-get.ROIs_lg(result, regions)
#knitr::kable(head(result))
#which(result$id=="") # --> annoation using "annotatr" package provide annotation for all region 

#load("qsea_sig_annot_2021Jan28.RData") # no file?


#gsub("[[:blank:]]","", unlist(strsplit(result_v2$id[1], ",")))
library(tidyverse)
#ggplot(data=result_v2) + geom_bar(mapping = aes(x=genes_promoter))
#sum(result_v2$genes_promoter, na.rm = T)/nrow(result_v2)

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

get.gene_type <- function(input, gene_info) {
  gene_type <- rep(NA, length(input))
  for (i in 1:length(input)) {
    gene_type_tem <- c()
    
    if (length(grep(",", input[i])) ==0) {
      gene_type_tem <- gene_info$type_of_gene[grep(input[i], gene_info$Symbol)]
      if (length(gene_type_tem) >0) gene_type[i] <- gene_type_tem
      
    } else {
      gene_symbol <- strsplit(input[i], "[,]")[[1]]
      gene_type_tem <- c()
      for (j in 1:length(gene_symbol)) {
        gene_type_tem <- c(gene_type_tem, gene_info$type_of_gene[grep(gsub("[[:blank:]]", "", gene_symbol[j]), gene_info$Symbol)])
      }
      gene_type[i] <- paste(gene_type_tem, collapse = ", ")
    }
  }
  return(gene_type)
}

res_gene_symbol <- get.gene_symbol(result)
identical(res_gene_symbol$ENTREZID, result$gene_region)
res_gene_symbol <- cbind(result, res_gene_symbol$SYMBOL)

# do the volcano plot
tmp<-res_gene_symbol[,c("TvN_log2FC", "TvN_adjPval", "res_gene_symbol$SYMBOL")] # annoation from the txbd package
#tmp<-res_gene_symbol[,c("TvN_log2FC", "TvN_adjPval","symbol")] # nnoation from the anotatr package
colnames(tmp) <- c("log2FC", "adjPval", "gene_symbol") 
#which(tmp$adjPval>0.05) # no genes
# calculate mean for duplicated gene symbols
library(tidyverse)
avg_tmp<- tmp %>% group_by(gene_symbol) %>% summarize(log2FC_mean = mean(log2FC),
                                                      adjPval_mean = mean(adjPval))


gene_info<- read.delim("Homo_sapiens.gene_info")
#unique(gene_info$type_of_gene)
gene_type <-get.gene_type(avg_tmp$gene_symbol, gene_info) 
avg_tmp <- cbind(avg_tmp, gene_type)

filted_tmp <- avg_tmp[complete.cases(avg_tmp),] # remove the row has NA in gene_symble and gene_type



Pvalue <- 0.01
Log2FC_cutoff <- 1.5
filted_tmp$diffexpressed <- "MILD"
filted_tmp$diffexpressed[filted_tmp$log2FC_mean > Log2FC_cutoff & filted_tmp$adjPval_mean <Pvalue] <- "UP"
filted_tmp$diffexpressed[filted_tmp$log2FC_mean < -Log2FC_cutoff & filted_tmp$adjPval_mean <Pvalue] <- "DOWN"
#filted_tmp$label[filted_tmp$diffexpressed != "MILD"] <- filted_tmp$gene_symbol[filted_tmp$diffexpressed != "MILD"]

index <- grep("ncRNA", filted_tmp$gene_type)
filted_tmp$label2 <- NA
filted_tmp$label2[index] <- filted_tmp$gene_symbol[index]
a<-filted_tmp[complete.cases(filted_tmp),]

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
