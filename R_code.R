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
    
    samples =data.frame(sample_name=c(paste0("EPI_L", file_EPI), paste0("ConDMSO_S", file_Con)),
                        file_name=c(paste0(folder_EPI, "EPI_L", file_EPI, "_pe.sorted.bam"), 
                                    paste0(folder_Con, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                        group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)
    
    qsea.get_DMR(samples, output_name = paste0(EPI_Dose, time))
  }
}

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


# Annotation, 
library(GenomicRanges)
sig <- isSignificant(qseaGLM, fdr_th=0.01) # No region was selected in EPI_The_002 vs EPI_The_008,but have sig region  EPI_Thee_002 vs control

library("rtracklayer")
gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
myGRanges <- as(gtfRangeData, "GRanges")

result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=list(the_ranges=myGRanges), norm_method="beta") # MAKE LIST 

save(qseaSet_blind, qseaGLM, sig, file = "midterm_outcome_2020Nov16.RData")

load("midterm_outcome_2020Nov16.RData")
load("myGRanges.RData")
data(annotation, package="MEDIPSData") --> ROIs
# test with new ROIs:
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #shorthand (for convenience)
txdb

#genome(txdb) <- "BSgenome.Hsapiens.UCSC.hg38" - could not run

seqlevels(txdb)
columns(txdb)

#keytypes(txdb)
#keys <- c("100033416", "100033417", "100033420")
#select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")

#transcript
transcript_reg <- transcripts(txdb)
prom_reg <- promoters(txdb, upstream=2000, downstream=2000) #2000 nu up TSS, and 2000 nu down TSS
exon_reg <- exons(txdb)
cds_reg <- cds(txdb)

#genes
#genes_txdb <- genes(txdb)
1613 genes were dropped because they have exons located on both strands of the same
reference sequence or on more than one reference sequence, so cannot be represented by a
single genomic range., but using  'single.strand.genes.only=FALSE' -- could nto get normal format outcome

genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb, upstream=2000, downstream=2000)
promoters_txdb
exon_reg <- exons(txdb)
cds_reg <- cds(genes_txdb)

# covragane range

exon_gene <- exonsBy(txdb, by="gene")
genes_cvg <- coverageByTranscript(genes_txdb, exon_gene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txBygene<- transcriptsBy(txdb, by="gene")
cds_gene <- cdsBy(txdb, "gene")


threeUTR <- threeUTRsByTranscript(txdb)
fiveUTR <- fiveUTRsByTranscript(txdb)

cds_gene2 <- cds(txdb, "gene")


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

# conver gene name
library(biomaRt)

genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("external_gene_name",'entrezgene_id'), values=names(genes),filters ='entrezgene_id', mart = mart)
names(genes) <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]
genes$gene_names <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]

# fix the naem of genome
genome(exon_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(prom_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(transcript_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(cds_reg) <- "BSgenome.Hsapiens.UCSC.hg38"

ROIs_2=list(transcript_reg, exon_reg, prom_reg, cds_reg)
names(ROIs_2)=c("transcript", "exon", "promoter", "coding_region")

result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=ROIs_2, norm_method="beta") # mistake

# error could not read file

result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=list(exon=exon_reg), norm_method="beta") # MAKE LIST 

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
## normalization

library(GenomicRanges)
sig=isSignificant(qseaGLM, fdr_th=.01)
data(annotation, package="MEDIPSData")
result=makeTable(qseaSet_blind, 
                 glm=qseaGLM, 
                 groupMeans=getSampleGroups(qseaSet), 
                 keep=sig, 
                 annotation=ROIs, 
                 norm_method="beta")
knitr::kable(head(result))


result2=makeTable(qseaSet_blind2, 
                 glm=qseaGLM2, 
                 groupMeans=getSampleGroups(qseaSet), 
                 keep=sig, 
                 annotation=ROIs, 
                 norm_method="beta")
knitr::kable(head(result2))

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
