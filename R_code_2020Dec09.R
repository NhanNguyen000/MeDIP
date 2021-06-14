# Using golden path not the Ensembl release (annotation)

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
#folder <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/Testing_2020Dec09/"
file_EPI_The <- 6794:6796
file_Con <- 10:12

samples =data.frame(sample_name=c(paste0(rep("EPI_The_T002_", 3), 1:3), paste0(rep("Con_DMSO_T002_", 3), 1:3)),
                    file_name=c(paste0(folder_EPI, "EPI_L", file_EPI_The, "_pe.sorted.bam"), 
                                paste0(folder_Con, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                    group=c(rep("EPI_The", 3), rep("Control", 3)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=samples, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 1:22),
                      window_size=500)
qseaSet=createQseaSet(sampleTable=samples, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 1),
                      window_size=500)
qseaSet

qseaSet <- addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)  
qseaSet <- addCNV(qseaSet, file_name="file_name",window_size=2e6, 
                  paired=TRUE, parallel=FALSE, MeDIP=TRUE)

qseaSet <- addLibraryFactors(qseaSet) # can add the normalization factors have been pre-computed 
qseaSet <- addPatternDensity(qseaSet, "CG", name="CpG")
qseaSet <- addOffset(qseaSet, enrichmentPattern = "CpG")

wd <- which(getRegions(qseaSet)$CpG_density >1 & getRegions(qseaSet)$CpG_density <15)
signal <- (15-getRegions(qseaSet)$CpG_density[wd]*0.55/15+0.25)
qseaSet_blind <- addEnrichmentParameters(qseaSet, enrichmentPattern="CpG",
                                         windowIdx=wd, signal=signal)


##Quality control
getOffset(qseaSet_blind, scale="fraction")
png("EP_matrix.png")
plotEPmatrix(qseaSet_blind) # returns enrichment profile coordinates for all depicted samples.
dev.off()

##Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
png("test2.png")
plotCNV(qseaSet_blind)
dev.off()

## PCA without specification
pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")

## PCA of  CpG Islands promoter --? 
#data(annotation, package="MEDIPSData") # ?
#CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE) # ? 
#pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom) # ? 
png("pca_test2.png")
col=c(rep("red",3), rep("green", 3))
plotPCA(pca_cgi, bgColor=col)
dev.off()

# Differential methylation analysis
design<-model.matrix(~group, getSampleTable(qseaSet_blind))
qseaGLM<-fitNBglm(qseaSet_blind, design, norm_method="beta")
qseaGLM<-addContrast(qseaSet_blind, qseaGLM, coef=2, name="TvN")

sig <- isSignificant(qseaGLM, fdr_th=0.01) # No region was selected in EPI_The_002 vs EPI_The_008,but have sig region  EPI_Thee_002 vs control


# Annotation:
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #shorthand (for convenience)
txdb

seqlevels(txdb)
columns(txdb)

#transcript
transcript_reg <- transcripts(txdb)
prom_reg <- promoters(txdb, upstream=2000, downstream=2000) #2000 nu up TSS, and 2000 nu down TSS
exon_reg <- exons(txdb)
cds_reg <- cds(txdb)

#genes
genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb, upstream=2000, downstream=2000)
promoters_txdb
exon_reg <- exons(txdb)
cds_reg <- cds(genes_txdb) # problems

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

#BiocManager::install("AnnotationHub")

library("AnnotationHub")

#BiocManager::install("annotatr")
#BiocManager::install("org.Hs.eg.db")
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

# conver gene name
library(biomaRt)

genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("external_gene_name",'entrezgene_id'), values=names(genes),filters ='entrezgene_id', mart = mart)
names(genes) <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]
genes$gene_names <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)] # problem

gene_names <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]
genes_v2 <- cbind(genes, gene_names) # problem
# fix the naem of genome
genome(exon_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(prom_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(transcript_reg) <- "BSgenome.Hsapiens.UCSC.hg38"
genome(cds_reg) <- "BSgenome.Hsapiens.UCSC.hg38"

ROIs_2=list(transcript_reg, exon_reg, prom_reg, cds_reg)
names(ROIs_2)=c("transcript", "exon", "promoter", "coding_region")

#ROIs_2=list(transcript_reg, exon_reg, prom_reg, cds_reg, genes_txdb, promoters_txdb) # fail
#names(ROIs_2)=c("transcript", "exon", "promoter", "coding_region", "gene", "gene_promoter")

result2 <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                     keep=sig, annotation=ROIs_2, norm_method="beta") 

save(qseaSet_blind, qseaGLM, sig, result2, file = "MeDIPresult_full_2020Dec17.RData")
#--> run in loof
#load("D:/TGX/GitHub/MeDIP/MeDIPresult_2020Dec17.RData")
load("D:/TGX/GitHub/MeDIP/MeDIPresult_full_2020Dec17.RData")





