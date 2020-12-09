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

folder_Con <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/Con_FlucDMSO/"
folder_EPI <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/EPI/"
folder <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/Testing/"
file_EPI <- 6794:6796
file_Con <- 10:12

# using sorted.bam files instead of.bam files --> no warning reduce efficiency, no bam index found --> increase time

samples =data.frame(sample_name=c(paste0(rep("EPI_The_T002_", 3), 1:3), paste0(rep("Con_DMSO_T002_", 3), 1:3)),
                       file_name=c(paste0(folder_EPI, "EPI_The_L", file_EPI, "_pe.sorted.bam"), 
                                   paste0(folder_Con, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                       group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)
# Note: EPI_*_pe.sorted.bamm has problem---------------
Reading bam alignment /ngs-data-2/analysis/NhanNguyen/MeDIP/EPI/EPI_The_L6794_pe.sorted.bam
Error in .io_bam(.scan_bamfile, file, reverseComplement, yieldSize(file),  : 
                   seqlevels(param) not in BAM header:
                   seqlevels: \u2018chr1\u2019, \u2018chr2\u2019, \u2018chr3\u2019, \u2018chr4\u2019, \u2018chr5\u2019, \u2018chr6\u2019, \u2018chr7\u2019, \u2018chr8\u2019, \u2018chr9\u2019, \u2018chr10\u2019, \u2018chr11\u2019, \u2018chr12\u2019, \u2018chr13\u2019, \u2018chr14\u2019, \u2018chr15\u2019, \u2018chr16\u2019, \u2018chr17\u2019, \u2018chr18\u2019, \u2018chr19\u2019, \u2018chr20\u2019, \u2018chr21\u2019, \u2018chr22\u2019
                 file: /ngs-data-2/analysis/NhanNguyen/MeDIP/EPI/EPI_The_L6794_pe.sorted.bam
                 index: /ngs-data-2/analysis/NhanNguyen/MeDIP/EPI/EPI_The_L6794_pe.sorted.bam
                 
                ## ---------------
samples =data.frame(sample_name=c(paste0(rep("EPI_The_T002_", 3), 1:3), paste0(rep("Con_DMSO_T002_", 3), 1:3)),
                    file_name=c(paste0(folder_EPI, "EPI_The_L", file_EPI, "_pe.bam"), 
                                paste0(folder_Con, "Cardiac_FlucDMSO_S", file_Con, "_pe.bam")),
                    group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)
## also probelm in >bam file

                 samples =data.frame(sample_name=c(paste0(rep("EPI_The_T002_", 3), 1:3), paste0(rep("Con_DMSO_T002_", 3), 1:3)),
                                     file_name=c(paste0(folder, "EPI_*_L", file_EPI, "_pe.sorted.bam"), 
                                                 paste0(folder, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                                     group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)
                 
qseaSet=createQseaSet(sampleTable=samples, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 1:22),
                      window_size=500)
qseaSet

qseaSet <- addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)  
qseaSet <- addCNV(qseaSet, file_name="file_name",window_size=2e6, 
                  paired=TRUE, parallel=FALSE, MeDIP=TRUE)

# warnings()
50:   RangedData objects are deprecated. Please migrate your code to use
GRanges or GRangesList objects instead. See IMPORTANT NOTE in
?RangedData

qseaSet <- addLibraryFactors(qseaSet) # can add the normalization factors have been pre-computed 
qseaSet <- addPatternDensity(qseaSet, "CG", name="CpG")
# Warning message:
In estimatePatternDensity(Regions = getRegions(qs), pattern = pattern,  :
                            Masks selected but not found in BSGenome: AGAPS, AMB. 
                          Consider using the .masked version of the package
                 
                                   
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

library(qsea)
library(BSgenome.Hsapiens.UCSC.hg19)
data(samplesNSCLC, package="MEDIPSData")
path=system.file("extdata", package="MEDIPSData")
samples_NSCLC$file_name=paste0(path,"/",samples_NSCLC$file_name )

qseaSet=createQseaSet(sampleTable=samples_NSCLC, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg19", 
                      chr.select=paste0("chr", 20:22), 
                      window_size=500)
qseaSet_0 <- qseaSet

qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)
qseaSet_1 <- qseaSet

qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE)
qseaSet_2 <- qseaSet

qseaSet=addLibraryFactors(qseaSet)
qseaSet_3 <- qseaSet

qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")
qseaSet_4 <- qseaSet

qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")
qseaSet_5 <- qseaSet

wd=which(getRegions(qseaSet)$CpG_density>1 &
           getRegions(qseaSet)$CpG_density<15)
signal=(15-getRegions(qseaSet)$CpG_density[wd])*.55/15+.25
qseaSet_blind=addEnrichmentParameters(qseaSet, enrichmentPattern="CpG",
                                      windowIdx=wd, signal=signal)

getOffset(qseaSet_blind, scale="fraction")
plotEPmatrix(qseaSet)

design=model.matrix(~group, getSampleTable(qseaSet_blind) )

qseaSet_blind2 <- qseaSet_blind
qseaSet_blind2@libraries$file_name[, "library_factor"] <- rep(1, 6)
qseaSet_blind2@libraries$file_name[, "offset"] <- rep(1, 6)
qseaSet_blind2@libraries

qseaGLM=fitNBglm(qseaSet_blind, design, norm_method="beta")
qseaGLM=addContrast(qseaSet_blind,qseaGLM, coef=2, name="TvN" )


qseaGLM2=fitNBglm(qseaSet_blind2, design, norm_method="beta")
qseaGLM2=addContrast(qseaSet_blind2,qseaGLM2, coef=2, name="TvN" )


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

# test


getNormMatrix(qseaSet, methods="beta" ,windows=getRegions(qseaSet),sampleIdx = row.names(design))



qseaSet_1=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)

qseaSet_1
qseaSet_2=addCNV(qseaSet_1, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE)

qseaSet_2

df_qseaSet_2_cnv <- data.frame(qseaSet_2@cnv)
df_qseaSet_1_cnv <- data.frame(qseaSet_1@cnv)

cnv_test <- makeGRangesFromDataFrame(df_qseaSet_2_cnv)
qseaSet_cnv_test=addCNV(qseaSet_1, cnv = cnv_test)

qseaSet_1_1 <- qseaSet_1
qseaSet_1_1@cnv <- cnv_test

qseaSet_3=addLibraryFactors(qseaSet_2)
qseaSet_3_test=addLibraryFactors(qseaSet_1)
qseaSet_cnv_test=addLibraryFactors(qseaSet_cnv_test)

qseaSet_cnv_test=addLibraryFactors(qseaSet_1_1)


qseaSet_4=addPatternDensity(qseaSet_3, "CG", name="CpG")
qseaSet_4_test=addPatternDensity(qseaSet_3_test, "CG", name="CpG")

findCNV(sampleTable=NULL,BSgenome,chr.select=c(20:22),
                  file_name="CNV_file_name",fragment_length=NULL, uniquePos=TRUE, 
                  paired=FALSE, mu =log2(c(1/2, 2/3, 1, 3/2,2,3)), window_size=1000000, 
                  normal_idx=NULL, plot_dir=NULL, MeDIP=FALSE,zygosity=NA,
                  parallel=FALSE)
a<-findCNV(sampleTable=qseaSet@sampleTable) # could not run - no function in qsea?

  
  
## Annotation using Hommer
perl /home/nnguyen/homer/configureHomer.pl -install
echo 'PATH=$PATH:/home/nnguyen/homer/bin/' >>~/.bash_profile # edit your ~/.bash_profile file to include the line: PATH=$PATH:/home/nnguyen/homer/bin/
source ~/.bash_profile
  
perl /home/nnguyen/homer/configureHomer.pl -list
perl /home/nnguyen/homer/configureHomer.pl -install hg38 # v6.4	human genome and annotation for UCSC hg38
annotatePeaks.pl test_2020Now05.txt hg38   > output0_2020Nov05.txt # need to add .gtf file or not? for the custom genes 
annotatePeaks.pl test1_2020Now05.txt hg38   > output1_2020Nov05.txt # need to add .gtf file or not? for the custom genes 


annotatePeaks.pl test_2020Now05.txt hg38 -gtf /ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz  > output0_2020Nov05_gtf.txt

annotatePeaks.pl test_2020Now05.txt hg38 -gtf /ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf  > output0_2020Nov05_gtf_release.txt
annotatePeaks.pl test1_2020Now05.txt hg38 -gtf /ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf  > output1_2020Nov05_gtf_release.txt


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

#-------------old code
library("MEDIPSData")
data("annotation",  package = "MEDIPSData")
result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), keep=sig, annotation=ROIs, norm_method="beta") # not strand --> annotation did not work well
knitr::kable(head(result))
write.table(result, "test_2020Oct09.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


# other code:
a<-result[,c(1:3)]
b<-cbind(a, c(1:nrow(a)))
b<-as.matrix(b)
colnames(b) <- c("Chr", "Start", "End", "PeakID")
write.table(b, "test_2020Oct05.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


# another way of annotation (2020Oct 09)
BiocManager::install("rtracklayer")


library("rtracklayer")
gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
myGRanges <- as(gtfRangeData, "GRanges")
a<- myGRanges$type
head(as.vector(a))
[1] "gene"       "transcript" "exon"       "exon"       "exon"      
[6] "transcript"
> length(as.vector(a))
[1] 2947461
> unique(as.vector(a))
[1] "gene"            "transcript"      "exon"            "CDS"            
[5] "start_codon"     "stop_codon"      "five_prime_utr"  "three_prime_utr"
[9] "Selenocysteine" 



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

# try: ChIPseeker: an R package for ChIP peak Annotation, Comparison and Visualization
https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#peak-annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)


## old code (before 22th Sep 2020) ------------------------------------------------------------
setwd("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/")

sample_ConDF2=data.frame(sample_name=c(paste0(rep("Con_DF2_T000_", 3), 1:3), paste0(rep("Con_DF2_T002_", 3), 1:3)),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/", list.files()),
                         group=c(rep("Normal", 3), rep("Tumor", 3)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:21),
                      window_size=100)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1.bam", "Con_DF2_T002_1.bam"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/", c("Con_DF2_T000_1.bam", "Con_DF2_T002_1.bam")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1", "Con_DF2_T002_1"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Sam_files/", 
                                          c("Con_DF2_T000_1.sorted.bam", "Con_DF2_T002_1.sorted.bam")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1", "Con_DF2_T002_1"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Sam_files/", 
                                          c("Con_DF2_T000_1.sorted.bam.bai", "Con_DF2_T002_1.sorted.bam.bai")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38",
                      window_size=50)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=50)

# Switch for parallel computing, using BiocParallel
library("BiocParallel")
register(MulticoreParam(workers=3))

# Import sequencing data
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=FALSE) # could not run parallel?? + error


qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=FALSE)


# estimate CNV information and add to qseaSet object
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE, normal_idx = "T000")
###Scaling: Estimate effective library size (normalization)
qseaSet=addLibraryFactors(qseaSet)

# Infer sequence pattern density values and add to qseaSet object
#we estimate the average CpG density per fragment for each genomic window.
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")

# Estimate background reads
#From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

##Quality control
getOffset(qseaSet, scale="fraction")
# returns enrichment profile coordinates for all depicted samples.
plotEPmatrix(qseaSet)

##Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
plotCNV(qseaSet)

## PCA of  CpG Islands promoter --? 
data(annotation, package="MEDIPSData") # ?
CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE) # ? 
pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom) # ? 
col=rep(c("red", "green"), 3)

# using sorted.bam files instead of.bam files --> no warning reduce efficiency, no bam index found --> increase time

sample_test=data.frame(sample_name=c(paste0(rep("EPI_The_T002_", 3), 1:3), paste0(rep("Con_DMSO_T002_", 3), 1:3)),
                       file_name=c(paste0(folder, "EPI_*_L", file_EPI, "_pe.sorted.bam"), paste0(folder, "Cardiac_FlucDMSO_S", file_Con, "_pe.sorted.bam")),
                       group=c(rep("EPI", 3), rep("Control", 3)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_test, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=500)
qseaSet=createQseaSet(sampleTable=sample_test, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 1:22),
                      window_size=500)

qseaSet

qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)  
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE)
qseaSet=addLibraryFactors(qseaSet) # can add the normalization factors have been pre-computed 
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

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

# Annotation, 
library(GenomicRanges)
sig <- isSignificant(qseaGLM, fdr_th=0.01) # No region was selected in EPI_The_002 vs EPI_The_008,but have sig region  EPI_Thee_002 vs control

library("MEDIPSData")
data("annotation",  package = "MEDIPSData")
result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), keep=sig, annotation=ROIs, norm_method="beta") # not strand --> annotation did not work well
knitr::kable(head(result))
write.table(result, "test_2020Oct09.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


# other code:
a<-result[,c(1:3)]
b<-cbind(a, c(1:nrow(a)))
b<-as.matrix(b)
colnames(b) <- c("Chr", "Start", "End", "PeakID")
write.table(b, "test_2020Oct05.txt", quote = FALSE, append = FALSE, sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


# another way of annotation (2020Oct 09)
BiocManager::install("rtracklayer")


library("rtracklayer")
gtfRangeData <- import.gff("/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz")
myGRanges <- as(gtfRangeData, "GRanges")
a<- myGRanges$type
head(as.vector(a))
[1] "gene"       "transcript" "exon"       "exon"       "exon"      
[6] "transcript"
> length(as.vector(a))
[1] 2947461
> unique(as.vector(a))
[1] "gene"            "transcript"      "exon"            "CDS"            
[5] "start_codon"     "stop_codon"      "five_prime_utr"  "three_prime_utr"
[9] "Selenocysteine" 



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

# try: ChIPseeker: an R package for ChIP peak Annotation, Comparison and Visualization
https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#peak-annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)


## old code (before 22th Sep 2020) ------------------------------------------------------------
setwd("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/")

sample_ConDF2=data.frame(sample_name=c(paste0(rep("Con_DF2_T000_", 3), 1:3), paste0(rep("Con_DF2_T002_", 3), 1:3)),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/", list.files()),
                         group=c(rep("Normal", 3), rep("Tumor", 3)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:21),
                      window_size=100)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1.bam", "Con_DF2_T002_1.bam"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Bowtie2/", c("Con_DF2_T000_1.bam", "Con_DF2_T002_1.bam")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1", "Con_DF2_T002_1"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Sam_files/", 
                                          c("Con_DF2_T000_1.sorted.bam", "Con_DF2_T002_1.sorted.bam")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

sample_ConDF2=data.frame(sample_name=c("Con_DF2_T000_1", "Con_DF2_T002_1"),
                         file_name=paste0("/ngs-data-2/analysis/NhanNguyen/MeDIP/Sam_files/", 
                                          c("Con_DF2_T000_1.sorted.bam.bai", "Con_DF2_T002_1.sorted.bam.bai")),
                         group=c(rep("T000", 1), rep("T002", 1)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38",
                      window_size=50)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=50)

# Switch for parallel computing, using BiocParallel
library("BiocParallel")
register(MulticoreParam(workers=3))

# Import sequencing data
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=FALSE) # could not run parallel?? + error


qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=FALSE)


# estimate CNV information and add to qseaSet object
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE, normal_idx = "T000")
###Scaling: Estimate effective library size (normalization)
qseaSet=addLibraryFactors(qseaSet)

# Infer sequence pattern density values and add to qseaSet object
#we estimate the average CpG density per fragment for each genomic window.
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")

# Estimate background reads
#From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

##Quality control
getOffset(qseaSet, scale="fraction")
# returns enrichment profile coordinates for all depicted samples.
plotEPmatrix(qseaSet)

##Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
plotCNV(qseaSet)

## PCA of  CpG Islands promoter --? 
data(annotation, package="MEDIPSData") # ?
CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE) # ? 
pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom) # ? 
col=rep(c("red", "green"), 3)



--------------------------------------------------------------- # Old code is here
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("MEDIPS")

# avaiblable genome
library("BSgenome")
available.genomes()



# The reference genome is hg38:
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# Play with toy data -----------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("MEDIPSData")

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library("MEDIPSData")

bam.file.hESCs.Input = system.file("extdata", "hESCs.Input.chr22.bam",
                                   package = "MEDIPSData")
bam.file.DE.Input = system.file("extdata", "DE.Input.chr22.bam",
                                package = "MEDIPSData")

# set parameters
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1e-3
extend=300
shift=0
ws=100
chr.select="chr22"

# Case study: Genome wide methylation and differential coverage between two conditions

# combine replicates
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam",
                                        package = "MEDIPSData")
hESCs_MeDIP = MEDIPS.createSet(file = bam.file.hESCs.Rep1.MeDIP, BSgenome = BSgenome,
                               extend = extend, shift = shift, uniq = uniq,
                               window_size = ws, chr.select = chr.select)

bam.file.hESCs.Rep2.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep2.chr22.bam",
                                        package = "MEDIPSData")
hESCs_MeDIP = c(hESCs_MeDIP, MEDIPS.createSet(file = bam.file.hESCs.Rep2.MeDIP, BSgenome = BSgenome,
                                              extend = extend, shift = shift, uniq = uniq,
                                              window_size = ws, chr.select = chr.select))

hESCs_MeDIP_v2 <- hESCs_MeDIP

# here we load the preprocessed lists of MeDIP-seq MEDIPS SETs available in the MEDIPSData package:
data(hESCs_MeDIP)
data(DE_MeDIP)

hESCs_Input = MEDIPS.createSet(file = bam.file.hESCs.Input, BSgenome = BSgenome,
                               extend = extend, shift = shift, uniq = uniq, window_size = ws,
                               chr.select = chr.select)
DE_Input = MEDIPS.createSet(file = bam.file.DE.Input, BSgenome = BSgenome,
                            extend = extend, shift = shift, uniq = uniq, window_size = ws,
                            chr.select = chr.select)
CS = MEDIPS.couplingVector(pattern = "CG", refObj = hESCs_MeDIP[[1]])


# Coverage, methylation profiles and differential coverage
mr.edgeR = MEDIPS.meth(MSet1 = DE_MeDIP, MSet2 = hESCs_MeDIP,
                       CSet = CS, ISet1 = DE_Input, ISet2 = hESCs_Input, p.adj = "bonferroni",
                       diff.method = "edgeR", MeDIP = T, CNV = F, minRowSum = 10)


# Differential coverage: selecting significant windows
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1,
                              adj = T, ratio = NULL, bg.counts = NULL, CNV = F)

# Merging neighboring significant windows
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC",
                                                     colnames(mr.edgeR.s))] > 0), ]

mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.gain,
                                       distance = 1)

# Extracting data at regions of interest
columns = names(mr.edgeR)[grep("counts", names(mr.edgeR))]
rois = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,
                         columns = columns, summarize = NULL)
rois.s = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,
                           columns = columns, summarize = "avg")

## Quality controls: MEDIPS provides three different quality controls -------------------------------
# run for the MeDIP-seq hESCs sample hESCs_Rep1_MeDIP.bam.

# Saturation analysis
sr = MEDIPS.saturation(file = bam.file.hESCs.Rep1.MeDIP, BSgenome = BSgenome,
                       uniq = uniq, extend = extend, shift = shift, window_size = ws,
                       chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
                       rank = FALSE)
MEDIPS.plotSaturation(sr)

# Correlation between samples
cor.matrix = MEDIPS.correlation(MSets = c(hESCs_MeDIP, DE_MeDIP, hESCs_Input, DE_Input), 
                                plot = T, method = "pearson") # could not run this code

# Sequence Pattern Coverage
cr = MEDIPS.seqCoverage(file = bam.file.hESCs.Rep1.MeDIP, pattern = "CG",
                        BSgenome = BSgenome, chr.select = chr.select, extend = extend,
                        shift = shift, uniq = uniq)
MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", cov.level = c(0,1, 2, 3, 4, 5))
MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="hist", t = 15, main="Sequence pattern coverage, histogram")


# CpG Enrichment
er = MEDIPS.CpGenrich(file = bam.file.hESCs.Rep1.MeDIP, BSgenome = BSgenome,
                      chr.select = chr.select, extend = extend, shift = shift,
                      uniq = uniq) # could not run this code
## Miscellaneous ---------------------
Read: 6. Comments on the experimental design and Input data
# Processing regions of interest
Instead of caluclating coverage and differential coverage at genome wide small windows,
it is also possible to perform targetd analyses of regions of interest (ROI's, e.g. exons, promoter regions, CpG islands etc.).'

# Export Wiggle Files
allows to export genome wide coverage profiles as wiggle files for visualization in common genome browsers.
#MEDIPS.exportWIG(Set = hESCs_MeDIP[[1]], file = "hESC.MeDIP.rep1.wig", format = "rpkm", descr = "")

# Merging MEDIPS SETs
Input.merged = MEDIPS.mergeSets(MSet1 = hESCs_Input, MSet2 = DE_Input,
                                name = "Input.hESCs.DE")

# Annotation of significant windows
anno.mart.gene = MEDIPS.getAnnotation(dataset = c("hsapiens_gene_ensembl"),
                                      annotation = c("GENE"), chr = "chr22")

mr.edgeR.s = MEDIPS.setAnnotation(regions = mr.edgeR.s, annotation = anno.mart.gene)
mr.edgeR.s.gain = MEDIPS.setAnnotation(regions = mr.edgeR.s.gain, annotation = anno.mart.gene)

# addCNV (copy number variation)
mr.edgeR = MEDIPS.addCNV(cnv.Frame = 10000, ISet1 = hESCs_Input,
                         ISet2 = DE_Input, results = mr.edgeR)

# Calibration Plot --? seem to be improtant:
MEDIPS.plotCalibrationPlot(CSet = CS, main = "Calibration Plot",
                           MSet = hESCs_MeDIP[[1]], plot_chr = "chr22", rpkm = TRUE,
                           + xrange = TRUE)
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

#setwd("/ngs-data-2/analysis/NhanNguyen/MeDIP/ConDF2/Alignments/")
#sample_ConDF2=data.frame(sample_name=c(paste0(rep("Con_DF2_T000_", 3), 1:3), paste0(rep("Con_DF2_T002_", 3), 1:3)),
#                        file_name=list.files(),
#                        group=c(rep("T000", 3), rep("T002", 3)), stringsAsFactors=FALSE)

#qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
#                      BSgenome="BSgenome.Hsapiens.UCSC.hg38",
#                      window_size=100)

setwd("D:/")
sample_ConDF2=data.frame(sample_name=c("1T000", "1T002"),
                         file_name=c("D:/Con_DF2_T000_1.bam", "D:/Con_DF2_T002_1.bam"),
                         group=c("T000", "T002"), stringsAsFactors=FALSE)


qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=100)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=100)


sample_ConDF2=data.frame(sample_name=c("T002", "T008"),
                         file_name=c("/ngs-data-2/analysis/NhanNguyen/MeDIP/bwa/IDA_The_002_p2.sorted.bam", 
                                     "/ngs-data-2/analysis/NhanNguyen/MeDIP/bwa/IDA_The_008_p2.sorted.bam"),
                         group=c("T002", "T008"), stringsAsFactors=FALSE)


qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", chr.select=paste0("chr", 20:22),
                      window_size=100)

# Switch for parallel computing, using BiocParallel
library("BiocParallel")
register(MulticoreParam(workers=3))

# Import sequencing data
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=FALSE)

qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=TRUE, MeDIP=TRUE, normal_idx = "T002")

sampleTable=getSampleTable(qseaSet)
Regions=getRegions(qseaSet)
fragment_length=NULL
fname_idx=which(names(sampleTable)=="file_name" )[1]
sam
coverage=unlist(bplapply(X=sampleTable[,fname_idx],FUN=getCoverage, 
                         Regions=Regions,fragment_length=fragment_length,
                         minMapQual=minMapQual, paired=paired, uniquePos=uniquePos, 
                         BPPARAM=BPPARAM ),FALSE, FALSE)

# toy data:
data(samplesNSCLC, package="MEDIPSData")
knitr::kable(samples_NSCLC)
path=system.file("extdata", package="MEDIPSData")
samples_NSCLC$file_name=paste0(path,"/",samples_NSCLC$file_name )
qseaSet=createQseaSet(sampleTable=samples_NSCLC, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38", 
                      chr.select=paste0("chr", 20:22), 
                      window_size=500)
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)
# estimate CNV information and add to qseaSet object
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=TRUE, MeDIP=TRUE, normal_idx = "T000")
###Scaling: Estimate effective library size (normalization)
qseaSet=addLibraryFactors(qseaSet)

# Infer sequence pattern density values and add to qseaSet object
#we estimate the average CpG density per fragment for each genomic window.
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")

# Estimate background reads
#From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

##Quality control
getOffset(qseaSet, scale="fraction")
# returns enrichment profile coordinates for all depicted samples.
plotEPmatrix(qseaSet)

##Exploratory Analysis: Plots a Heatmap-like Overview of the CNVs
plotCNV(qseaSet)

## PCA of  CpG Islands promoter --? 
data(annotation, package="MEDIPSData") # ?
CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE) # ? 
pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom) # ? 
col=rep(c("red", "green"), 3)

# run toy data ---------------------------------------------
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#BiocManager::install("GenomicRanges")
#BiocManager::install("MEDIPSData")
data(samplesNSCLC, package="MEDIPSData")
knitr::kable(samples_NSCLC)

path=system.file("extdata", package="MEDIPSData")
samples_NSCLC$file_name=paste0(path,"/",samples_NSCLC$file_name )

# run qsea
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg19)

qseaSet=createQseaSet(sampleTable=samples_NSCLC, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg19", 
                      chr.select=paste0("chr", 20:22), 
                      window_size=500)
qseaSet

qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE)

qseaSet
##Normalization
#Copy Number Variation (CNV), and account for their influence on MeDIP read density
qseaSet=addCNV(qseaSet, file_name="file_name",window_size=2e6, 
               paired=TRUE, parallel=FALSE, MeDIP=TRUE)
###Scaling
qseaSet=addLibraryFactors(qseaSet)
###Estimating model parameters for transformation to absolute methylation values
#we estimate the average CpG density per fragment for each genomic window.
qseaSet=addPatternDensity(qseaSet, "CG", name="CpG")
#From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

## another exaple data set: 
data(tcga_luad_lusc_450kmeth, package="MEDIPSData")

wd=findOverlaps(tcga_luad_lusc_450kmeth, getRegions(qseaSet), select="first")
signal=as.matrix(mcols(tcga_luad_lusc_450kmeth)[,rep(1:2,3)])
# This function analyses the dependency of enrichment on a sequence pattern, based on a subset of windows for which the signal is known.
qseaSet=addEnrichmentParameters(qseaSet, enrichmentPattern="CpG", 
                                windowIdx=wd, signal=signal)

wd=which(getRegions(qseaSet)$CpG_density>1 &
           getRegions(qseaSet)$CpG_density<15)
signal=(15-getRegions(qseaSet)$CpG_density[wd])*.55/15+.25
qseaSet_blind=addEnrichmentParameters(qseaSet, enrichmentPattern="CpG", 
                                      windowIdx=wd, signal=signal)


##Quality control
getOffset(qseaSet, scale="fraction")
plotEPmatrix(qseaSet)
##Exploratory Analysis
plotCNV(qseaSet)
## PCA of  CpG Islands promoter
data(annotation, package="MEDIPSData")
CGIprom=intersect(ROIs[["CpG Island"]], ROIs[["TSS"]],ignore.strand=TRUE)
pca_cgi=getPCA(qseaSet, norm_method="beta", ROIs=CGIprom)
col=rep(c("red", "green"), 3)

##Differential Methylation Analysis Differential Methylation Analysis in QSEA is based on generalized linear models (GLMs),
design=model.matrix(~group, getSampleTable(qseaSet) )
qseaGLM=fitNBglm(qseaSet, design, norm_method="beta")
qseaGLM=addContrast(qseaSet,qseaGLM, coef=2, name="TvN" )

##Annotating, Exploring and Exporting Results
library(GenomicRanges)
sig=isSignificant(qseaGLM, fdr_th=.01)

result=makeTable(qseaSet, 
                 glm=qseaGLM, 
                 groupMeans=getSampleGroups(qseaSet), 
                 keep=sig, 
                 annotation=ROIs, 
                 norm_method="beta")


plotPCA(pca_cgi, bg=col, main="PCA plot based on CpG Island Promoters")
knitr::kable(head(result))

# To assess the enrichment of differentially methylated regions within genomic annotations
sigList=list(gain=isSignificant(qseaGLM, fdr_th=.1,direction="gain"),
             loss=isSignificant(qseaGLM, fdr_th=.1,direction="loss"))
roi_stats=regionStats(qseaSet, subsets=sigList, ROIs=ROIs)
knitr::kable(roi_stats)
roi_stats_rel=roi_stats[,-1]/roi_stats[,1]
x=barplot(t(roi_stats_rel)*100,ylab="fraction of ROIs[%]",
          names.arg=rep("", length(ROIs)+1),  beside=TRUE, legend=TRUE, 
          las=2, args.legend=list(x="topleft"), 
          main="Feature enrichment Tumor vs Normal DNA methylation")
text(x=x[2,],y=-.15,labels=rownames(roi_stats_rel), xpd=TRUE, srt=30, cex=1, adj=c(1,0))

# If we are interested in a particular genomic region, it can be depicted in a genome browser 
plotCoverage(qseaSet, samples=getSampleNames(qseaSet), 
             chr="chr20", start=38076001, end=38090000, norm_method="beta", 
             col=rep(c("red", "green"), 3), yoffset=1,space=.1, reorder="clust", 
             regions=ROIs["TFBS"],regions_offset=.5, cex=.7 ) 
#Parallelization
# A large part of the run time is required for processing the alignment files. 
# These steps can be parallelized using the BiocParallel package:
  
library("BiocParallel")
register(MulticoreParam(workers=3))
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=TRUE)

# annotation for the real data:
load("D:/TGX/GitHub/MeDIP/myGRanges")

# another way of annotation:
library(BSgenome.Hsapiens.UCSC.hg38)
test = BSgenome.Hsapiens.UCSC.hg38
gene_body = genes(test)
