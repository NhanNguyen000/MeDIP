suppressPackageStartupMessages(library('GenomicFeatures'))

gffFile <- "/ngs-data/analysis/hecatos/NhanNguyen/Genome/bwa_genome_CRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz"
txdb <- makeTxDbFromGFF(file=gffFile, format= "gtf")
txdb2 <- makeTxDbFromGFF(file=gffFile)

#-------------------------------

# Genic annotation

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #shorthand (for convenience)
txdb

seqlevels(txdb)
columns(txdb)

#keytypes(txdb)
#keys <- c("100033416", "100033417", "100033420")
#select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")

transcript_reg <- transcripts(txdb)
prom_reg <- promoters(txdb, upstream=2000, downstream=400)
exon_reg <- exons(txdb)
cds_reg <- cds(txdb)


ROIs=list(transcript_reg, exon_reg, prom_reg, cds_reg)
names(ROIs)=c("transcript", "exon", "promoter", "coding_region")


## CpG annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")

library(AnnotationHub)
ah = AnnotationHub()
ah
unique(ah$dataprovider)
unique(ah$species)
unique(ah$rdataclass)

dm <- query(ah, c("GRanges", "UCSC", "Homo sapiens" )) # chose which data provider?

epiFiles <- query(ah, "EpigenomeRoadMap")
unique(epiFiles$species)
unique(epiFiles$genome)


