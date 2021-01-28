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
EPI_Tox <- seq(6815, 6835)

The_samples <- data.frame(sample_name = c(paste0("EPI_The_L", EPI_The), paste0("ConDMSO_S", ConDMSO)),
                          file_name = c(paste0(folder_EPI, "EPI_L", EPI_The, "_pe.sorted.bam"), 
                                        paste0(folder_Con, "Cardiac_FlucDMSO_S", ConDMSO, "_pe.sorted.bam")),
                          group = c(rep("EPI_The", length(EPI_The)), rep("Control", length(ConDMSO))), stringsAsFactors = FALSE)
qsea.get_DMR(The_samples, output_name = "EPI_The_allSamples")

Tox_samples <- data.frame(sample_name = c(paste0("EPI_Tox_L", EPI_Tox), paste0("ConDMSO_S", ConDMSO)),
                          file_name = c(paste0(folder_EPI, "EPI_L", EPI_Tox, "_pe.sorted.bam"), 
                                        paste0(folder_Con, "Cardiac_FlucDMSO_S", ConDMSO, "_pe.sorted.bam")),
                          group = c(rep("EPI_Tox", length(EPI_Tox)), rep("Control", length(ConDMSO))), stringsAsFactors = FALSE)
qsea.get_DMR(Tox_samples, output_name = "EPI_Tox_allSamples")

## Annotation: ---------------------------------------------------
get.ROIs_lg <- function(input, regions){
  region_temp <- matrix(NA, ncol = length(regions), nrow=nrow(input))
  colnames(region_temp) <- regions
  output <- cbind(input, region_temp)
  
  for (i in 1: nrow(output)) {
    for (region in regions) {
      if (length(grep(unlist(strsplit(region, "_"))[2], output$id[i])) >0) output[i, region] <- 1
    }
  }
  return(output)
}


load("ROIs_2021Jan19.RData")
load("ROIs_2_2021Jan19.RData")
library(GenomicRanges)
library(qsea)


qsea_outcome <- list.files()[grepl("qsea_outcome_" ,list.files())]
qsea_sig_annot <- list()
for(qsea_result in qsea_outcome) {
  load(qsea_result)
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                      keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")
  result_v2<-get.ROIs_lg(result, regions)
  
  result_name <- str_sub(qsea_result, start=1, end = -7)
  qsea_sig_annot[[result_name]] <- list("sig" = sig, "annot" = result, "annot_ROIs" = result_v2)
}
save(qsea_sig_annot, file = "qsea_sig_annot_2021Jan28.RData")

get.QC_plots <- function(qseaSet_blind) {
  plist <- list()
  
  plist[["EPprofile"]] <- getOffset(qseaSet_blind, scale="fraction") # enrichment profile
  plist[["EPmatrix"]] <- plotEPmatrix(qseaSet_blind) # enrichment matrix
  plist[["CNV"]] <- plotCNV(qseaSet_blind) # a Heatmap-like Overview of the CNVs
  pca_cgi<-getPCA(qseaSet_blind, norm_method="beta")
  plist[["PCA"]] <- plotPCA(pca_cgi)
  
  return(plist)
}

QC_plots <- get.QC_plots(qseaSet_blind)
for(qsea_result in qsea_outcome) {
  load(qsea_result)
  QC_plots[[str_sub(qsea_result, start=1, end = -7)]] <- get.QC_plots(qseaSet_blind)
}
save(QC_plots, file = "qsea_QC_plots_2021Jan28.RData")