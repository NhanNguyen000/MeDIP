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

folder <- "/ngs-data-2/analysis/NhanNguyen/MeDIP/Hepatic/"

data_files <- as.data.frame(matrix(NA, ncol = 3, nrow = 24))
colnames(data_files) <- c("file_name", "Condition", "Time")
data_files$file_name <- c(paste0("Cyclosporin_", c(784, 785, 786, seq(796, 801, 1), seq(949, 954, 1)), "_pe.sorted.bam"),
                          paste0("Control_", c(421, 422, 423), "_pe.sorted.bam"), 
                          paste0("ConDMSO_", c(454, 455, 456), "_pe.sorted.bam"), 
                          paste0("ConDMSO_", c(457, 458, 459), "_pe.sorted.bam"))
data_files$Condition <- c(rep("Cyc_000", 3), rep("Cyc_The", 6), rep("Cyc_Tox", 6), rep("Control", 9))
data_files$Time <- c(rep("000", 3), rep("072", 3), rep("168", 3), 
                     rep("072", 3), rep("168", 3), 
                     rep("000", 3), rep("072", 3), rep("168", 3))

data_number <- as.data.frame(matrix(NA, nrow = 3, ncol = 3))
rownames(data_number) <-c("000", "072", "168")
colnames(data_number) <- c("Control", "Cyc_The", "Cyc_Tox")

for(time_point in c("072", "168")) {
  for(Condition in c("Cyc_The", "Cyc_Tox")) {
    samples<-data.frame(sample_name=c(data_files$file_name[data_files$Condition == Condition & data_files$Time== time_point],
                                      data_files$file_name[data_files$Condition == "Control" & data_files$Time== time_point]),
                        file_name=paste0(folder,
                                         c(data_files$file_name[data_files$Condition == Condition & data_files$Time== time_point],
                                           data_files$file_name[data_files$Condition == "Control" & data_files$Time== time_point])),
                        group=c(rep("Cyc", 3), rep("Control", 3)), stringsAsFactors=FALSE)
    qsea.get_DMR(samples, output_name = paste0(Condition, time_point))
  }
}

#install.packages("tidyverse")
library(tidyveres)

data_The <- data_files %>% subset(!Condition == "Cyc_Tox")
The_samples <- data.frame(sample_name = data_The$file_name, 
                          file_name = paste0(folder, data_The$file_name),
                          group = c(rep("Cyc_The", 9), rep("Control", 9)), stringsAsFactors = FALSE)
qsea.get_DMR(The_samples, output_name = "Cyc_The_allSamples")

data_Tox <- data_files %>% subset(!Condition == "Cyc_The")
Tox_samples <- data.frame(sample_name = data_Tox$file_name,
                          file_name = paste0(folder, data_Tox$file_name),
                          group = c(rep("Cyc_Tox", 9), rep("Control", 9)), stringsAsFactors = FALSE)
qsea.get_DMR(Tox_samples, output_name = "Cyc_Tox_allSamples")

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


qsea_outcome <- list.files()[grepl("qsea_outcome_Cyc" ,list.files())]
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

## Check the gene region: ---------------------------------------------------
load("ROIs_2021Jan19.RData")
load("ROIs_2_2021Jan19.RData")
library(GenomicRanges)

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
qsea_outcome <- list.files()[grep("qsea_outcome_Cyc",list.files())]
for (i in 1:length(qsea_outcome)) {
  load(qsea_outcome[i])
  sig <- isSignificant(qseaGLM, fdr_th=0.01)
  result <- makeTable(qseaSet_blind, glm=qseaGLM, groupMeans=getSampleGroups(qseaSet_blind), 
                      keep=sig, annotation=c(ROIs, ROIs_2), norm_method="beta")
  qsea_outcome_region[[qsea_outcome[i]]] <-get.sum_region(result) 
}

save(qsea_outcome_region, file = "qsea_outcome_region_2021April.RData")
