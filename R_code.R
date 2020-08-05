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
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/ngs-data-2/analysis/NhanNguyen/MeDIP/ConDF2/Alignments/")
sample_ConDF2=data.frame(sample_name=c(paste0(rep("Con_DF2_T000_", 3), 1:3), paste0(rep("Con_DF2_T002_", 3), 1:3)),
                        file_name=list.files(),
                        group=c(rep("T000", 3), rep("T002", 3)), stringsAsFactors=FALSE)

qseaSet=createQseaSet(sampleTable=sample_ConDF2, 
                      BSgenome="BSgenome.Hsapiens.UCSC.hg38",
                      window_size=100)
# Switch for parallel computing, using BiocParallel
library("BiocParallel")
register(MulticoreParam(workers=3))

# Import sequencing data
qseaSet=addCoverage(qseaSet, uniquePos=TRUE, paired=TRUE, parallel=TRUE)

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

