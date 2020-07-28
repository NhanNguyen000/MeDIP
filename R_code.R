if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("MEDIPS")

# avaiblable genome
library("BSgenome")
available.genomes()



# The reference genome is hg38:
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# Play with toy data -----------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("MEDIPSData")

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library("MEDIPSData")
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam",
                                        package = "MEDIPSData")
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
hESCs_MeDIP = MEDIPS.createSet(file = bam.file.hESCs.Rep1.MeDIP,
                               BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
                               window_size = ws, chr.select = chr.select)

bam.file.hESCs.Rep2.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep2.chr22.bam",
                                        package = "MEDIPSData")
hESCs_MeDIP = c(hESCs_MeDIP, MEDIPS.createSet(file = bam.file.hESCs.Rep2.MeDIP,
                                              BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
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
                                      + annotation = c("GENE"), chr = "chr22")

mr.edgeR.s = MEDIPS.setAnnotation(regions = mr.edgeR.s, annotation = anno.mart.gene)

# addCNV (copy number variation)
mr.edgeR = MEDIPS.addCNV(cnv.Frame = 10000, ISet1 = hESCs_Input,
                         ISet2 = DE_Input, results = mr.edgeR)

# Calibration Plot --? seem to be improtant:
MEDIPS.plotCalibrationPlot(CSet = CS, main = "Calibration Plot",
                           MSet = hESCs_MeDIP[[1]], plot_chr = "chr22", rpkm = TRUE,
                           + xrange = TRUE)
