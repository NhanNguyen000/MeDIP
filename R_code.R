#source("http://bioconductor.org/biocLite.R")
#biocLite("MEDIPS")

# check requirement package
packageDescription("BSgenome")
packageDescription("gtools")
packageDescription("edgeR")
packageDescription("DNAcopy")
packageDescription("Rsamtools")
packageDescription("rtracklayer")

# selected genome:
library("BSgenome")
available.genomes()

#We mapped the short reads against the human genome build hg38. Therefore, we download and install this genome build:
source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

# data toy
source("http://bioconductor.org/biocLite.R")
#biocLite("MEDIPSData")

# Load the MEDIPS package
library(MEDIPS)
# we mapped the short reads against the human genome build hg38 -> load the hg38 library:
#library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

# load the data with bam files
library("MEDIPSData")
bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam",
                                        package = "MEDIPSData")
bam.file.hESCs.Input = system.file("extdata", "hESCs.Input.chr22.bam",
                                   package = "MEDIPSData")
bam.file.DE.Input = system.file("extdata", "DE.Input.chr22.bam",
                                package = "MEDIPSData")

# set parametersL
# The reference genome is hg38:
#BSgenome="BSgenome.Hsapiens.UCSC.hg38"
BSgenome="BSgenome.Hsapiens.UCSC.hg19"

#MEDIPS will replace all reads which map to exactly the same start and end
#positions on the same strand by only one representative:
uniq=TRUE

#All reads will be extended to a length of 300nt according to the given strand information:
extend=300

# As an alternative to the extend parameter, the shift parameter can be specified. Here, the reads are not extended but shifted by the specified number
# of nucleotides with respect to the given strand infomation. One of the two
# parameters extend or shift has to be 0.
shift=0

#The genome will be divided into adjacent windows of length 100nt and all further
#calculations (short read coverage, differential coverage between conditions etc.)
#will be calculated for each window.
ws=100

#In this manual, we are going to process only one chromosome (i.e. chr22) in
#order to exemplify a typical workflow. Therefore, we specify the chr.select
#parameter. Please note, the example bam files contain only data for chr22
#anyway. Therefore, it is not necessary to explicitly specify this chromosome in
#this example.
chr.select="chr22"


#  Quality controls
# MEDIPS provides three different quality controls
# 5.1 Saturation analysis - could not run this codeline
sr = MEDIPS.saturation(file = bam.file.hESCs.Rep1.MeDIP, BSgenome = BSgenome,
                       uniq = uniq, extend = extend, shift = shift, window_size = ws,
                       chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
                       rank = FALSE)
Check code in: https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/MEDIPS/inst/doc/MEDIPS.pdf
