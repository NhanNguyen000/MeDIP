#install.packages("/project/modifiers/MEDIPS/qsea_0.0.8.tar.gz")
#detach("package:qsea", unload=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(qsea)
library(limma)
library(GenomicRanges)
date=format(Sys.time(), "%Y%m%d")
dir="/project/epitreat/user/lienhard/"


###############
##Annotations##
###############

refseq_tab=read.table("/project/42/references/refseq/RefSeq_GRCh37_hg19_known_genes_150812.txt")
tss=refseq_tab$V5
tss[refseq_tab$V4=="-"]=refseq_tab$V6[refseq_tab$V4=="-"]
tss_reg=GRanges(refseq_tab$V3, IRanges(tss, width=1), strand=refseq_tab$V4, gene_name=refseq_tab$V13, refseq=refseq_tab$V2)
prom_reg=GRanges(refseq_tab$V3, IRanges(tss-2000, width=4000), strand=refseq_tab$V4, gene_name=refseq_tab$V13, refseq=refseq_tab$V2)

tfbs_tab=read.table("/project/42/references/TFBS/ENCODERegTFBSv3_140424_anno.txt", header=T, stringsAsFactors=F)
tfbs_tab$nextID=tfbs_tab$upstreamTSS_id
down=as.numeric(tfbs_tab$upstreamTSS_dist)+as.numeric(tfbs_tab$downstreamTSS_dist)<0
down[is.na(down)]=F
sum(down)
sum(!down)
tfbs_tab$nextID[down]=tfbs_tab$downstreamTSS_id[down]
tfbs_reg=GRanges(tfbs_tab$chrom, IRanges(tfbs_tab$chromStart, tfbs_tab$chromEnd), name=tfbs_tab$name, nextGene=tfbs_tab$nextID)

transcr_factors=unique(tfbs_tab$name)
tfbs=list()
for(tf in transcr_factors){
	idx=tfbs_tab$name==tf
	tfbs[[tf]]=GRanges(tfbs_tab$chrom[idx], IRanges(tfbs_tab$chromStart[idx], tfbs_tab$chromEnd[idx]), name=tfbs_tab$name[idx], nextGene=tfbs_tab$nextID[idx])
}

CGI_tab=read.table("/project/42/references/hg19/model-based-cpg-islands-hg19_rafalab.jhsph.edu.txt", header=T)
CGI_reg=GRanges(CGI_tab$chr, IRanges(CGI_tab$start, CGI_tab$end), obsExp=CGI_tab$obsExp)
gb_reg=GRanges(refseq_tab$V3, IRanges(refseq_tab$V5,refseq_tab$V6), strand=refseq_tab$V4, gene_name=refseq_tab$V13, refseq=refseq_tab$V2)
exon_tab=read.table("/project/42/references/refseq/RefSeq_GRCh37_hg19_known_genes_150812.tab", header=T)
exon_reg=GRanges(exon_tab$chr, IRanges(exon_tab$start,exon_tab$end),  gene_name=exon_tab$gene, refseq=exon_tab$transcript, exon_nr=exon_tab$exon_nr)

intron_reg=setdiff(gb_reg, exon_reg, ignore.strand=TRUE)
ol=findOverlaps(intron_reg, gb_reg,select="first")
which(is.na(ol))
intron_reg$gene_name=gb_reg$gene_name[ol]


#CGIprom_reg=intersect(CGI_reg,prom_reg , ignore.strand=TRUE)
CGIprom_reg=CGI_reg[overlapsAny(CGI_reg,prom_reg , ignore.strand=TRUE)]
CGItfbs_reg=intersect(CGI_reg,tfbs_reg , ignore.strand=TRUE)
CGItfbsprom_reg=intersect(CGIprom_reg,tfbs_reg , ignore.strand=TRUE)
ol=findOverlaps(CGIprom_reg, prom_reg,select="first")
CGIprom_reg$prom=prom_reg$gene_name[ol]
CGIdistal_reg=CGI_reg[!overlapsAny(CGI_reg,prom_reg , ignore.strand=TRUE)]

ol=findOverlaps(CGItfbsprom_reg, prom_reg,select="first")
promTFBS_reg=intersect(tfbs_reg,prom_reg , ignore.strand=TRUE)


CGItfbsprom_reg$prom=prom_reg$gene_name[ol]

promNoTFBS_reg=setdiff(prom_reg,tfbs_reg , ignore.strand=TRUE)
CGINoTFBS_reg=setdiff(CGI_reg,tfbs_reg , ignore.strand=TRUE)
CGINoprom_reg=setdiff(CGI_reg,prom_reg , ignore.strand=TRUE)
CGIshore=GRanges(seqnames(CGI_reg), IRanges(end=start(CGI_reg), width=2000))
CGIshore=reduce(c(CGIshore,GRanges(seqnames(CGI_reg), IRanges(start=end(CGI_reg), width=2000)) ))
CGIshore_prom=GRanges(seqnames(CGIprom_reg), IRanges(end=start(CGIprom_reg), width=2000))
CGIshore_prom=reduce(c(CGIshore_prom,GRanges(seqnames(CGIprom_reg), IRanges(start=end(CGIprom_reg), width=2000)) ))
CGIshore_distal=GRanges(seqnames(CGIdistal_reg), IRanges(end=start(CGIdistal_reg), width=2000))
CGIshore_distal=reduce(c(CGIshore_distal,GRanges(seqnames(CGIdistal_reg), IRanges(start=end(CGIdistal_reg), width=2000)) ))


GRVenn<-function(x,y){
	ol=sum(overlapsAny(x,y))	
	c(length(x), length(y), ol)
}
PRC=reduce(c(tfbs[["SUZ12"]],tfbs[["EZH2"]] ))
PRCbig=reduce(c(tfbs[["SUZ12"]], tfbs[["EZH2"]]))
start(PRCbig)=start(PRC)-500
end(PRCbig)=end(PRC)-500

#ol=sapply(tfbs, GRVenn, y=PRC)
#sort(ol[3,]/ol[1,])
#GRVenn(tfbs[["SUZ12"]], tfbs[["EZH2"]])
olCGI=sapply(tfbs, GRVenn, y=CGI_reg)
sort(olCGI[3,]/olCGI[1,])#fraction CGI overlap
olProm=sapply(tfbs, GRVenn, y=prom_reg)
sort(olProm[3,]/olProm[1,])#fraction Promoter overlap


nonPRC2tfbs=setdiff(reduce(tfbs_reg), PRCbig)
distal_prom
enh=read.table("/project/42/references/Ensemblv82/ens82_grch37_vista_enhancers.tsv",sep="\t", header=T)
vistaEnh=GRanges(enh[,1], IRanges(enh[,2], enh[,3]))
A549_reg=read.table("/project/42/references/Ensemblv82/ens82_grch37_A549_features.txt", sep="\t", header=T, stringsAsFactors=F)
#A549_reg=A549_reg[A549_reg$Chromosome.Name %in% 1:22,]
A549_peaks=list()
for(n in unique(A549_reg$Feature.Type)){
f=A549_reg$Feature.Type==n
message(n)
A549_peaks[[n]]=GRanges(paste0("chr",A549_reg[f,1]), IRanges(A549_reg[f,2],A549_reg[f,3]) ) 
}
A549_reg=sort(GRanges(paste0("chr",A549_reg[,1]), IRanges(A549_reg[,2],A549_reg[,3]), feature=A549_reg[,4] ) )

ROIs=list(gb_reg, exon_reg,intron_reg,vistaEnh,A549_reg,tfbs_reg,PRC, nonPRC2tfbs,prom_reg,promTFBS_reg, promNoTFBS_reg, CGI_reg, CGIshore,CGIprom_reg, CGIshore_prom, CGIdistal_reg, CGIshore_distal, CGItfbs_reg,CGINoTFBS_reg,CGItfbsprom_reg)
names(ROIs)=c("gene body", "exon", "intron", "enhancer", "ENCODE_A549","TFBS","PRC2 TFBS","non PRC2 TFBS",  "promoter","TFBS at promoter","nonTFBS promoter", "CGI","CGI shores" ,"promoter CGI","promoter CGI shores" ,"distal CGI","distal CGI shores" , "CGI TFBS","nonTFBS CGI","CGI TFBS at promoter")


save(ROIs,tfbs,A549_peaks, file=paste0(dir,"RData/epitreat_",date,"_annotation.RData"))

##############################
# methyl seq target regionen #
##############################
#CpG 


targets=read.table(file="/project/epitreat_data/MethylSeq/docs/MethylSeq.all.244K_wID_conc.bed", header=F, sep="\t", stringsAsFactors=FALSE)
names(targets)=c("seqnames","start", "end", "probes")
targets=targets[targets$seqnames %in%paste0("chr", 1:22),]
tReg=makeGRangesFromDataFrame(targets,seqinfo=seqinfo(getRegions(epitreatQS)))
tReg=reduce(tReg)


########################
# 450k target regionen #
########################
#perl -F',' -lane 'if ($.>7){print "$F[11]\t$F[12]\t$F[21]\t$F[23]\t$F[24]"}' /project/epitreat/user/lienhard/docs/HumanMethylation450_15017482_v1-2.csv>/project/epitreat/user/lienhard/docs/HumanMethylation450_15017482_v1-2.tab
probes=read.table(file="/project/epitreat/user/lienhard/docs/HumanMethylation450_15017482_v1-2.tab", sep="\t",header=T, stringsAsFactors=FALSE)
#probes=probes[probes$CHR %in% 1:22,]
iReg=GRanges(paste0("chr",probes$CHR), IRanges(probes$MAPINFO, width=2), name=UCSC_RefGene_Name)
g=getRegions(epitreatQS)
iWd=g[overlapsAny(g, iReg)]

targets=list(Illumina450k=iWd, MethylSeq=tReg)
save(targets, file=paste0(dir,"RData/epitreat_",date,"_methylation_assay_targets.RData"))
library(FDb.InfiniumMethylation.hg19)
hm450.hg19 <- getPlatform(platform='HM450', genome='hg19')


