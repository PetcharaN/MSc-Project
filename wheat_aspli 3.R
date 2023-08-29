##################################################
#ASpli quick start tutorial using wheat BAM files#
##################################################

#Requirements
library(ASpli)
library(GenomicFeatures)

#Set WD
setwd("~/Project/Wheat ASpli")

# gff preprocessing
core_clock_genes <- "gff/core_clock_genes.gff3"
genomeTxDb <- makeTxDbFromGFF(core_clock_genes)

#feature extraction
features <- binGenome(genomeTxDb)
#* Number of extracted Genes = 78
#* Number of extracted Exon Bins = 553
#* Number of extracted intron bins = 497
#* Number of extracted trascripts = 165
#* Number of extracted junctions = 436
#* Number of AS bins (not include external) = 73
#* Number of AS bins (include external) = 73
#* Classified as: 
#  ES bins = 9	(12%)
#IR bins = 16	(22%)
#Alt5'ss bins = 13	(18%)
#	Alt3'ss bins = 20	(27%)
#Multiple AS bins = 15	(21%)
#classified as:
#  ES bins = 3	(20%)
#IR bins = 2	(13%)
#Alt5'ss bins = 2	(13%)
#			Alt3'ss bins = 7	(47%)


#bams and target file

#Set wd to bam folder
setwd("~/Project/Wheat ASpli/filtered_bams")

BAMfiles <- list("sort_trim_S1_filter2.bam.core.bam",
"sort_trim_S2_filter2.bam.core.bam",
"sort_trim_S3_filter2.bam.core.bam",
"sort_trim_S4_filter2.bam.core.bam",
"sort_trim_S5_filter2.bam.core.bam",
"sort_trim_S6_filter2.bam.core.bam",
"sort_trim_S7_filter2.bam.core.bam",
"sort_trim_S8_filter2.bam.core.bam",
"sort_trim_S9_filter2.bam.core.bam",
"sort_trim_S10_filter2.bam.core.bam",
"sort_trim_S11_filter2.bam.core.bam",
"sort_trim_S12_filter2.bam.core.bam",
"sort_trim_S13_filter2.bam.core.bam",
"sort_trim_S14_filter2.bam.core.bam",
"sort_trim_S15_filter2.bam.core.bam",
"sort_trim_S16_filter2.bam.core.bam",
"sort_trim_S17_filter2.bam.core.bam",
"sort_trim_S18_filter2.bam.core.bam",
"sort_trim_S19_filter2.bam.core.bam",
"sort_trim_S20_filter2.bam.core.bam",
"sort_trim_S21_filter2.bam.core.bam",
"sort_trim_S22_filter2.bam.core.bam",
"sort_trim_S23_filter2.bam.core.bam",
"sort_trim_S24_filter2.bam.core.bam",
"sort_trim_S25_filter2.bam.core.bam",
"sort_trim_S26_filter2.bam.core.bam",
"sort_trim_S27_filter2.bam.core.bam",
"sort_trim_S28_filter2.bam.core.bam",
"sort_trim_S29_filter2.bam.core.bam",
"sort_trim_S30_filter2.bam.core.bam",
"sort_trim_S31_filter2.bam.core.bam",
"sort_trim_S32_filter2.bam.core.bam",
"sort_trim_S33_filter2.bam.core.bam",
"sort_trim_S34_filter2.bam.core.bam",
"sort_trim_S35_filter2.bam.core.bam",
"sort_trim_S36_filter2.bam.core.bam",
"sort_trim_S37_filter2.bam.core.bam",
"sort_trim_S38_filter2.bam.core.bam",
"sort_trim_S39_filter2.bam.core.bam",
"sort_trim_S40_filter2.bam.core.bam",
"sort_trim_S41_filter2.bam.core.bam",
"sort_trim_S42_filter2.bam.core.bam",
"sort_trim_S43_filter2.bam.core.bam",
"sort_trim_S44_filter2.bam.core.bam",
"sort_trim_S45_filter2.bam.core.bam",
"sort_trim_S46_filter2.bam.core.bam",
"sort_trim_S47_filter2.bam.core.bam",
"sort_trim_S48_filter2.bam.core.bam",
"sort_trim_S49_filter2.bam.core.bam",
"sort_trim_S50_filter2.bam.core.bam",
"sort_trim_S51_filter2.bam.core.bam",
"sort_trim_S52_filter2.bam.core.bam",
"sort_trim_S53_filter2.bam.core.bam",
"sort_trim_S54_filter2.bam.core.bam")

#set wd back to project folder
setwd("~/Project/Wheat ASpli")

targets <- read.table("targets_wheat_5.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

#To find the mBAM merged conditions
getConditions(targets)

mBAMs <- read.table("mBAM_4.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

#wd back to bam folder
setwd("~/Project/Wheat ASpli/filtered_bams")

gbcounts <- gbCounts(features=features, targets=targets,minReadLength = 100, maxISize = 50000)
gbcounts

#Object of class ASpliCounts 
#Gene counts: 78 genes analysed. Access using countsg(object) 
#Gene RD: 78 genes analysed. Access using rdsg(object) 
#Bin counts: 975 bins analysed. Access using countsb(object) 
#Bin RD: 975 bins analysed. Access using rdsb(object) 
#Junction counts: 2938 junctions analysed. Access using countsj(object) 

asd <- jCounts(counts=gbcounts, features=features, minReadLength = 100)
asd

#Object of class ASpliAS 
#IR PIR:  440 intron bins analysed. 348 intron bins passed the filters. Access using irPIR(object) 
#ES PSI: 312 exon bins analysed. 20 exon bins passed the filters.  Access using esPSI(object) 
#AltSS PSI: 42 exon bins analysed. 32 exon bins passed the filters.  Access using altPSI(object) 
#Junctions PIR: 462 junctions analysed. Access using junctionsPIR(object) 
#Junctions PJU: 2767 junctions analysed. Access using junctionsPJU(object) 

#Differential gene expression and bin usage signal estimation:
gb  <- gbDUreport(gbcounts, contrast = c(24,4,8,12,16,20))
gb

#Differential junction usage analysis
jdur <- jDUreport(asd, contrast=c(24,4,8,12,16,20))
jdur

#Bin and junction signal integration:
sr <- splicingReport(gb, jdur, counts=gbcounts)

#Summary of integration of splicing signals along genomic-regions.
is <- integrateSignals(sr,asd)

#Export results:
exportIntegratedSignals(is,sr=sr,
                        output.dir = "wheat_aspli_24",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs,makeGraphs = T, bforce = T)
y

exportSplicingReports( sr,
                       output.dir="sr_24",
                       openInBrowser = FALSE,
                       maxBinFDR = 0.2,
                       maxJunctionFDR = 0.2 )
