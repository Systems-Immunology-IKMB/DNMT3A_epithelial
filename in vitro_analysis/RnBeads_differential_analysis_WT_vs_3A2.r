#path to data.files
data.dir <- "/home/user/Projects/Antonella_EPIC/"

idat.dir <- file.path(data.dir, "iDats/")

analysis.dir <- file.path(data.dir, "Differential_analysis/WT_vs_3A2/")

report.dir <- file.path(analysis.dir, "reports/")

sample.annotation <- file.path(data.dir, "Annotation/annotation_file_WT_vs_3A2.csv")
data.type<-"idat.dir"

#source("https://bioconductor.org/biocLite.R")
#biocLite("RnBeads")
#
#
#

library(RnBeads)
library(LOLA)
library(GO.db)

help("rnb.options")
rnb.options(analysis.name = "example")
rnb.options(logging = TRUE)
#rnb.options(email = "")
rnb.options(assembly = "hg19")
rnb.options(identifiers.column = "Sample_ID")
rnb.options2xml(pretty=TRUE)
#colors.category = c("darkblue", "darkred")
rnb.options(min.group.size = 1)

###### modules
rnb.options(import = TRUE)
rnb.options(preprocessing = TRUE)
rnb.options(qc = TRUE)
rnb.options(inference = FALSE)
rnb.options(exploratory = FALSE)
rnb.options(differential = TRUE)

##############some selected options (see more at help("rnb.options"))

#rnb.options(import.default.data.type = "infinium.idat.dir") # NGS: "bs.bed.dir"
rnb.options(import.table.separator = ",") # for TABdel: "\t"
rnb.options(import.gender.prediction = FALSE) # only for array datasets

#normalization
rnb.options(normalization = NULL)
rnb.options(normalization.method = "wm.dasen")
rnb.options(normalization.background.method = "none")
rnb.options(normalization.plot.shifts = TRUE)

rnb.options(qc.boxplots = TRUE)
rnb.options(qc.barplots = TRUE)
rnb.options(qc.negative.boxplot = TRUE)
rnb.options(qc.snp.distances = TRUE)
rnb.options(qc.snp.boxplot = TRUE)
rnb.options(qc.snp.barplot = TRUE)
rnb.options(qc.sample.batch.size = 50)
rnb.options(qc.coverage.plots = FALSE)
rnb.options(qc.coverage.threshold.plot = 1:10)
rnb.options(qc.coverage.histograms = FALSE)
rnb.options(qc.coverage.violins = FALSE)

rnb.options(filtering.whitelist = NULL)
rnb.options(filtering.blacklist = NULL)
rnb.options(filtering.snp = "3")
rnb.options(filtering.cross.reactive = FALSE)
rnb.options(filtering.greedycut = TRUE)
rnb.options(filtering.greedycut.pvalue.threshold = 0.05)
rnb.options(filtering.greedycut.rc.ties = "row")
rnb.options(filtering.sex.chromosomes.removal = TRUE)
rnb.options(filtering.missing.value.quantile = 0.8)
rnb.options(filtering.coverage.threshold = 5)
rnb.options(filtering.low.coverage.masking = FALSE)
rnb.options(filtering.high.coverage.outliers = FALSE)
rnb.options(filtering.deviation.threshold = 0)

rnb.load.annotation.from.db("ensembleRegBuildBPall")
rnb.load.annotation.from.db("ensembleRegBuildBPproximal")
rnb.load.annotation.from.db("ensembleRegBuildBPdistal")
rnb.options(region.types=c("promoters","genes","cpgislands","tiling","ensembleRegBuildBPall","ensembleRegBuildBPproximal","ensembleRegBuildBPdistal")) 


rnb.options(differential.site.test.method = "ttest") #limma"
rnb.options(differential.comparison.columns = c("Sample_Group"))
#rnb.options(differential.comparison.columns.all.pairwise = c("CTRL_Vs_MSA_Npos",
#                                                            "ExVivoNpos_Vs_InVitroNpos",
#                                                           "Npos_Vs_Nneg"))

rnb.options(covariate.adjustment.columns = NULL)
#rnb.options(columns.pairing = c("Timepoint"="Patient"))
rnb.options(differential.adjustment.sva = TRUE)
rnb.options(differential.adjustment.celltype = FALSE)
rnb.options(differential.enrichment.go = TRUE)
rnb.options(differential.enrichment.lola = TRUE)
rnb.options(differential.variability = FALSE)
rnb.options(differential.report.sites = TRUE)

rnb.options(export.to.bed = TRUE)
rnb.options(export.to.trackhub = c("bigBed","bigWig"))
rnb.options(export.to.csv = TRUE) #default = FALSE #creates big tables
rnb.options(export.to.ewasher = FALSE)
rnb.options(export.types = "sites")

rnb.options(disk.dump.big.matrices = TRUE)
rnb.options(enforce.memory.management = TRUE)
rnb.options(enforce.destroy.disk.dumps = TRUE)

##############################################################################

##### Run Analysis
rnb.run.analysis(dir.reports=report.dir,
                 sample.sheet=sample.annotation,
                 data.dir=idat.dir,
                 data.type=data.type)


##############################################################################
