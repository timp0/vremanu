library(tidyverse)

##Illumina kraken files:
ill.krak="/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"
ill.krak.CARD="/atium/Data/NGS/Aligned/170112_VRE/Illumina_CARD"
statdir="/atium/Data/NGS/Aligned/170112_VRE/"
nano.krak="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/report"

samp.list = c("VRE10_withhuman", "VRE10_nanopore")
krak.report=list()

ill.krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

nano.krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

nano.krak.report=bind_rows(lapply(samp.list, function(x) {tbl_df(data.frame(read.delim(
                                                                 file.path(nano.krak, paste0("VRE10_nanopore.kraken.output")),
                                                                 header=F, col.names=nano.krak.report.cnames), sample=x))}))

ill.krak.report=bind_rows(lapply(samp.list, function(x) {tbl_df(data.frame(read.delim(
                                                                 file.path(ill.krak, paste0("VRE10_withhuman.kraken.output")),
                                                                 header=F, col.names=ill.krak.report.cnames), sample=x))}))

#ill.krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

#nano.krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

g.ill.filt = filter(ill.krak.report, rank.code=="G")
s.ill.filt = filter(ill.krak.report, rank.code=="S")

g.nano.filt = filter(nano.krak.report, rank.code=='G') #weirdly it thinks that it's rank.c.ode instead of rank.code
s.nano.filt = filter(nano.krak.report, rank.code=='S')

##okay, so performing either Bray-Curtis or spearman requires the same number of reads for each sample, so will probably compare ncbi.tax.id. If they are the same, append to a new list, if not, then we shall just ignore it and move on
##compared.vec <- matrix(,,3)
grouped <- match(g.nano.filt$ncbi.tax.id, g.ill.filt$ncbi.tax.id) #nano is longer
nonagrouped<-na.omit(grouped)
compared.vec <- matrix(,length(nonagrouped),3)
nanoval <- which(!is.na(grouped), arr.ind=TRUE)

for (i in 1:length(nonagrouped)) {
    compared.vec[i,1] <- g.ill.filt$ncbi.tax.id[nonagrouped[i]]
    compared.vec[i,2] <- g.ill.filt$num.reads.taxon[nonagrouped[i]]
    compared.vec[i,3] <- g.nano.filt$num.reads.taxon[nanoval[i]]
}

cor.test(compared.vec[,2],compared.vec[,3], method='pearson') # test with pearson first





#samp.list=tibble(filey=file.path(ill.krak, paste0("VRE", 1:10, "_CARD.kraken.output")),
#                 label=paste0("VRE", 1:10, "Illumina"))

#plot_stacked(samp.list, plotdir, 'stack_ill_CARD')

#nano.krak="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/report"
#nano.krak.CARD="/atium/Data/NGS/Aligned/170112_VRE/Nanopore_CARD"

#samp.list=tibble(filey=file.path(nano.krak, paste0("VRE", c("NC", 3, 5, 6, 7, 10), "_nanopore_CARD.kraken.output")),
#                 label=paste0("VRE", c("NC", 3, 5, 6, 7, 10), "Nanopore"))

#plot_stacked(samp.list, plotdir, 'stack_nano_CARD')
