library(tidyverse)

##Illumina kraken files:
ill.krak="/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"
plotdir="/atium/Data/NGS/Aligned/161110_VREkraken"

nano.krak="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/report"
ref.krak = "/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"
#samp.list=paste0("VRE", 1:10,"_withhuman") #for Illumina
samp.list = c("NC_nanopore","VRE3_nanopore","VRE5_nanopore", "VRE6_nanopore", "VRE7nanopore","VRE10nanopore")
#samp.list = c("ENFM_withhuman","ENFS_withhuman","KPN_withhuman","VRENC_withhuman","VREPC1_withhuman","VREPC2_withhuman","VREPC3_withhuman")

krak.report=list()
krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

krak.report=bind_rows(lapply(samp.list, function(x) {tbl_df(data.frame(read.delim(
                                                         file.path(nano.krak, paste0(x, ".kraken.output")),
                                                         header=F, col.names=krak.report.cnames), sample=x))}))

##Trim left side white space in sci.name
krak.report$sci.name=trimws(krak.report$sci.name)

##ok - for first pass, let's try filtering by "species" or S in rank.code
##per.reads is of total reads, not the subsection of bacteria or just classified, as opposed to how pavian does it
##We also want .clade reads because that takes all the subspecies/strains - otherwise it will undercall because of subspecies

##Sort for top 3?
krak.top=filter(krak.report, rank.code=="S") %>%
    group_by(sample) %>%
    top_n(n=2, wt=num.reads.clade) %>%
    ungroup %>%
    select(sci.name, ncbi.tax.id) %>%
    distinct

##ok - of these top hits, get the max num reads per
krak.report %>%
    filter(ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
    group_by(sci.name) %>%
    summarise(big.reads=max(num.reads.clade))


##ok - per group other bacteria calc
##ncbi.tax.id 2 is all bacteria
plot.res=filter(krak.report, ncbi.tax.id %in% c(2, krak.top$ncbi.tax.id)) %>%
    group_by(sample) %>%
    summarise(sci.name="Other bacteria",
              num.reads.clade=num.reads.clade[ncbi.tax.id==2]-sum(num.reads.clade[ncbi.tax.id %in% krak.top$ncbi.tax.id]))

##per group human calc (skip for now)

##per group non bac/non human hits
plot.res=bind_rows(plot.res,
                   filter(krak.report, rank.code=="D"& ncbi.tax.id!=2) %>%
                   select(sample, sci.name, num.reads.clade))

##per group unclass hits
##Unclass just means save U
plot.res=bind_rows(plot.res,
                   filter(krak.report, rank.code=="U") %>%
                   select(sample, sci.name, num.reads.clade))


##Finally our actual hits of interest
plot.res=bind_rows(plot.res,
                   filter(krak.report, ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
                   select(sample, sci.name, num.reads.clade))

##Sort by rank order, also generate % reads per samp
plot.res=plot.res %>%
    group_by(sample) %>%
    arrange(sample, desc(num.reads.clade)) %>%
    mutate(per.reads=round(100*num.reads.clade/sum(num.reads.clade), 2))

##Double check per.reads adds up to 100
plot.res %>%
    group_by(sample) %>%
    summarize(n=sum(per.reads))
##Pretty close - rounding errors only

##fix x-axis order so VRE10 on the end
plot.res$sample=factor(plot.res$sample, levels=samp.list)


##Better color scheme still
##Human reads
##reorder so that blocks are all sorted in same order - some order we define, like unclass, other bac, archae, virus,
##then our bac of interest
##Italize bacteria
pdf(file.path(plotdir, "stack2.pdf"), width=11, height=8.5)

##ggplot(plot.res, aes(x=sample, y=num.reads.clade, fill=sci.name))+geom_bar(stat='identity')+theme_bw()
ggplot(plot.res, aes(x=sample, y=per.reads, fill=sci.name))+geom_bar(color="black", stat='identity')+
    scale_fill_hue()+theme_bw()+theme(legend.position="bottom")

dev.off()


