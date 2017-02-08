library(tidyverse)

##Illumina kraken files:
ill.krak="/atium/Data/NGS/Aligned/170112_VRE/Illumina_CARD"
plotdir="/atium/Data/NGS/Aligned/170112_VRE/compare_CARD"

nano.krak="/atium/Data/NGS/Aligned/170112_VRE/Nanopore_CARD"
ref.krak = "/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"

comp.krak="/atium/Data/NGS/Aligned/170112_VRE/compare_CARD"
samp.list=c("VRENC_CARD","VRENC_nanopore_CARD", "VRE7_CARD", "VRE7_nanopore_CARD","VRE10_CARD","VRE10_nanopore_CARD")
samp.list=paste0("VRE", 1:10,"_CARD") #for Illumina
samp.list = c("VRENC_nanopore_CARD","VRE3_nanopore_CARD","VRE5_nanopore_CARD", "VRE6_nanopore_CARD", "VRE7_nanopore_CARD","VRE10_nanopore_CARD")
#samp.list = c("ENFM_withhuman","ENFS_withhuman","KPN_withhuman","VRENC_withhuman","VREPC1_withhuman","VREPC2_withhuman","VREPC3_withhuman")

krak.report=list()
krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

krak.report=bind_rows(lapply(samp.list, function(x) {tbl_df(data.frame(read.delim(
                                                         file.path(ill.krak, paste0(x, ".kraken.output")),
                                                         header=F, col.names=krak.report.cnames), sample=x))}))

##Trim left side white space in sci.name
krak.report$sci.name=trimws(krak.report$sci.name)

##ok - for first pass, let's try filtering by "species" or S in rank.code
##per.reads is of total reads, not the subsection of bacteria or just classified, as opposed to how pavian does it
##We also want .clade reads because that takes all the subspecies/strains - otherwise it will undercall because of subspecies


##Sort for top 3?
krak.top=filter(krak.report,rank.code!= "U", ncbi.tax.id!=1, ncbi.tax.id!=1187, ncbi.tax.id!=969, ncbi.tax.id!=3065, num.reads.taxon!=0 ) %>%
    group_by(sample) %>%
    top_n(n=5, wt=num.reads.taxon) %>%
    ungroup %>%
    select(sci.name, ncbi.tax.id) %>%
    distinct

##ok - of these top hits, get the max num reads per
krak.report %>%
    filter(ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
    group_by(sci.name) %>%
    summarise(big.reads=max(num.reads.clade))

#krak.report %>%
#    filter(ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
#    group <- by(sci.name) %>%
#        summarise(big.reads=max(num.reads.clade))


##ok - per group other bacteria calc
##ncbi.tax.id 2 is all bacteria
plot.res=filter(krak.report, ncbi.tax.id %in% c(2,krak.top$ncbi.tax.id)) %>%
    group_by(sample) %>%
    summarise(sci.name="other AMR",
              num.reads.clade=sum(num.reads.clade[ncbi.tax.id!=0])-sum(num.reads.clade[ncbi.tax.id %in% krak.top$ncbi.tax.id]))
#              per.reads=100*(sum(per.reads[ncbi.tax.id!=0])-sum(per.reads[ncbi.tax.id %in% krak.top.withper$ncbi.tax.id]))/sum(per.reads[ncbi.tax.id!=0]))

##per group human calc (skip for now)

##20129 is vanA, 20371 is vanB, 1073 is KPC (general)


plot.res=bind_rows(plot.res,
                   filter(krak.report, ncbi.tax.id==20129) %>%
                   select(sample, sci.name, num.reads.taxon))

plot.res=bind_rows(plot.res,
                      filter(krak.report, ncbi.tax.id==20371) %>%
                      select(sample, sci.name, num.reads.taxon))

plot.res=bind_rows(plot.res,
                      filter(krak.report, ncbi.tax.id==1073) %>%
                      select(sample, sci.name, num.reads.taxon))

##per group unclass hits
##Unclass just means save U
plot.res=bind_rows(plot.res,
                   filter(krak.report, rank.code=="U") %>%
                   select(sample, sci.name, num.reads.clade))


##Finally our actual reads of interest
plot.res=bind_rows(plot.res,
                   filter(krak.report, ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
                   select(sample, sci.name, num.reads.clade))

##okay, overpowered by all of the unidentified reads, let's just try without the unidentified
plot.res_old <- plot.res[ which(plot.res$rank.code!="U"),]

##Sort by rank order, also generate % reads per samp
plot.res_old=plot.res_old%>%
    group_by(sample) %>%
    arrange(sample, desc(num.reads.taxon)) %>%
    mutate(per.reads=round(100*num.reads.taxon/sum(num.reads.taxon), 2))

plot.res[is.na(plot.res)]<- 0

plot.res=plot.res %>%
    group_by(sample) %>%
    arrange(sample, desc(num.reads.clade)) %>%
##    mutate(per.reads=100*num.reads.clade/sum(num.reads.clade))
        mutate(per.reads=round(100*num.reads.clade/sum(num.reads.clade), 2))

##Double check per.reads adds up to 100
plot.res%>%
    group_by(sample) %>%
    summarize(n=sum(per.reads))
##Pretty close - rounding errors only

##fix x-axis order so VRE10 on the end
#plot.res_old$sample=factor(plot.res_old$sample, levels=samp.list)

plot.res$sample=factor(plot.res$sample,levels=samp.list)
##Better color scheme still
##Human reads
##reorder so that blocks are all sorted in same order - some order we define, like unclass, other bac, archae, virus,
##then our bac of interest
##Italize bacteria
pdf(file.path(plotdir, "heat_ill_CARD_top10_wVRE.pdf"), width=24, height=15)

ggplot(plot.res, aes(sample, sci.name)) + geom_tile(aes(fill=per.reads), color= "white") +
    scale_fill_gradient(low="white", high="black") + ylab("AMR genes") + xlab("Samples") + theme(panel.background = element_rect(fill='lightgoldenrod'))
#ggplot(plot.res, aes(x=sample, y=per.reads, fill=sci.name))+geom_bar(color="black", stat='identity')+
#    scale_fill_hue()+theme_bw()+theme(legend.position="bottom")

dev.off()


