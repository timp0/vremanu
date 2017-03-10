library(tidyverse)

##Illumina kraken files:
ill.krak="/atium/Data/NGS/Aligned/170112_VRE/Illumina_CARD"
plotdir="/atium/Data/NGS/Aligned/170112_VRE/compare_CARD"

nano.krak="/atium/Data/NGS/Aligned/170112_VRE/Nanopore_CARD"
ref.krak = "/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"

comp.krak="/atium/Data/NGS/Aligned/170112_VRE/compare_CARD"
#samp.list=c("VRENC_CARD","VRENC_nanopore_CARD", "VRE7_CARD", "VRE7_nanopore_CARD","VRE10_CARD","VRE10_nanopore_CARD")
samp.list=paste0("VRE", 1:10,"_CARD") #for Illumina
#samp.list = c("VRENC_nanopore_CARD","VRE3_nanopore_CARD","VRE5_nanopore_CARD", "VRE6_nanopore_CARD", "VRE7_nanopore_CARD","VRE10_nanopore_CARD")
#samp.list = c("ENFM_withhuman","ENFS_withhuman","KPN_withhuman","VRENC_withhuman","VREPC1_withhuman","VREPC2_withhuman","VREPC3_withhuman")

krak.report=list()
krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

krak.report=bind_rows(lapply(samp.list, function(x) {tbl_df(data.frame(read.delim(
                                                         file.path(ill.krak, paste0(x, ".kraken.output")),
                                                         header=F, col.names=krak.report.cnames), sample=x))}))

##Trim left side white space in sci.name
krak.report$sci.name=trimws(krak.report$sci.name)

##So, let's look at major AMR groups - from CARD, the determinants of antibiotic resistance - listing ncbi.tax.id's for the genes that show up in any of the ten samples atm
##antibiotic efflux genes - 969
##antibiotic target protection protein - 3107
##gene modulating beta-lactam resistance - 2046
##antibiotic resistance gene cluster - 2496
##tetracyline - 1133
##aminocoumarin resistance - 2034
##protein modulating permeability to antibiotic - 173
##antibiotic inactivation enzyme - 2175
##fosfomycin resistance gene -- 77
##antibiotic resistant gene variant or mutant - 886
##glycopeptide - 2965
##aminoglycoside - 2577
##antibiotic target replacement protein - 2415
##mupirocin - 1729
##gene modulating antibiotic efflux - 629
##peptide antibiotic resistance proteins - 3320
##polymyxin (colistin in this group) - 1388
##altering cell wall charge - 2799

##these were not in the samples but of the major subgroups of CARD
##macrolide resistance protein - 161
##beta-lactam resistance protein - not in database?, but there's already a group of genes that are beta-lactam resistant
##streptogramin - ARO not in database
##rifamycin - 3484
##sulfonamide - ARO not in database
##gene involved in self-resistance to antibiotic - 1113
##antibiotic target modifying enzyme - ARO not in database
##phenicol - ARO not in database
##lincosamide resistance - ARO not in database
##linezolid - ARO not in database
##AMR via molecular bypass - 1253
##AMR screening microarray - ARO not in database
##fluroquinoline - ARO not in database
##gene involved in antibiotic sequestration - 658
##elfamycin - ARO not in database
##fusidic acid - ARO not in database
##lipopeptide - ARO not in database
##isoniazid - ARO not in database
##pyrazinamide - ARO not in database
##polyamine - ARO not in database
##triclosan - ARO not in database
##oxazolidinone - ARO not in database
##pleuromutilin - ARO not in database
##nitrofuratoin - ARO not in database
##gene conferring resistance via absence - ARO not in database
##nybomycin - ARO not in database
##nucleoside - ARO not in database
##diaminopyrimidine - ARO not in database


##ok - for first pass, let's try filtering by "species" or S in rank.code
##per.reads is of total reads, not the subsection of bacteria or just classified, as opposed to how pavian does it
##We also want .clade reads because that takes all the subspecies/strains - otherwise it will undercall because of subspecies


##Sort for top 3?
krak.top=filter(krak.report,ncbi.tax.id==969| ncbi.tax.id==3107 | ncbi.tax.id==2046 | ncbi.tax.id==2496 | ncbi.tax.id==1133 | ncbi.tax.id==2034 | ncbi.tax.id==173 | ncbi.tax.id==2175 | ncbi.tax.id==77 |ncbi.tax.id==886 |ncbi.tax.id==886 | ncbi.tax.id==2965 | ncbi.tax.id==2577 | ncbi.tax.id==2415 | ncbi.tax.id==1729 | ncbi.tax.id==629 | ncbi.tax.id==3320 | ncbi.tax.id==1388 | ncbi.tax.id==2799 | ncbi.tax.id==161 | ncbi.tax.id == 3484 | ncbi.tax.id==1113 | ncbi.tax.id == 1253 | ncbi.tax.id==658) %>%
    group_by(sample) %>%
#    top_n(n=5, wt=num.reads.taxon) %>%
#    ungroup %>%
    select(sample,sci.name, ncbi.tax.id, num.reads.clade) #%>%
#    distinct

##ok - of these top hits, get the max num reads per
#krak.report %>%
#    filter(ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
#    group_by(sci.name) %>%
#    summarise(big.reads=max(num.reads.clade))

#krak.report %>%
#    filter(ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
#    group <- by(sci.name) %>%
#        summarise(big.reads=max(num.reads.clade))


##ok - per group other bacteria calc
##ncbi.tax.id 2 is all bacteria
#plot.res=filter(krak.report, ncbi.tax.id %in% c(2,krak.top$ncbi.tax.id)) %>%
#    group_by(sample) %>%
#    select(sample, sci.name, num.reads.clade)
#    summarise(sci.name="other AMR",
#              num.reads.clade=sum(num.reads.clade[ncbi.tax.id!=0])-sum(num.reads.clade[ncbi.tax.id %in% krak.top$ncbi.tax.id]))
#              per.reads=100*(sum(per.reads[ncbi.tax.id!=0])-sum(per.reads[ncbi.tax.id %in% krak.top.withper$ncbi.tax.id]))/sum(per.reads[ncbi.tax.id!=0]))

##per group human calc (skip for now)

##20129 is vanA, 20371 is vanB, 1073 is KPC (general)


#plot.res=bind_rows(plot.res,
#                   filter(krak.report, ncbi.tax.id==20129) %>%
#                   select(sample, sci.name, num.reads.taxon))

#plot.res=bind_rows(plot.res,
#                      filter(krak.report, ncbi.tax.id==20371) %>%
#                      select(sample, sci.name, num.reads.taxon))

#plot.res=bind_rows(plot.res,
#                      filter(krak.report, ncbi.tax.id==1073) %>%
#                      select(sample, sci.name, num.reads.taxon))

##per group unclass hits
##Unclass just means save U
#plot.res=bind_rows(plot.res,
#                   filter(krak.report, rank.code=="U") %>%
#                   select(sample, sci.name, num.reads.clade))


##Finally our actual reads of interest
#plot.res=bind_rows(plot.res,
#                   filter(krak.report, ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
#                   select(sample, sci.name, num.reads.clade))

##okay, overpowered by all of the unidentified reads, let's just try without the unidentified
#plot.res_old <- plot.res[ which(plot.res$rank.code!="U"),]

##Sort by rank order, also generate % reads per samp
#plot.res_old=plot.res_old%>%
#    group_by(sample) %>%
#    arrange(sample, desc(num.reads.taxon)) %>%
#    mutate(per.reads=round(100*num.reads.taxon/sum(num.reads.taxon), 2))

#krak.top[is.na(krak.top)]<- 0

plot.res=krak.top %>%
    group_by(sample) %>%
    arrange(sample, desc(num.reads.clade)) %>%
##    mutate(per.reads=100*num.reads.clade/sum(num.reads.clade))
        mutate(per.reads=round(100*num.reads.clade/sum(num.reads.clade), 2))

##plot.res[is.na(plot.res)] <- 0
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
pdf(file.path(plotdir, "heat_ill_CARD_grouped_new.pdf"), width=15, height=8.5)

ggplot(plot.res, aes(sample, sci.name)) + geom_tile(aes(fill=per.reads)) +
    scale_fill_gradient(low="blue", high="red") + ylab("AMR genes") + xlab("Samples") #+ theme(panel.background = element_rect(fill='lightgoldenrod'))
#ggplot(plot.res, aes(x=sample, y=per.reads, fill=sci.name))+geom_bar(color="black", stat='identity')+
#    scale_fill_hue()+theme_bw()+theme(legend.position="bottom")

dev.off()


