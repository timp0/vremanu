library(tidyverse)

plot_stacked <- function(samp.list, plotdir, namey) {

    ##column names
    krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")

    ##Load report data and combine with sample name embedded
    krak.report=bind_rows(lapply(samp.list$filey, function(x) {read_tsv(x, col_names=krak.report.cnames) %>% mutate(file=x)})) %>%
        mutate(sample=samp.list$label[pmatch(file, samp.list$filey, duplicates.ok=T)])

    ##Trim left side white space in sci.name
    krak.report$sci.name=trimws(krak.report$sci.name)
    
    ##ok - for first pass, let's try filtering by "species" or S in rank.code
    ##per.reads is of total reads, not the subsection of bacteria or just classified, as opposed to how pavian does it
    ##We also want .clade reads because that takes all the subspecies/strains - otherwise it will undercall because of subspecies
    
    ##Sort for top 2
    krak.top=filter(krak.report, (rank.code=="S" & sci.name!="Homo sapiens")) %>%
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
    plot.res=bind_rows(plot.res,
                       filter(krak.report, sci.name=="Homo sapiens") %>%
                       select(sample, sci.name, num.reads.clade))
    
    
    ##per group non bac/non human hits
    ##Note - we assume no eukaryotes (2759) in kraken report but human, because that's what kraken database was
    plot.res=bind_rows(plot.res,
                       filter(krak.report, rank.code=="D", ncbi.tax.id!=2, ncbi.tax.id!=2759) %>%
                       select(sample, sci.name, num.reads.clade))
    
    ##per group unclass hits
    ##Unclass just means save U
    #plot.res=bind_rows(plot.res,
    #                   filter(krak.report, rank.code=="U") %>%
    #                   select(sample, sci.name, num.reads.clade))
    
    
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
    plot.res$sample=factor(plot.res$sample, levels=samp.list$label)
    
    ##Fix y-axis order
    yorder=c("unclassified", "Other bacteria", "Homo sapiens", "Archaea", "Viruses", krak.top$sci.name)
    plot.res$sci.name=factor(plot.res$sci.name, levels=yorder)

    plot.res=plot.res %>%
        arrange(sample, sci.name)

    
    ##Better color scheme still
    ##Human reads
    ##reorder so that blocks are all sorted in same order - some order we define, like unclass, other bac, archae, virus,
    ##then our bac of interest

    ##Italize bacteria
    pdf(file.path(plotdir, paste0(namey, '.pdf')), width=11, height=8.5)
    
    print(ggplot(plot.res, aes(x=sample, y=per.reads, fill=sci.name))+geom_bar(color="black", stat='identity')+
        scale_fill_hue()+theme_bw()+theme(legend.position="bottom"))
    
    dev.off()

    write_csv(plot.res, file.path(plotdir, paste0(namey, '.csv')))
    
    
}



##Illumina kraken files:
ill.krak="/atium/Data/NGS/Aligned/161110_VREkraken/VRE_Illumina_withhuman"
plotdir="~/Dropbox/Data/Nanopore/161110_vrepaper"

samp.list=tibble(filey=file.path(ill.krak, paste0("VRE", 1:10, "_withhuman.kraken.output")),
                 label=paste0("VRE", 1:10, "Illumina"))

plot_stacked(samp.list, plotdir, 'stack_ill1')


nano.krak="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/report"

samp.list=tibble(filey=file.path(nano.krak, paste0("VRE", c("NC", 3, 5, 6, 7, 10), "_nanopore.kraken.output")),
                 label=paste0("VRE", c("NC", 3, 5, 6, 7, 10), "Nanopore"))

plot_stacked(samp.list, plotdir, 'stack_nano1')

