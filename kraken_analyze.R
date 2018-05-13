library(tidyverse)
library(googlesheets)
#library(phyloseq)
library(vegan)


workdir="~/Data/kraken"
plotdir="~/Dropbox/Data/Nanopore/180303_vre"
dir.create(plotdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

read_kraken_report <- function(samp.list) {
    ##column names
    krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")
    
    ##Load report data and combine with sample name embedded
    krak.report=bind_rows(lapply(samp.list$filey, function(x) {read_tsv(x, col_names=krak.report.cnames) %>% mutate(file=x)})) %>%
        mutate(sample=samp.list$sample[pmatch(file, samp.list$filey, duplicates.ok=T)])
    
    ##Trim left side white space in sci.name
    krak.report$sci.name=trimws(krak.report$sci.name)

    return(krak.report)
}

filter_krak_report <- function(krak.report, numtop=2) {
    
    ##ok - for first pass, let's try filtering by "species" or S in rank.code
    ##per.reads is of total reads, not the subsection of bacteria or just classified, as opposed to how pavian does it
    ##We also want .clade reads because that takes all the subspecies/strains - otherwise it will undercall because of subspecies
    
    ##Sort for top 2 not human
    krak.top=filter(krak.report, (rank.code=="S" & sci.name!="Homo sapiens")) %>%
        group_by(sample) %>%
        top_n(n=numtop, wt=num.reads.clade) %>%
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
    
    
    ##Finally our actual hits of interest
    plot.res=bind_rows(plot.res,
                       filter(krak.report, ncbi.tax.id %in% krak.top$ncbi.tax.id) %>%
                       select(sample, sci.name, num.reads.clade))
    
    ##Sort by rank order, also generate % reads per samp
    plot.res=plot.res %>%
        group_by(sample) %>%
        arrange(sample, desc(num.reads.clade)) %>%
        mutate(per.reads=round(100*num.reads.clade/sum(num.reads.clade), 2))


    ##fix x-axis order so VRE10 on the end
    plot.res$sample=factor(plot.res$sample, levels=samp.list$sample)

    
    ##Fix y-axis order
    yorder=c("unclassified", "Other bacteria", "Homo sapiens", "Archaea", "Viruses", krak.top$sci.name)
    plot.res$sci.name=factor(plot.res$sci.name, levels=yorder)

   
    return(plot.res)
    
}



plot_stacked <- function(plot.res, plotdir, namey) {
    krak.report=read_kraken_report(samp.list)    
    
    plot.res=filter_krak_report(krak.report)    

    plot.res=plot.res %>%
        arrange(sample, sci.name)
   
    ##Italize bacteria(?)
    pdf(file.path(plotdir, paste0(namey, '.pdf')), width=11, height=8.5)
    
    print(ggplot(plot.res, aes(x=sample, y=per.reads, fill=sci.name))+geom_bar(color="black", stat='identity')+
        scale_fill_hue()+theme_bw()+theme(legend.position="bottom"))
    
    dev.off()

    write_csv(plot.res, file.path(plotdir, paste0(namey, '.csv')))
    
    
}

##system(paste0("gunzip ~/Data/kraken/*report.gz"))

if (TRUE) {
    ##Ill first
    
    samp.list=dataloc %>%
        select(sample, k31.ill.aws) %>%
        rename(krak.file=k31.ill.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))

    plot_stacked(samp.list, plotdir, 'stack_k31_ill1b')

    krak.report=read_kraken_report(samp.list)    
    
    plot.res=filter_krak_report(krak.report, numtop=10)    
    
    try=plot.res %>%
        select(sample, sci.name, num.reads.clade) %>%
        spread(sci.name, num.reads.clade, fill=0)
            
    z=as.matrix(try[,-1])

    vegdist(z)
    diversity(z)

    ##Ok - figure out ways to plot alpha diversity, bray-curtis similarity, clustering to see if things from nanopore
    ##And illumina show same result?

    ##Figure out comparsion besides plot

    ##Should I bother trying to use centrifuge results?  Not clear
    
    ##Bracken
    samp.list=dataloc %>%
        select(sample, k31.ill.aws) %>%
        rename(krak.file=k31.ill.aws) %>%
        mutate(filey=paste0(krak.file, "_bracken.report"))
    brack.report=read_kraken_report(samp.list)

    plot.res=filter_krak_report(brack.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k31_ill_bracken1a')

    ##Ill 24 k-mer
    samp.list=dataloc %>%
        select(sample, k24.ill.aws) %>%
        rename(krak.file=k24.ill.aws) %>%
        mutate(filey=paste0(krak.file, "_bracken.report"))
    brack.report=read_kraken_report(samp.list)
    
    plot.res=filter_krak_report(brack.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k24_ill_bracken1a')

    samp.list=dataloc %>%
        select(sample, k24.ill.aws) %>%
        rename(krak.file=k24.ill.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    krak.report=read_kraken_report(samp.list)
    
    plot.res=filter_krak_report(krak.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k24_ill1')


    ##Nanopore    
    samp.list=dataloc %>%
        filter(!is.na(k31.nano.aws)) %>%
        select(sample, k31.nano.aws) %>%
        rename(krak.file=k31.nano.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    krak.report=read_kraken_report(samp.list)

    plot.res=filter_krak_report(krak.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k31_nano1')


    samp.list=dataloc %>%
        filter(!is.na(k31.nano.aws)) %>%
        select(sample, k31.nano.aws) %>%
        rename(krak.file=k31.nano.aws) %>%
        mutate(filey=paste0(krak.file, "_bracken.report"))
    krak.report=read_kraken_report(samp.list)

    plot.res=filter_krak_report(krak.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k31_bracken_nano1')


    
    ##K24 nanopore
    samp.list=dataloc %>%
        filter(!is.na(k24.nano.aws)) %>%
        select(sample, k24.nano.aws) %>%
        rename(krak.file=k24.nano.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    krak.report=read_kraken_report(samp.list)

    plot.res=filter_krak_report(krak.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k24_nano1')

    ##Rerun bracken on nanopore k24
    ##K24 nanopore
    samp.list=dataloc %>%
        filter(!is.na(k24.nano.aws)) %>%
        select(sample, k24.nano.aws) %>%
        rename(krak.file=k24.nano.aws) %>%
        mutate(filey=paste0(krak.file, "_bracken.report"))
    krak.report=read_kraken_report(samp.list)

    plot.res=filter_krak_report(krak.report)
    
    plot_stacked(plot.res, plotdir, 'stack_k24_bracken_nano1')

    
}

