library(tidyverse)
library(googlesheets)
library(phyloseq)

workdir="~/Data/kraken"
plotdir="~/Data/plots"
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



plot_stacked <- function(samp.list, plotdir, namey) {

    krak.report=read_kraken_report(samp.list)
    
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
    plot.res$sample=factor(plot.res$sample, levels=samp.list$sample)
    
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

##system(paste0("gunzip ~/Data/kraken/*report.gz"))

if (TRUE) {
    ##Ill first
    
    samp.list=dataloc %>%
        select(sample, k31.ill.aws) %>%
        rename(krak.file=k31.ill.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    
    
    plot_stacked(samp.list, plotdir, 'stack_k31_ill1')
    
    samp.list=dataloc %>%
        select(sample, k24.ill.aws) %>%
        rename(krak.file=k24.ill.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    
    plot_stacked(samp.list, plotdir, 'stack_k24_ill1')
    
    samp.list=dataloc %>%
        filter(!is.na(k31.nano.aws)) %>%
        select(sample, k31.nano.aws) %>%
        rename(krak.file=k31.nano.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    
    plot_stacked(samp.list, plotdir, 'stack_k31_nano1')
    
    samp.list=dataloc %>%
        filter(!is.na(k24.nano.aws)) %>%
        select(sample, k24.nano.aws) %>%
        rename(krak.file=k24.nano.aws) %>%
        mutate(filey=paste0(krak.file, ".report"))
    
    plot_stacked(samp.list, plotdir, 'stack_k24_nano1')
    
}

if (TRUE) {



    ##take first 10 for now
    
    kraken.report.list=paste(paste0(dataloc$k31.ill.aws[1:10], ".report"), collapse=" ")

    system(paste0("/bin/bash -c ", shQuote(paste0("source activate qiime2-2017.9; ",
                                                  "kraken-biom ", kraken.report.list, " -o ~/Data/kraken/try.biom ",
                                                  "--otu_fp ~/Data/kraken/try.outfp"))))

    system(paste0("/bin/bash -c ", shQuote(paste0("source activate qiime2-2017.9; ",
                                                  "qiime tools import --input-path ~/Data/kraken/try.biom ",
                                                  "--type 'FeatureTable[Frequency]' ",
                                                  "--source-format BIOMV210Format ",
                                                  "--output-path try.qza"))))

    

    z=import_biom("~/Data/kraken/try.biom")

}
