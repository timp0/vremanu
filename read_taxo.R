library(tidyverse)
library(cowplot)
library(googlesheets)
library(animation)


workdir="~/Data/kraken"
plotdir="~/Dropbox/timplab_data/vremanu/180524_vre"
dir.create(plotdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

dataloc=dataloc %>%
    filter(!is.na(nanopore.fasta))

##Extract times
if (TRUE) {
    for (i in 1:nrow(dataloc)) {
        system(paste0("~/vremanu/time_extract.py -i ", dataloc$nanopore.raw[i], " -o ", file.path(plotdir, paste0(dataloc$sample[i], ".csv.gz"))))
    }
}

##Also need to get kraken reports loaded and matched up 
if (TRUE) {

    keeptop=5
    
    i=5
    
    times=read_csv(file.path(plotdir, paste0(dataloc$sample[i], ".csv.gz"))) %>%
        separate(readname, into=c("nothing", "id"), sep="[@ ]", extra="drop") %>%
        select(-nothing) %>%
        rename(read.id=id)

    
    krak.cols=c("class.state", "read.id", "tax.id", "seq.length", "lca.map")
    krak=read_tsv(dataloc$k31.nano.aws[i], col_names=krak.cols) %>%
        filter(class.state=="C") ##Only care about classified reads

    krak=left_join(krak, times)

       
    ##Load report
    krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "tax.id", "sci.name")
    krak.report=read_tsv(paste0(dataloc$k31.nano.aws[i], ".report"), col_name=krak.report.cnames) %>%
        mutate(sci.name=trimws(sci.name))

    krak.top=filter(krak.report, (rank.code=="S")) %>%
        top_n(n=keeptop, wt=num.reads.clade) %>%
        ungroup %>%
        select(sci.name, tax.id) %>%
        distinct                          

    krak.clade=krak.report %>%
        select(-per.reads, -num.reads.clade, -num.reads.taxon) %>%
        mutate(subidx=cumsum(rank.code!="-")) %>%
        group_by(subidx) %>%
        mutate(clade.name=sci.name[1]) %>%
        ungroup()     
    
    krak=left_join(krak, krak.clade) %>%
        select(-lca.map, -subidx, -class.state)

    top.clades=krak.report %>%
        filter(rank.code=="S") %>%
        arrange(-num.reads.clade) %>%
        slice(1:keeptop)

    krak$clade.name=ifelse(krak$clade.name %in% top.clades$sci.name, krak$clade.name, "Other")
    
    krak=krak %>%
        mutate(minute=start.time/60) %>%
        arrange(start.time) %>%
        mutate(time.block=floor(start.time/600))

    ##Let's plot!

    krak=krak %>%
        group_by(clade.name) %>%
        arrange(start.time) %>%
        mutate(cumread=1:n()) %>%
        ungroup()

    
    pdf(file.path(plotdir, "tax_time.pdf"), width=11, height=8.5)
    print(ggplot(krak, aes(x=minute, y=cumread, color=clade.name))+geom_step()+scale_x_continuous()+theme_bw()+xlim(0,1000))
    dev.off()
        
    pdf(file.path(plotdir, "tax_pie.pdf"))
    print(ggplot(krak, aes(x=factor(1), fill=clade.name))+geom_bar(color="black")+coord_polar(theta="y"))
    dev.off()
    
    tstep=100
    tmax=1000

    ani.options('autobrowse'=FALSE)
    
    saveGIF({
        for (i in seq(from=1, to=1001, by=100)) {
            temp.dat=krak %>%
                filter(minute<i)
            print(ggplot(temp.dat, aes(x=factor(1), fill=clade.name))+geom_bar(color="black")+coord_polar(theta="y"))
        }}, movie.name=file.path(plotdir, "pietry.gif"))

    
}
    

##Original VRE6 pie_chart movie
if (FALSE) {

    timefile="/mithril/Data/Nanopore/Analysis/160530_newVRE/160526_VRE6.called_2dhq.tout.csv.gz"
    
    vre6.times=read_csv(gzfile(timefile)) %>%
        mutate(read.id=gsub("@", "", sapply(strsplit(readname, split=" "), function(x) {x[1]})))
    
    krakfile="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/VRE6_nanopore.kraken"
    
    krak.cols=c("class.state", "seq.id", "tax.id", "seq.length", "lca.map")
    
    vre6.krak=read_tsv(krakfile, col_names=krak.cols)
    
    ##column names
    krak.report.cnames=c("per.reads", "num.reads.clade", "num.reads.taxon", "rank.code", "ncbi.tax.id", "sci.name")
    
    reportfile="/atium/Data/NGS/Aligned/161110_VREkraken/Nanopore_kraken/report/VRE6_nanopore.kraken.output"
    
    ##Load report data and combine with sample name embedded
    krak.report=read_tsv(reportfile, col_name=krak.report.cnames) %>%
        mutate(sci.name=trimws(sci.name))
    
    ##Get top 5
                                        #krak.top=filter(krak.report, (rank.code=="S" & sci.name!="Homo sapiens")) %>%
    krak.top=filter(krak.report, (rank.code=="S")) %>%
        top_n(n=5, wt=num.reads.clade) %>%
        ungroup %>%
        select(sci.name, ncbi.tax.id) %>%
        distinct                          
    
    ##Ok - pseudocode because too much talking
    ##Find row index of the top 5
    ##Check followin indexes for *non -* - apparently - indicates a sub species hit
    ##Get those ncbi.taxi.ids
    ##bundle them under the global sci.name
    
    krak.top.clade=krak.report %>%
        mutate(subidx=cumsum(rank.code!="-")) %>%
        group_by(subidx) %>%
        mutate(clade.name=sci.name[1]) %>%
        ungroup() %>%
        filter(clade.name %in% krak.top$sci.name)
    
    ##filter those reads
    vre6.krak.filt=vre6.krak %>%
        filter(tax.id %in% krak.top.clade$ncbi.tax.id) %>%
        mutate(tax.name=krak.top.clade$clade.name[pmatch(tax.id, krak.top.clade$ncbi.tax.id, duplicates.ok=T)])
    
    ##Get times for those reads
    vre6.krak.filt=vre6.krak.filt %>%
        mutate(time=vre6.times$end.time[pmatch(seq.id, vre6.times$read.id)])
    
    
    ##Let's plot!
    
    plot.dat=vre6.krak.filt %>%
        group_by(tax.name) %>%
        arrange(time) %>%
        mutate(cumread=1:n()) %>%
        ungroup() %>%
        mutate(minute=time/60)
    
    plotdir="~/Dropbox/Data/Nanopore/161110_vrepaper"
    
    pdf(file.path(plotdir, "tax_time.pdf"))
    print(ggplot(plot.dat, aes(x=minute, y=cumread, color=tax.name))+geom_step()+scale_x_continuous()+theme_bw()+xlim(0,1000))
    dev.off()
    
    
    pdf(file.path(plotdir, "tax_pie.pdf"))
    print(ggplot(plot.dat, aes(x=factor(1), fill=tax.name))+geom_bar(color="black")+coord_polar(theta="y"))
    dev.off()
    
    ##ok - let's get dangerous and try to animate
    
    library(animation)
    theme_set(theme_bw())
    
    tstep=100
                                        #tmax=max(plot.dat$minute)
    tmax=1000
    
    saveVideo({
        for (i in seq(from=1, to=1001, by=100)) {
            temp.dat=plot.dat %>%
                filter(minute<i)
            print(ggplot(temp.dat, aes(x=factor(1), fill=tax.name))+geom_bar(color="black")+coord_polar(theta="y"))        
        }}, movie.name = "/home/timp/Dropbox/Temp/pietry.mp4")
    
}

    
    
