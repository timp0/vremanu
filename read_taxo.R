library(tidyverse)
library(ggplot2)

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
    for (i in tstep:tstep:tmax) {
        temp.dat=plot.dat %>%
            filter(minute<i)
        print(ggplot(temp.dat, aes(x=factor(1), fill=tax.name))+geom_bar(color="black")+coord_polar(theta="y"))        
    }}, movie.name = "/home/timp/Dropbox/Temp/pietry.mp4")
   

