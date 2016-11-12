library(Biostrings)
library(ggplot2)
library(dplyr)
library(ggExtra)


plotdir="~/Dropbox/Data/Nanopore/161110_vrepaper"

blast.cols=c("query.id", "subject.id", "per.identity", "align.len", "mismatch", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")

card.dir="/atium/Data/Nanopore/Analysis/161102_vrecontext/CARD"

card.fa=readDNAStringSet(file.path(card.dir, "nucleotide_fasta_protein_homolog_model.fasta"))
card.tab=read.csv(file.path(card.dir, "aro.csv"), stringsAsFactors=F)

workdir="/atium/Data/Nanopore/Analysis/161102_vrecontext/VRE7db"

blast.res=read.delim(file.path(workdir, "blastry2"), stringsAsFactors=F, header=F, col.names=blast.cols)

nano.reads=readDNAStringSet(file.path(workdir, "160223_VRE7_recall_2dhq.fa"))

blast.res$aro=sapply(strsplit(blast.res$query.id, split="\\|"), function(x) {x[5]})

blast.tbl=tbl_df(blast.res)

if (FALSE) {

    ##Ok - let's plot e value and bit score for each aro
    pdf(file.path(plotdir, "evalue.pdf"))
    
    ggplot(blast.tbl, aes(x=aro, y=evalue))+geom_point(alpha=.1)+theme_bw()+scale_y_log10()
    ggplot(blast.tbl, aes(x=aro, y=bit.score))+geom_point(alpha=.1)+theme_bw()
    
    dev.off()
    
    blast.grp=group_by(blast.tbl, aro)
    
    aro.stat=summarize(blast.grp, eval.min=min(evalue), eval.max=max(evalue), eval.mean=mean(evalue),
                       bit.min=min(bit.score), bit.max=max(bit.score), bit.mean=mean(bit.score))
    
    pdf(file.path(plotdir, "arostats.pdf"))
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10()
    ggMarginal(z)
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10(limits=c(1e-29, 1))
    ggMarginal(z)
    
    dev.off()

}
    
    ##By inspection, lots of values about 1e-3 - asssume those are BS, throw them out

if (TRUE) {

    blast.tbl=filter(blast.tbl, evalue<1e-3)
    blast.grp=group_by(blast.tbl, aro)
    
    ##ok - from 2102 possible, now only have 659 hits, replot

    pdf(file.path(plotdir, "evalue_filt.pdf"))
    
    ggplot(blast.tbl, aes(x=aro, y=evalue))+geom_point(alpha=.1)+theme_bw()+scale_y_log10()
    ggplot(blast.tbl, aes(x=aro, y=bit.score))+geom_point(alpha=.1)+theme_bw()
    
    dev.off()
    
    aro.stat=summarize(blast.grp, eval.min=min(evalue), eval.max=max(evalue), eval.mean=mean(evalue),
                       bit.min=min(bit.score), bit.max=max(bit.score), bit.mean=mean(bit.score))
    
    pdf(file.path(plotdir, "arostats_filt.pdf"))
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10()
    ggMarginal(z)
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10(limits=c(1e-29, 1))
    ggMarginal(z)
    
    dev.off()

    ##Still a lot of hits that may or may not be real
    
}    


gene.hits=unique(blast.tbl$query.id)
gene.fa=card.fa[pmatch(gene.hits, names(card.fa))]

#this.gene=filter(blast.tbl, query.id==gene.hits[i])

setwd(workdir)

i=1



writeXStringSet(gene.fa[i], "temp.fa", append=FALSE, format="fasta")

system(paste("blastn", "-db", "VRE7blast" , "-query", "temp.fa", "-gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue .01 -outfmt '7 std qseq sseq stitle' -out", "temp.blast.tsv"))

system(paste0("~/Code/mview/bin/mview -in blast ", "temp.blast.tsv", " -html head -coloring identity -moltype dna >",
              file.path(plotdir, "temp.html")))
#system(paste0("~/Code/mview/bin/mview -in blast ", temp.blast, " -out fasta >", temp.align.fasta))



#y=pmatch(this.gene$subject.id, names(nano.reads), duplicates.ok=T)


##Need a alignment of some kind of blast hits to look at and see if I believe them/they are significant.  Binned for the two peaks I see in bit score - try <500, <1000, >1000 as bins
    
    
                                        #z=pmatch(x=blast.res$aro, table=card.tab$Accession, duplicates.ok=T)




