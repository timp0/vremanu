library(Biostrings)
library(tidyverse)
library(ggExtra)
library(GenomicAlignments)

plotdir="~/Dropbox/Data/Nanopore/161110_vrepaper"

card.dir="/atium/Data/Nanopore/Analysis/161102_vrecontext/CARD"

card.fa=readDNAStringSet(file.path(card.dir, "nucleotide_fasta_protein_homolog_model.fasta"))
card.tab=read_csv(file.path(card.dir, "aro.csv"))

nano.reads=readDNAStringSet(file.path(workdir, "160223_VRE7_recall_2dhq.fa"))

workdir="/atium/Data/Nanopore/Analysis/161102_vrecontext/VRE7_try"


if (TRUE) {

    ##BLAST reads against CARD

    #blast.cols=c("query.id", "subject.id", "per.identity", "align.len", "mismatch", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
    
    blast.cols=c("query.id", "subject.id", "per.identity", "align.len", "query.len", "subject.len", "evalue", "bit.score")
    
    blast.res=tbl_df(read.delim(file.path(workdir, "nanoblast2"), stringsAsFactors=F, header=F, col.names=blast.cols)) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
        mutate(per.sub=align.len/subject.len)
    
    ##Ok - let's plot e value and bit score for each aro
    pdf(file.path(plotdir, "evalue2.pdf"))
    
    ggplot(blast.res, aes(x=aro, y=evalue))+geom_point(alpha=.1)+theme_bw()+scale_y_log10()
    ggplot(blast.res, aes(x=aro, y=bit.score))+geom_point(alpha=.1)+theme_bw()
    
    dev.off()
    
    aro.stat=blast.res %>%
        group_by(aro) %>%
        summarize(eval.min=min(evalue), eval.max=max(evalue), eval.mean=mean(evalue),
                       bit.min=min(bit.score), bit.max=max(bit.score), bit.mean=mean(bit.score))
    
    pdf(file.path(plotdir, "arostats2.pdf"))
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10()
    ggMarginal(z)
    
    z=ggplot(aro.stat, aes(x=eval.min, y=bit.max))+geom_point(alpha=.3)+theme_bw()+scale_x_log10(limits=c(1e-29, 1))
    ggMarginal(z)
    
    dev.off()

    ##ok - filter

    aro.hits=filter(blast.res, evalue==0) %>%
        group_by(aro) %>%
        summarize(n=n(), mean.id=mean(per.identity), min.id=min(per.identity), mean.al=mean(per.sub), bit.mean=mean(bit.score))
    
    pdf(file.path(plotdir, "blast_aro_hits.pdf"))
    
    z=ggplot(aro.hits, aes(x=mean.id, y=mean.al))+geom_point(alpha=.3)+theme_bw()
    ggMarginal(z)

    dev.off()
    
    pdf(file.path(plotdir, "blast_aro_hits_id.pdf"))
    
    ggplot(blast.res, aes(x=as.numeric(factor(aro)), y=per.identity))+geom_density_2d()+theme_bw()
    ggplot(blast.res, aes(x=as.numeric(factor(aro)), y=per.sub))+geom_density_2d()+theme_bw()
    ggplot(blast.res, aes(x=as.numeric(factor(aro)), y=bit.score))+geom_density_2d()+theme_bw()
    
    dev.off()
    
    ##Still lots of false positives?

}
    

if (FALSE) {

    ##LAST load
    
    last.tsv=c("score", "subject.id", "s.start", "align.len", "s.strand", "s.length", "query.id", "q.start", "align.len_2", "q.strand", "q.len", "aln", "mismap")
    
    last.res=read_tsv(file.path(workdir, "lastcard.tsv"), comment="#", col_names=last.tsv) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
        mutate(per.sub=align.len/s.length)
    
    last.blast=c("query.id", "subject.id", "per.identity", "align.len", "mismatch", "gap.opens", "q.start", "q.end", "s.start", "s.end")
    
    last.blast.res=read_tsv(file.path(workdir, "lastblast.tsv"), col_names=last.blast) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) 
        
    
    last.hits=last.res %>%
        group_by(aro) %>%
        summarize(n=n(), score=mean(score), per.sub=mean(per.sub)) %>%
        mutate(blast.match=(aro %in% aro.hits$aro))
    
    
    pdf(file.path(plotdir, "last_aro_hits.pdf"))
        
    ggplot(last.hits, aes(x=score, y=per.sub, color=blast.match, size=n))+geom_point(alpha=.3)+theme_bw()    
    
    dev.off()

    last.blast.hits=last.blast.res %>%
        group_by(aro) %>%
        summarize(n=n(), per.id=mean(per.identity), align.len=mean(align.len)) %>%
        mutate(blast.match=(aro %in% aro.hits$aro))

    pdf(file.path(plotdir, "last_blast_aro_hits.pdf"))

    ggplot(last.blast.hits, aes(x=per.id, y=align.len, color=blast.match, size=n))+geom_point(alpha=.3)+theme_bw()
   

    dev.off()

    
}


if (TRUE) {
    ##Ok - trying bwa alignment of nanopore to CARD.  seems ok to me to do this.

    ##Load in BAM.  Plot reads per gene aligned, mapq score distro

    tags=ScanBamParam(tag="NM", what=c("flag", "qname"))
    
    bwa.card=tbl_df(readGAlignments(file.path(workdir, "vre7_card.sorted.bam"), param=tags)) %>%
        mutate(subject.id=as.character(seqnames)) %>%
        select(subject.id, qwidth, width, qname, NM) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
        mutate(s.length=width(card.fa[pmatch(subject.id, names(card.fa), duplicates.ok=T)])) %>%
        mutate(per.sub=width/s.length, per.identity=1-(NM/width))

    bwa.hits=bwa.card %>%
        group_by(aro) %>%
        summarize(n=n(), per.id=mean(per.identity), per.sub=mean(per.sub)) %>%
        mutate(blast.match=factor((aro %in% aro.hits$aro)+(aro %in% last.hits$aro)))
    
    
    pdf(file.path(plotdir, "bwa_aro_hits.pdf"))
        
    ggplot(bwa.hits, aes(x=per.id, y=per.sub, color=blast.match, size=n))+geom_point(alpha=.3)+theme_bw()    
    
    dev.off()    
    
}

##Try last bam file

if (TRUE) {

    tags=ScanBamParam(tag="NM", what=c("flag", "qname"))
    
    last.card=tbl_df(readGAlignments(file.path(workdir, "last.sorted.bam"), param=tags))%>%
        mutate(subject.id=as.character(seqnames)) %>%
        select(subject.id, qwidth, width, qname, NM) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
        mutate(s.length=width(card.fa[pmatch(subject.id, names(card.fa), duplicates.ok=T)])) %>%
        mutate(per.sub=width/s.length, per.identity=1-(NM/width))

    last.bam.hits=last.card %>%
        group_by(aro) %>%
        summarize(n=n(), per.id=mean(per.identity), per.sub=mean(per.sub)) %>%
        mutate(blast.match=factor((aro %in% aro.hits$aro)+(aro %in% bwa.hits$aro)))
    
    
    pdf(file.path(plotdir, "last_bam_aro_hits.pdf"))
        
    ggplot(last.bam.hits, aes(x=per.id, y=per.sub, color=blast.match, size=n))+geom_point(alpha=.3)+theme_bw()    
    
    dev.off()    
 
}

##Output tables 
if (TRUE) {

    ##Table document
    aro.hits=aro.hits %>%
        mutate(last.match=aro %in% last.bam.hits$aro) %>%
        mutate(bwa.match=aro %in% bwa.hits$aro) %>%
        mutate(aro.idx=pmatch(aro, card.tab$Accession)) %>%
        mutate(aro.name=card.tab$Name[aro.idx]) %>%
        mutate(aro.desc=card.tab$Description[aro.idx])

    write_csv(aro.hits, path=file.path(plotdir, "blastres.csv"))

    ##Table document
    last.bam.hits=last.bam.hits %>%
        mutate(blast.match=aro %in% aro.hits$aro) %>%
        mutate(bwa.match=aro %in% bwa.hits$aro) %>%
        mutate(aro.idx=pmatch(aro, card.tab$Accession)) %>%
        mutate(aro.name=card.tab$Name[aro.idx]) %>%
        mutate(aro.desc=card.tab$Description[aro.idx])

    write_csv(last.bam.hits, path=file.path(plotdir, "lastres.csv"))

    ##Table document
    bwa.hits=bwa.hits %>%
        mutate(blast.match=aro %in% bwa.hits$aro) %>%
        mutate(last.match=aro %in% last.bam.hits$aro) %>%
        mutate(aro.idx=pmatch(aro, card.tab$Accession)) %>%
        mutate(aro.name=card.tab$Name[aro.idx]) %>%
        mutate(aro.desc=card.tab$Description[aro.idx])

    write_csv(bwa.hits, path=file.path(plotdir, "bwares.csv"))

}    
