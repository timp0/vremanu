library(tidyverse)
library(googlesheets)
library(cowplot)
library(jsonlite)

##aws s3 sync  s3://timpawsanalysis/171012_vre/kraken/ ~/Data/kraken/
workdir="~/Data/kraken"
plotdir="~/Dropbox/timplab_data/vremanu/180513_vre7"
dir.create(plotdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

##Ok, basically blast for KPC
##There are lots of KPC variants
##Let's just try to blast all genes first and see what that gets us

##I think we just need protein_homolog model
##First make blastdb

cardhomdb="/mithril/Data/NGS/Reference/card/201/nuchomolog"

vre7dat=dataloc %>%
    filter(sample=="VRE7")


if (FALSE) {

   system(paste0("makeblastdb -in /mithril/Data/NGS/Reference/card/201/nucleotide_fasta_protein_homolog_model.fasta -out ", cardhomdb, " -dbtype nucl -input_type fasta"))
    
}

##Now blast

if (FALSE) {

    system(paste0("blastn -db ", cardhomdb, " -query ", vre7dat$nanopore.fasta, " -out ", file.path(plotdir, "cardblast2"),
                  " -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen evalue bitscore\" -num_threads=8"))
}


##Extract times for VRE7
if (FALSE) {
    system(paste0("~/vremanu/time_extract.py -i ", vre7dat$nanopore.raw, " -o ", file.path(plotdir, "timetest.csv.gz")))
}

if (TRUE) {

    ##Split times readname by spa ce, keep just first one
    times=read_csv(file.path(plotdir, "timetest.csv.gz")) %>%
        separate(readname, into=c("nothing", "id"), sep="[@ ]", extra="drop") %>%
        select(-nothing) %>%
        rename(read.id=id)


    
}

if (TRUE) {

    ##Card parents from florian
    parent=read_tsv(file.path(plotdir, 'card-parents.tsv')) %>%
        rename(aro=ID)
    
}


if (TRUE) {
        
    blast.cols=c("query.id", "subject.id", "per.identity", "query.len", "subject.len", "align.len",
                 "query.start", "query.end", "subject.start", "subject.end",
                 "mismatch", "gap.opens", "evalue", "bit.score")

    ##load in data
    blastres=read_tsv(file.path(plotdir, "cardblast2"), col_names=blast.cols)

    blastres=blastres %>%
        mutate(subcov=align.len/subject.len) %>%
        filter((per.identity>80) & (subcov>.5))
    
    blastres= blastres %>%
        separate(subject.id, sep="\\|", into=c("gb", "ncbi", "strand", "coord", "aro", "name"), remove=F) %>%
        select(-gb, -ncbi, -strand, -coord) %>%
        rename(read.id=query.id)

    blastres=left_join(blastres, parent)
    
}

##Now abricate

if (FALSE) {

    system(paste0("/bin/bash -c ", shQuote(paste0("source activate tseemann; abricate --db card ", vre7dat$nanopore.fasta, " >",
                                                  file.path(plotdir, "abricate7.tsv")))))
}

if (FALSE) {

    abr.card=read_tsv(file.path(plotdir, "abricate7.tsv")) %>%
        rename(per.identity=`%IDENTITY`, per.coverage=`%COVERAGE`) %>%
        filter((per.identity>80) & (per.coverage>50)) %>%
        select(-COVERAGE_MAP, -DATABASE) %>%
        rename(Name=GENE)

    ##Not working as well as blast tbh
    z=left_join(abr.card, parent)    
    
}




if (FALSE) {
   
    readhits=blastres %>%
        group_by(query.id) %>%
        slice(1)
    
    kpchits=readhits %>%
        filter(grepl("KPC", subject.id)) %>%
        arrange(-query.len)

    favkpc=kpchits %>%
        filter(grepl("KPC-1$", subject.id)) %>%
        arrange(-query.len)

    write(favkpc$query.id, file=file.path(plotdir, "kpcreads.txt"))
    
    ##Get out just those hits

    kpc.gff=favkpc %>%
        mutate(source="cardblast", feature="CDS", frame="0") %>%
        mutate(strand=ifelse(subject.end>subject.start, yes="+", no="-")) %>%
        mutate(attrib=paste0("Name=KPC-1;gene=KPC-1;per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue)) %>%
        select(seqname=query.id, source=source, feature=feature, start=query.start, end=query.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)

    write_tsv(kpc.gff, file.path(plotdir, "kpcblast.gff"), col_names=F)
    
    system(paste0("seqtk subseq ", vre7dat$nanopore.fasta, " ", file.path(plotdir, "kpcreads.txt"),
           " >", file.path(plotdir, "vre7kpc.fa")))
            
}

if (TRUE) {

    perread=blastres %>%
        group_by(read.id) %>%
        slice(1)           
    
    perhit=perread %>%
        group_by(Parents) %>%
        summarise(hits=n()) %>%
        arrange(-hits)
   
    write_csv(perhit, file.path(plotdir, "perhit2.csv"))

    perread=left_join(perread, times) %>%
        mutate(minutes=start.time/60)
    
    myplots=perread %>%
        group_by(Parents) %>%
        filter(n() > 10) %>% ##Remove those with less than 10 hits
        arrange(start.time) %>%
        mutate(cumread=1:length(start.time)) %>%
        do(
            plots=ggplot(data=.) +
                geom_step(aes(x=minutes, y=cumread))+
                ggsave(file.path(plotdir, "hits", paste0(.$Parents[1], ".pdf")))
        )
           
}
