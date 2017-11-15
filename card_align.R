library(tidyverse)
library(googlesheets)

rawdir="~/Data/vre"
workdir="~/Data/kraken"

dir.create(workdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

##Setup outfilenames
dataloc$nano.card=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(workdir, paste0(dataloc$sample, ".nano.card")))
dataloc$ill.card=file.path(workdir, paste0(dataloc$sample, ".ill.card"))

setwd("/home/ubuntu/")

if (TRUE) {
    ##card minimap

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        system(paste0("/bin/bash -c ", shQuote(paste0("source activate mappers; ",
                                                      "bowtie2 -p 30 --local -x ~/carddb/card_homolog -1 ", dataloc$illumina.r1.aws[i],
                                                      " -2 ", dataloc$illumina.r2.aws[i], " | samtools view -b - | samtools sort - -o ",
                                                      paste0(dataloc$ill.card[i], ".bam")))))
        



        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("/bin/bash -c ", shQuote(paste0("source activate mappers; ",
                                                          "minimap2 -a -x map-ont -t 30 ~/carddb/card_homolog.mmi ",
                                                          dataloc$nanopore.fa.aws[i], " | samtools view -b - | samtools sort -o ",
                                                          paste0(dataloc$nano.card[i], ".bam")))))
        }       
    }

}

if (TRUE) {

    ##ariba

    
    for (i in 1:dim(dataloc)[1]) {                    
        system(paste0("/bin/bash -c ", shQuote(paste0("source activate ariba; ",
                                                      "ariba run ~/ariba-work/outcard.prepareref --threads 30 ",
                                                      dataloc$illumina.r1.aws[i], " ",
                                                      dataloc$illumina.r2.aws[i], " ",
                                                      paste0(dataloc$ill.card[i], ".dir/")))))
                                                             
    }       
}
    
if (TRUE) {
    ##card blast

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        #system(paste0("/bin/bash -c ", shQuote(paste0("source activate mappers; ",
        #"bowtie2 -p 30 --local -x ~/carddb/card_homolog -1 ", dataloc$illumina.r1.aws[i],
        #" -2 ", dataloc$illumina.r2.aws[i], " | samtools view -b - | samtools sort - -o ",
        #paste0(dataloc$ill.card[i], ".bam")))))

        ##Blast reads against CARD
        cat 160223_VRE7_recall_2dhq.fa | parallel --block 100k --recstart '>' --pipe \
        blastn -db CARDblast -query - -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 \
        -outfmt '"6 qseqid sseqid pident length qlen slen evalue bitscore"' > nanoblast2
        
        


        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("/bin/bash -c ", shQuote(paste0("source activate mappers; ",
                                                          "minimap2 -a -x map-ont -t 30 ~/carddb/card_homolog.mmi ",
                                                          dataloc$nanopore.fa.aws[i], " | samtools view -b - | samtools sort -o ",
                                                          paste0(dataloc$nano.card[i], ".bam")))))
        }       
    }

}



#fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
#gs_ws_new(fullsheet, ws_title="VRE_171011", input=dataloc)

#setwd(workdir)

#system("gzip *krak*")

#system("aws s3 sync ~/Data/kraken/ s3://timpawsanalysis/171012_vre/kraken/")
