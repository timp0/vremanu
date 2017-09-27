library(tidyverse)
library(googlesheets)

workdir="~/Data/vre"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_0905")

dataloc$illumina.r1.aws=file.path(workdir, basename(dataloc$illumina.r1))
dataloc$illumina.r2.aws=file.path(workdir, basename(dataloc$illumina.r2))

dataloc$nanopore.fa.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(workdir, basename(dataloc$nanopore.fasta)))

if (TRUE) {
    ##kraken k-mer 31
    system(paste0("aws s3 sync --delete s3://timpawsanalysis/170913_krakendb.masked.fullstd/ ~/krakendb/"))

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ --paired ", dataloc$illumina.r1.aws[i], " ", dataloc$illumina.r2.aws[i], " >", file.path(workdir, paste0(dataloc$sample[i], ".ill.k31.krak"))))
        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ ", dataloc$nanopore.fa.aws[i], " >", file.path(workdir, paste0(dataloc$sample[i], ".nano.k31.krak"))))
        }       
    }

}

if (TRUE) {
    ##kraken k-mer 24
    system(paste0("aws s3 sync --delete s3://timpawsanalysis/170913_krakendb.masked.fullk24/ ~/krakendb/"))

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ --paired ", dataloc$illumina.r1.aws[i], " ", dataloc$illumina.r2.aws[i], " >", file.path(workdir, paste0(dataloc$sample[i], ".ill.k24.krak"))))
        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ ", dataloc$nanopore.fa.aws[i], " >", file.path(workdir, paste0(dataloc$sample[i], ".nano.k24.krak"))))
        }       
    }
       
   
}


