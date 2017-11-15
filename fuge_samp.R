library(tidyverse)
library(googlesheets)

rawdir="~/Data/vre"
workdir="~/Data/fuge"

dir.create(workdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

dataloc$illumina.r1.aws=file.path(rawdir, basename(dataloc$illumina.r1))
dataloc$illumina.r2.aws=file.path(rawdir, basename(dataloc$illumina.r2))

dataloc$nanopore.fa.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(rawdir, basename(dataloc$nanopore.fasta)))


##Setup outfilenames
dataloc$fuge.nano.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(workdir, paste0(dataloc$sample, ".nano.fuge")))
dataloc$fuge.ill.aws=file.path(workdir, paste0(dataloc$sample, ".ill.fuge"))

setwd("/home/ubuntu/")

for (i in 1:dim(dataloc)[1]) {                    
    ##illumina first
    system(paste0("centrifuge -p 16 -x ~/centrifuge/indices/p+h+v -1 ",dataloc$illumina.r1.aws[i], " -2 ",
                  dataloc$illumina.r2.aws[i], " -S ", dataloc$fuge.ill.aws[i], " --report-file ",
                  paste0(dataloc$fuge.ill.aws[i], ".report")))
    
    system(paste0("centrifuge-kreport -x ~/centrifuge/indices/p+h+v ", dataloc$fuge.ill.aws[i],
                  " >", paste0(dataloc$fuge.ill.aws[i], ".kreport")))
    
    ##nanopore next
    if (!is.na(dataloc$nanopore.fa.aws[i])) {

        system(paste0("centrifuge -p 16 -x ~/centrifuge/indices/p+h+v -U ",dataloc$nanopore.fa.aws[i], 
                      " -S ", dataloc$fuge.nano.aws[i], " --report-file ",
                      paste0(dataloc$fuge.nano.aws[i], ".report")))
        
        system(paste0("centrifuge-kreport -x ~/centrifuge/indices/p+h+v ", dataloc$fuge.nano.aws[i],
                      " >", paste0(dataloc$fuge.nano.aws[i], ".kreport")))
        
    }       
}




