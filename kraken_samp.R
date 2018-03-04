library(tidyverse)
library(googlesheets)

rawdir="~/Data/vre"
workdir="~/Data/kraken"

dir.create(workdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_0905")

dataloc$illumina.r1.aws=file.path(rawdir, basename(dataloc$illumina.r1))
dataloc$illumina.r2.aws=file.path(rawdir, basename(dataloc$illumina.r2))

dataloc$nanopore.fa.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(rawdir, basename(dataloc$nanopore.fasta)))


##Setup outfilenames
dataloc$k31.nano.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(workdir, paste0(dataloc$sample, ".nano.k31.krak")))
dataloc$k31.ill.aws=file.path(workdir, paste0(dataloc$sample, ".ill.k31.krak"))

setwd("/home/ubuntu/")

if (TRUE) {
    ##kraken k-mer 31
    unlink("~/krakendb", recursive=TRUE)
    system(paste0("aws s3 sync --quiet --delete s3://timpawsanalysis/170913_krakendb.masked.fullstd/ ~/krakendb/"))

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        system(paste0("~/kraken/kraken --threads 30 --db krakendb/ --paired ", dataloc$illumina.r1.aws[i],
                      " ", dataloc$illumina.r2.aws[i], " >", dataloc$k31.ill.aws[i]))
        system(paste0("~/kraken/kraken-translate --db krakendb/ ", dataloc$k31.ill.aws[i], " >",
                      dataloc$k31.ill.aws[i], ".labels"))
        system(paste0("~/kraken/kraken-report --db krakendb/ ", dataloc$k31.ill.aws[i], " >",
                      dataloc$k31.ill.aws[i], ".report"))
        
        system(paste0("python ~/Bracken/src/est_abundance.py -i ", dataloc$k31.ill.aws[i], ".report -k ~/krakendb/KMER_DISTR.ill.TXT -o ",
                      dataloc$k31.ill.aws[i], ".bracken.txt"))



        
        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("~/kraken/kraken --threads 30 --db krakendb/ ", dataloc$nanopore.fa.aws[i], " >", dataloc$k31.nano.aws[i]))
            system(paste0("~/kraken/kraken-translate --db krakendb/ ", dataloc$k31.nano.aws[i], " >",
                          dataloc$k31.nano.aws[i], ".labels"))
            system(paste0("~/kraken/kraken-report --db krakendb/ ", dataloc$k31.nano.aws[i], " >",
                          dataloc$k31.nano.aws[i], ".report"))

            system(paste0("python ~/Bracken/src/est_abundance.py -i ", dataloc$k31.nano.aws[i], ".report -k ",
                          "~/krakendb/KMER_DISTR.nano.TXT -o ",
                          dataloc$k31.nano.aws[i], ".bracken.txt"))

        }       
    }

}

dataloc$k24.nano.aws=ifelse(is.na(dataloc$nanopore.fasta), NA, file.path(workdir, paste0(dataloc$sample, ".nano.k24.krak")))
dataloc$k24.ill.aws=file.path(workdir, paste0(dataloc$sample, ".ill.k24.krak"))



if (TRUE) {
    ##kraken k-mer 24
    unlink("~/krakendb", recursive=TRUE)
    system(paste0("aws s3 sync --quiet --delete s3://timpawsanalysis/170913_krakendb.masked.fullk24/ ~/krakendb/"))

    for (i in 1:dim(dataloc)[1]) {                    
        ##illumina first
        system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ --paired ", dataloc$illumina.r1.aws[i],
                      " ", dataloc$illumina.r2.aws[i], " >", dataloc$k24.ill.aws[i]))
        system(paste0("~/kraken/kraken-translate --db krakendb/ ", dataloc$k24.ill.aws[i], " >",
                      dataloc$k24.ill.aws[i], ".labels"))
        system(paste0("~/kraken/kraken-report --db krakendb/ ", dataloc$k31.ill.aws[i], " >",
                      dataloc$k24.ill.aws[i], ".report"))

        system(paste0("python ~/Bracken/src/est_abundance.py -i ", dataloc$k24.ill.aws[i], ".report -k ~/krakendb/KMER_DISTR.ill.TXT -o ",
                      dataloc$k24.ill.aws[i], ".bracken.txt"))

        
        
        
        ##nanopore next
        if (!is.na(dataloc$nanopore.fa.aws[i])) {
            system(paste0("~/kraken/kraken --threads 30 --db ~/krakendb/ ", dataloc$nanopore.fa.aws[i], " >", dataloc$k24.nano.aws[i]))
            system(paste0("~/kraken/kraken-translate --db krakendb/ ", dataloc$k24.nano.aws[i], " >",
                          dataloc$k24.nano.aws[i], ".labels"))
            system(paste0("~/kraken/kraken-report --db krakendb/ ", dataloc$k24.nano.aws[i], " >",
                      dataloc$k24.nano.aws[i], ".report"))

            system(paste0("python ~/Bracken/src/est_abundance.py -i ", dataloc$k24.nano.aws[i], ".report -k ",
                          "~/krakendb/KMER_DISTR.nano.TXT -o ",
                          dataloc$k24.nano.aws[i], ".bracken.txt"))

        }       
    }
       
   
}

#fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
#gs_ws_new(fullsheet, ws_title="VRE_171011", input=dataloc)

setwd(workdir)


system("aws s3 sync ~/Data/kraken/ s3://timpawsanalysis/180303_vre/kraken/")
