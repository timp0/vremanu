library(tidyverse)
library(googlesheets)
library(Biostrings)

workdir="/mithril/Data/Nanopore/Analysis/170905_VRE"
outdir="~/Dropbox/Data/Nanopore/170905_VRE"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws=1)

setwd(workdir)


dataloc=dataloc %>%
    mutate(nanopore.fasta=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "nanopore.raw", paste0(sample, ".nanopore.2D.fasta"))))



if (TRUE) {   
    dir.create(file.path(workdir, "nanopore.raw"))

    for (i in 1:dim(dataloc)[1]) {
    
        if (!is.na(dataloc$nanopore.raw[i])) {
            
            ##Run Poretools            
            ##Are there two flowcells?
            two.flow=(!is.na(dataloc$nanopore.raw2[i]))
            system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27; poretools fasta --type 2D ", dataloc$nanopore.raw[i], " >", dataloc$nanopore.fasta[i]))))
            
            if (two.flow) {
                system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27; poretools fasta --type 2D ", dataloc$nanopore.raw2[i], " >>", dataloc$nanopore.fasta[i]))))
            }
        }
    }
}

##Plot/aquire stats on nanopore fasta/q
if (TRUE) {

    dataloc=dataloc %>%
        mutate(nanopore.n50=NA) %>%
        mutate(nanopore.yield=NA)
    
    
    for (i in 1:dim(dataloc)[1]) {
        if (dataloc$nanopore.fasta[i] != "NA") {
            fa=readDNAStringSet(dataloc$nanopore.fasta[i])
            dataloc$nanopore.n50[i]=N50(width(fa))
            dataloc$nanopore.yield[i]=sum(as.numeric(width(fa)))/1e6
        }     
    }
}

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
gs_ws_new(fullsheet, ws_title="VRE_0905", input=dataloc)
