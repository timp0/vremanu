library(tidyverse)


##ok - taking vcf output from gingr

vcfs=read_tsv("~/Dropbox/Temp/parsnp.vcf", comment="##")

##ok - for each sample column, need to look at delta for all other sample columns, then sum to get that square on the table

##Interesting names are names(vcfs)[10:22]

comp=as_data_frame(t(combn(names(vcfs[10:22]),2))) %>%
    mutate(val=colSums(abs(vcfs[V1]-vcfs[V2])))

comp=rbind(comp, data_frame(V1=comp$V2, V2=comp$V1, val=comp$val))

write_csv(spread(comp, V2, val), "~/Dropbox/Temp/outsnp.csv")




