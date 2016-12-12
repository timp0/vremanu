library(tidyverse)


##ok - taking vcf output from gingr

vcfs=read_tsv("/atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/P_2016_12_02_130357018963/out.vcf")
##vcfs=read_tsv("/atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/P_2016_12_05_231905287037/SNP_3ref.vcf", comment="##")

##ok - for each sample column, need to look at delta for all other sample columns, then sum to get that square on the table

##Interesting names are names(vcfs)[10:22]

comp=as_data_frame(t(combn(names(vcfs[10:22]),2))) %>%
    mutate(val=colSums(abs(vcfs[comp$V1]-vcfs[comp$V2])))

comp=rbind(comp, data_frame(V1=comp$V2, V2=comp$V1, val=comp$val))

write_csv(spread(comp, V2, val), "/atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/P2016_12_05_231905287035/outsnp_3ref.csv")




##ok - we need opposite comparisons to make a square table - fastest, but probably not best way:
y=z[,c(2,1,3)]

x=rbind(z,y)

x=distinct(x)

spread(x, V2, val)


for (i in 10:21) {
    for (j in (i+1):22) {
        print(paste0(names(vcfs)[i], "-", names(vcfs)[j]))
        
    }
}
