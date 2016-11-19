#!/bin/bash

#what have I done so far with KPC stuff
# align against all of the KPC genes
bowtie2 -q -x ./CARD_database/KPC -1 /mithril/Data/NGS/Raw/160120_VRE710/VRE7_S1_L001_R1_001.fastq.gz -2 /mithril/Data/NGS/Raw/160120_VRE710/VRE7_S1_L001_R2_001.fastq.gz -S ~/VRE_KPC/VRE7_KPC.sam

samtools view -bS ~/VRE_KPC/VRE7_KPC.sam > ~/VRE_KPC/VRE7_KPC.bam

samtools index ~/VRE_KPC/VRE7_KPC.bam

samtools view -b -F 4 ~/VRE7_KPC.bam > VRE7_mapped.bam

samtools sort VRE7_mapped.bam > VRE7_mapped.sorted.bam

bedtools bamtofastq -i VRE7_mapped.bam -fq VRE7_KPC.fastq 

#convert to fasta

clustalo -i ~/VRE_KPC/VRE7_KPC.fasta --distmat-out=~/VRE_KPC/VRE7KPC.mat --guidetree-out=./VRE7KPC.dnd --force --full 

#get the mito_tree.R file from timp_nanopore, Yunfan's clustal_view.R file

Rscript clustal_view.R

samtools sort VRE7_mapped.bam > VRE7_mapped.sorted.bam

samtools mpileup -uf ~/CARD_database/KPC.fasta VRE7_mapped.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > KPC.cons.fq

samtools mpileup -ugf ~/CARD_database/KPC.fasta VRE7_mapped.sorted.bam | bcftools call -vcOv -o VRE7.vcf

#try alignment of just KPC2, build consensus
bowtie2 -q -x KPC2 -1 /mithril/Data/NGS/Raw/160120_VRE710/VRE7_S1_L001_R1_001.fastq.gz -2 /mithril/Data/NGS/Raw/160120_VRE710/VRE7_S1_L001_R2_001.fastq.gz -S ~/VRE_KPC/VRE7_KPC2.sam

samtools view -bS VRE7_KPC2.sam > VRE7_KPC2.bam

samtools index VRE7_KPC2.bam
samtools sort VRE7_KPC2.bam >VRE7_KPC2.sorted.bam







