#!/bin/bash

##From Steph
#blastn -query /home/shao4/CARD_database/AR-genes.fa -db ./VRE_fasta/sample -out ./VRE_fasta/VRE_blast -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 7

##Read database
#makeblastdb -parse_seqids -in 160223_VRE7_recall_2dhq.fa -out VRE7blast -dbtype nucl -input_type fasta

##CARD database
makeblastdb -in ../CARD/nucleotide_fasta_protein_homolog_model.fasta -out CARDblast -dbtype nucl -input_type fasta
lastdb -uNEAR -R01 CARDlast ../CARD/nucleotide_fasta_protein_homolog_model.fasta 

##Blast CARD against reads
blastn -db VRE7blast -query ../CARD/nucleotide_fasta_protein_homolog_model.fasta -out blastry2 -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 6 -num_threads=8
##Blast reads against CARD
cat 160223_VRE7_recall_2dhq.fa | parallel --block 100k --recstart '>' --pipe \
    blastn -db CARDblast -query - -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 \
    -outfmt '"6 qseqid sseqid pident length qlen slen evalue bitscore"' > nanoblast2


#last-train -P10 CARDlast 160223_VRE7_recall_2dhq.fa  > CARD.par
#cp /atium/Data/Nanopore/Analysis/161202_last/r73.par .

lastal -P10 -p r73.par CARDlast 160223_VRE7_recall_2dhq.fa | last-split -m1e-6 > lastcard.maf

maf-convert sam lastcard.maf | \
    samtools view -T ../CARD/nucleotide_fasta_protein_homolog_model.fasta -b - | \
    samtools sort - -o last.sorted.bam

maf-convert blasttab lastcard.maf >lastblast.tsv
maf-convert tab lastcard.maf >lastcard.tsv

bwa mem -t 8 -x ont2d ../CARD/nucleotide_fasta_protein_homolog_model.fasta ../151004_VRE10_recall_called_2dhq.fastq.gz | samtools view -S -b - | samtools sort - -o vre7_card.sorted.bam

~/Code/timp_nanopore/metagenomics/fa2gb.py -i nucleotide_fasta_protein_homolog_model.fasta >nucleotide_fasta_protein_homolog_model.gb

##CARD BLAST (RGI)
python ~/Code/rgi_card/rgi.py -t read -i 151004_VRE10_recall_called_2dhq.fastq.gz -o try -e 0


##Kraken CARD (Florian)

#/usr/local/kraken/kraken --db ~/Code/card-krakendb/ ~/Code/card-krakendb/library/* > test.kraken
#/usr/local/kraken/kraken-report --db ~/Code/card-krakendb/ test.kraken > test.report

/usr/local/kraken/kraken --db ~/Code/card-krakendb/ 160223_VRE7_recall_2dhq.fa >VRE7.card.kraken
/usr/local/kraken/kraken-report --db ~/Code/card-krakendb/ VRE7.card.kraken >VRE7.card.kraken.report
