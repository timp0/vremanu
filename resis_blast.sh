#!/bin/bash

##From Steph
#blastn -query /home/shao4/CARD_database/AR-genes.fa -db ./VRE_fasta/sample -out ./VRE_fasta/VRE_blast -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 7

makeblastdb -parse_seqids -in 160223_VRE7_recall_2dhq.fa -out VRE7blast -dbtype nucl -input_type fasta
makeblastdb -in ../CARD/nucleotide_fasta_protein_homolog_model.fasta -out CARDblast -dbtype nucl -input_type fasta
lastdb -uNEAR -R01 CARDlast ../CARD/nucleotide_fasta_protein_homolog_model.fasta 


blastn -db VRE7blast -query ../CARD/nucleotide_fasta_protein_homolog_model.fasta -out blastry2 -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 6 -num_threads=8
blastn -db CARDblast -query 160223_VRE7_recall_2dhq.fa -out nanoblast1 -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 6 -num_threads 8

#last-train -P10 CARDlast 160223_VRE7_recall_2dhq.fa  > CARD.par
#cp /atium/Data/Nanopore/Analysis/161202_last/r73.par .

lastal -P10 -p r73.par CARDlast 160223_VRE7_recall_2dhq.fa | last-split -m1e-6 > lastcard.maf

#maf-convert sam lastcard.maf | \
#    samtools view -T card/mithril/Data/NGS/Reference/lambda/lambda.fasta -b - | \
#    samtools sort - ${outdir}/${prefix}.sorted


cat 160223_VRE7_recall_2dhq.fa | parallel --block 100k --recstart '>' --pipe \
    blastn -db CARDblast -query - -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 \
    -outfmt '"6 qseqid sseqid pident length qlen slen evalue bitscore"' > nanoblast2

cat 1gb.fasta | parallel --block 100k --recstart '>' --pipe blastp -evalue 0.01 -outfmt 6 -db db.fa -query - > results

bwa mem -t 8 -x ont2d CARD/nucleotide_fasta_protein_homolog_model.fasta 151004_VRE10_recall_called_2dhq.fastq.gz | samtools view -S -b - | samtools sort - -o vre7_card.sorted.bam

~/Code/timp_nanopore/metagenomics/fa2gb.py -i nucleotide_fasta_protein_homolog_model.fasta >nucleotide_fasta_protein_homolog_model.gb


python ~/Code/rgi_card/rgi.py -t read -i 151004_VRE10_recall_called_2dhq.fastq.gz -o try -e 0
