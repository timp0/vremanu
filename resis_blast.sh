#!/bin/bash

##From Steph
#blastn -query /home/shao4/CARD_database/AR-genes.fa -db ./VRE_fasta/sample -out ./VRE_fasta/VRE_blast -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 7

makeblastdb -parse_seqids -in 160223_VRE7_recall_2dhq.fa -out VRE7blast -dbtype nucl -input_type fasta

blastn -db VRE7blast -query ../CARD/nucleotide_fasta_protein_homolog_model.fasta -out blastry2 -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 -outfmt 6 -num_threads=8
