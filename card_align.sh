#!/bin/bash

##CARD Alignment

if [ "$1" == "mappers" ]; then

   conda create --name mappers
   source activate mappers
   conda install -c bioconda bowtie2 minimap2 samtools blast
   source deactivate

fi

if [ "$1" == "get.card" ]; then

    cd ~
    mkdir -p carddb
    cd carddb
    wget https://card.mcmaster.ca/download/0/broadstreet-v1.2.1.tar.bz2
    tar -xvjf broadstreet-v1.2.1.tar.bz2
    rm *bz2

    source activate mappers
    minimap2 -t 30 -x map-ont -d card_homolog.mmi nucleotide_fasta_protein_homolog_model.fasta
    bowtie2-build --threads 30 nucleotide_fasta_protein_homolog_model.fasta card_homolog
    makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -out CARDblast -dbtype nucl -input_type fasta
    source deactivate

fi


