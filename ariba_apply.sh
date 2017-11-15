#!/bin/bash



if [ "$1" == "get.card" ]; then
    cd ~
    source activate ariba
    mkdir -p ariba-work
    cd ariba-work
    ariba getref card out.card
    ariba prepareref -f out.card.fa -m out.card.tsv outcard.prepareref
    source deactivate
    

fi

