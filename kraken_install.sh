#!/bin/bash

##using an ubuntu 16.04 lts r4.4xlarge instance, with full upgrade as of 090517
##also emacs, htop, build-essential

if [ "$1" == "kraken.make" ]; then

    mkdir -p /home/ubuntu/src
    cd /home/ubuntu/src/
    git clone https://github.com/DerrickWood/kraken.git
    cd kraken
    mkdir -p /home/ubuntu/kraken/
    ./install_kraken.sh /home/ubuntu/kraken/
fi
   
if [ "$1" == "jellyfish.make" ]; then

    mkdir -p /home/ubuntu/src/jellyfish
    cd /home/ubuntu/src/jellyfish
    wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
    tar -xzf jellyfish-1.1.11.tar.gz
    cd jellyfish-1.1.11
    ./configure
    make
    sudo make install
fi

if [ "$1" == "conda" ]; then

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
    ##Set forever
    echo 'export PATH=/home/ubuntu/miniconda3/bin:$PATH' >> .bashrc
    ##set for this (don't know why but source of .bashrc doesn't work)
    export PATH=/home/ubuntu/miniconda3/bin:$PATH

    conda install -y r-essentials 
    
fi



if [ "$1" == "kraken.db" ]; then
    export LD_LIBRARY_PATH="/usr/local/lib"
    cd /home/ubuntu/
    mkdir -p /home/ubuntu/krakendb
    ~/kraken/kraken-build --download-taxonomy --db krakendb
    ~/kraken/kraken-build --download-library bacteria --db krakendb
    ~/kraken/kraken-build --download-library plasmids --db krakendb
    ~/kraken/kraken-build --download-library viral --db krakendb
    ~/kraken/kraken-build --download-library archaea --db krakendb
    ~/kraken/kraken-build --download-library human --db krakendb
    ~/kraken/kraken-build --threads 14 --build --jellyfish-hash-size 12800M --db krakendb 
    
fi



if [ "$1" == "clean" ]; then

    cd /home/ubuntu
    rm -Rf kraken
    rm -Rf src
    rm -Rf kraken.db

fi


