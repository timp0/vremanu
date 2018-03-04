#!/bin/bash
##using an ubuntu 16.04 lts r4.4xlarge instance, with full upgrade as of 090517
##also emacs, htop, build-essential, ess, 

if [ "$1" == "conda" ]; then

    sudo apt install awscli
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b
    ##Set forever
    echo 'export PATH=/home/ubuntu/miniconda3/bin:$PATH' >> .bashrc
    ##set for this (don't know why but source of .bashrc doesn't work)
    export PATH=/home/ubuntu/miniconda3/bin:$PATH
   
fi


if [ "$1" == "r.pkgs" ]; then
    
    sudo apt install ess r-base
    sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev default-jdk
    sudo R CMD javareconf
    sudo Rscript ~/vremanu/r_install.R
    
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


if [ "$1" == "kraken.make" ]; then

    mkdir -p /home/ubuntu/src
    cd /home/ubuntu/src/
    git clone https://github.com/DerrickWood/kraken.git
    cd kraken
    mkdir -p /home/ubuntu/kraken/
    ./install_kraken.sh /home/ubuntu/kraken/

    conda install -c bioconda blast
fi

if [ "$1" == "bracken.make" ]; then

    cd ~
    git clone https://github.com/jenniferlu717/Bracken.git
 
    sudo cpan -i List::MoreUtils Parallel::ForkManager   
  

fi



if [ "$1" == "kraken.db.prep" ]; then
    export LD_LIBRARY_PATH="/usr/local/lib"
    cd /home/ubuntu/
    
    rm -R ~/krakendb/
    
    mkdir -p /home/ubuntu/krakendb
    ~/kraken/kraken-build --download-taxonomy --db krakendb
    ~/kraken/kraken-build --download-library bacteria --db krakendb
    ~/kraken/kraken-build --download-library plasmid --db krakendb
    ~/kraken/kraken-build --download-library viral --db krakendb
    ~/kraken/kraken-build --download-library archaea --db krakendb
    ~/kraken/kraken-build --download-library human --db krakendb

    #filter with Dustmasker and convert low complexity regions to N's with Sed (skipping headers)
    for i in `find krakendb \( -name '*.fna' -o -name '*.ffn' \)`
    do
	echo $i
	dustmasker -in $i -infmt fasta -outfmt fasta | sed -e '/>/!s/a\|c\|g\|t/N/g' > tempfile
	mv tempfile $i
    done
    
    ##~/kraken/kraken-build --threads 30 --build --jellyfish-hash-size 25600M --db krakendb
    
fi

if [ "$1" == "kraken.db.k31" ]; then
    export LD_LIBRARY_PATH="/usr/local/lib"
    cd /home/ubuntu/

    ~/kraken/kraken-build --threads 30 --build --jellyfish-hash-size 25600M --db krakendb

    aws s3 sync krakendb/ s3://timpawsanalysis/170913_krakendb.masked.fullstd/ 

fi

if [ "$1" == "kraken.db.k24" ]; then
    export LD_LIBRARY_PATH="/usr/local/lib"
    cd /home/ubuntu/

    ~/kraken/kraken-build --threads 30 --build --kmer-len 24 --jellyfish-hash-size 25600M --db krakendb

    aws s3 sync krakendb/ s3://timpawsanalysis/170913_krakendb.masked.fullk24/
    
fi

if [ "$1" == "kraken.db.tx31" ]; then
    cd ~    
    rm -R ~/krakendb/
    aws s3 sync s3://timpawsanalysis/170913_krakendb.masked.fullstd/ krakendb/
fi

if [ "$1" == "kraken.db.tx24" ]; then
    cd ~    
    rm -R ~/krakendb/    
    aws s3 sync s3://timpawsanalysis/170913_krakendb.masked.fullk24/ krakendb/
fi


if [ "$1" == "bracken.k24" ]; then
    
    find -L krakendb/library  -name "*.fna" -o name "*.fa" -o -name "*.fasta" -exec cat {} +  >allseqs.fa    
    ~/kraken/kraken --db=krakendb --fasta-input --threads=30 \
	   <( find -L krakendb/library -name "*.fna" -exec cat {} + )  > krakendb/database.kraken
    perl Bracken/count-kmer-abundances.pl --threads 30 --db=krakendb/ \
	 --read-length=2000 krakendb/database.kraken > krakendb/database2kmer.kraken_cnts

    perl Bracken/count-kmer-abundances.pl --threads 30 --db=krakendb/ \
	 --read-length=250 krakendb/database.kraken > krakendb/database250mer.kraken_cnts

    python Bracken/generate_kmer_distribution.py -i krakendb/database2kmer.kraken_cnts -o krakendb/KMER_DISTR.nano.TXT
    python Bracken/generate_kmer_distribution.py -i krakendb/database250mer.kraken_cnts -o krakendb/KMER_DISTR.ill.TXT

    aws s3 sync krakendb/ s3://timpawsanalysis/170913_krakendb.masked.fullk24/
    
fi

if [ "$1" == "bracken.k31" ]; then

    ~/kraken/kraken --db=krakendb --fasta-input --threads=30 \
	   <( find -L krakendb/library -name "*.fna" -exec cat {} + )  > krakendb/database.kraken
    perl Bracken/count-kmer-abundances.pl --threads 30 --db=krakendb/ \
	 --read-length=2000 krakendb/database.kraken > krakendb/database2kmer.kraken_cnts

    perl Bracken/count-kmer-abundances.pl --threads 30 --db=krakendb/ \
	 --read-length=250 krakendb/database.kraken > krakendb/database250mer.kraken_cnts

    python Bracken/generate_kmer_distribution.py -i krakendb/database2kmer.kraken_cnts -o krakendb/KMER_DISTR.nano.TXT
    python Bracken/generate_kmer_distribution.py -i krakendb/database250mer.kraken_cnts -o krakendb/KMER_DISTR.ill.TXT

    aws s3 sync krakendb/ s3://timpawsanalysis/170913_krakendb.masked.fullstd/

    

fi




if [ "$1" == "phylo" ]; then

    conda create -n phylo r python=3.6
    source activate phylo
    conda install r-essentials
    Rscript ~/vremanu/r_install.R
    conda install kraken-biom bioconductor-phyloseq
    
    source deactivate  

fi



if [ "$1" == "clean.db" ]; then
    echo "Cleaning db"
    cd /home/ubuntu
    rm -Rf krakendb

fi


if [ "$1" == "clean.all" ]; then
    echo "Cleaning all"
    cd /home/ubuntu
    rm -Rf kraken
    rm -Rf src
    rm -Rf kraken.db

fi


