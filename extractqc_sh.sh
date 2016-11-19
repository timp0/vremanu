#!/bin/bash

rawdir=/mithril/Data/Nanopore/oxford

for namey in 151004_VRE10_recall 160223_VRE7_recall
	     
do

    filey=`ls ${rawdir}/${namey}/*called.tar.gz`
    #samp
    samp=${filey##*/}
    #echo $filey
    
    samp="${samp%%\.tar.gz}"
    echo $samp

    python3.5 /home/shao4/timp_nanopore/timp_nanopore/oxford/oxford_qc.py -t -1 -i ${filey}

done
