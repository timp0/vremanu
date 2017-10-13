#!/bin/bash


if [ "$1" == "kraken" ]; then
    mkdir -p ~/Data/kraken
    aws s3 sync s3://timpawsanalysis/171012_vre/kraken/ ~/Data/kraken/

fi

