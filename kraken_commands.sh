#kraken database commands

#make AWS instance (use the AMI Timp has created)

#mount database from AWS S3

sudo mkfs -t ext4 /dev/xvdb
sudo mount /dev/xvbd data
sudo chmod 777 data
mkdir data/std
sync s3://timpawsanalysis/160309_krakenstd/ data/std

kraken-build --download-library human --db data/std
kraken-build --rebuild --db data/std #okay, this step will take quite a bit of time sadness
