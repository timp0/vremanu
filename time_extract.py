#! /usr/bin/env python3

class h5read:
    
    def __init__(self):
        self.fastq=[]
        self.rlength=-1
        self.avgqual=-1
        self.stime=-1
        self.etime=-1
        self.duration=-1
        self.rname=[]
        #self.script=[]
        
        self.fqkey='/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
        self.tevkey='/Analyses/Basecall_1D_000/BaseCalled_template/Events'
        self.cevkey='/Analyses/Basecall_1D_000/BaseCalled_complement/Events'

        
    def check(self, hdf):
        ##Check to see if things are in hdf file
        if self.fqkey in hdf:
            self.fqcheck=True
        else:
            self.fqcheck=False

    def extract(self, hdf):
        import statistics

        #self.script=hdf['/UniqueGlobalKey/tracking_id'].attrs.get('exp_script_name')
        
        self.fastq = hdf[self.fqkey][()]
        ##Python3 wants explicit for ascii, assumed unicode
        self.rname=self.fastq.decode('ascii').split('\n')[0]
        
        qual=self.fastq.decode('ascii').split('\n')[3]
        self.rlength=len(qual)

        qualval = [ord(x) - 33 for x in qual.rstrip('\n')]
        self.avgqual = statistics.mean(qualval)

        #This was too look for the raw events
        #subev=[self.evkey+"/"+x for x in hdf[self.evkey]]
        ##instead get times this way
        self.stime=hdf[self.tevkey].attrs.get('start_time')

        self.etime=hdf[self.cevkey].attrs.get('start_time')+hdf[self.cevkey].attrs.get('duration')
                    
        self.duration=self.etime-self.stime

        

class runsum:
        
    def __init__(self, prefix):
        self.nreads=0
        self.qual=[]
        self.rname=[]
        self.rlength=[]
        self.stime=[]
        self.duration=[]
        self.etime=[]

        self.tout=gzip.open(prefix, 'wt')
        self.tout.write(','.join(['readname', 'length', 'avgqual', 'start.time', 'end.time', 'duration'])+'\n') 

        #stats        
        self.totyield=0
        self.avglen=0
        self.avgqual=0
        self.minlen=0
        self.maxlen=0
        self.minqual=0
        self.maxqual=0

    def addread(self, read):
        self.nreads += 1
        self.rlength.append(read.rlength)
        self.qual.append(read.avgqual)
        self.stime.append(read.stime)
        self.duration.append(read.duration)
        self.etime.append(read.etime)
        self.rname.append(read.rname)
        
        delim=','
        self.tout.write(delim.join([read.rname,str(read.rlength),str(read.avgqual), str(read.stime), str(read.etime), str(read.duration)])+'\n')

    def close(self):
        self.tout.close()


##Main 

##HDF format
import h5py
##Os IO stuff
import os
import shutil
##Tarball 
import tarfile
##File recognize
import glob
##shell util
import shutil
##Random number
import random
##Gzip lib
import gzip
##Stats
import statistics
##Get args
import argparse
##RE
import re
##Tempfile directory generate
import tempfile


parser = argparse.ArgumentParser( description='Extract times/meta of r7.3 fast5')
parser.add_argument('--input', '-i', type=str, required=True, help='input location of tarball containing fast5 files')
parser.add_argument('--outtime', '-o', type=str, required=True, help='output meta info file')

args=parser.parse_args()


##Random number dir location for tar
tmpdir=str(random.randint(1e6, 9e6))+'/'

os.mkdir(tmpdir)
prefix=os.path.basename(args.input)
prefix=prefix.replace('.tar.gz','')
print(prefix)
    
##open tarball to random dir
tarball=tarfile.open(name=args.input, mode='r')

files=0

for member in tarball.getmembers():
    if member.isreg():
        if ("pass" in member.name) and (".fast5" in member.name):
            member.name=os.path.basename(member.name)
            tarball.extract(member,tmpdir)
            files += 1
            if files > 50:
                break



all2d=runsum(args.outtime)
filelist=glob.glob(tmpdir+'*fast5')

nreads=0

for filename in filelist:
            
    try: 
        ##Load file
        hdf = h5py.File(filename, 'r')
        
        nreads += 1

        read2d=h5read()
        read2d.check(hdf)
                    
        if (read2d.fqcheck):
            read2d.extract(hdf)

            all2d.addread(read2d)
           
        hdf.close()


        
    except Exception as e:
        print(e)
        print(filename+" failed to open!!")
        
    

all2d.close()


##Clear tmpdir
shutil.rmtree(tmpdir)



