#! /usr/bin/env python3


def lq_scat(prefix, data): 
    ##From MATPLOTLIB EXAMPLES
    import matplotlib
    matplotlib.use('PDF')

    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    
    nullfmt = NullFormatter()         # no labels
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    
    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))
    
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # the scatter plot:
    ##Unlikely initial limits
    llimit=[1e6,0.0]
    qlimit=[60.0,0.0]

    for dset in data:
        if len(dset.rlength)>0:
            axScatter.scatter(dset.rlength, dset.qual, c=dset.col, alpha=.3)
            llimit[0]=min(llimit[0], dset.minlen)
            llimit[1]=max(llimit[1], dset.maxlen)
            qlimit[0]=min(qlimit[0], dset.minqual)
            qlimit[1]=max(qlimit[1], dset.maxqual)
    
    ##now flatten nice limits
    llimit=[float(math.floor(llimit[0])), float(math.ceil(llimit[1]))]
    qlimit=[float(math.floor(qlimit[0])), float(math.ceil(qlimit[1]))]

    axScatter.set_xlim(llimit)
    axScatter.set_xlabel('Length (bp)')
    axScatter.set_ylim(qlimit)
    axScatter.set_ylabel('PHRED quality')
    
    ldiv= (llimit[1]-llimit[0])/100
    qdiv= (qlimit[1]-qlimit[0])/100

    xbins = np.arange(llimit[0]-ldiv, llimit[1]+ldiv, ldiv)
    ybins = np.arange(qlimit[0]-qdiv, qlimit[1]+qdiv, qdiv)
    
    for dset in data:
        if len(dset.rlength)>0:
            axHistx.hist(dset.rlength, bins=xbins, histtype='step', color=dset.col)
            axHisty.hist(dset.qual, bins=ybins, orientation='horizontal', histtype='step', color=dset.col)

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    
    plt.savefig((prefix+'_qc.pdf'), format='pdf')

class h5read:
    
    def __init__(self, rtype):
        self.fastq=[]
        self.rlength=-1
        self.avgqual=-1
        self.stime=-1
        self.etime=-1
        self.duration=-1
        self.rtype=rtype
        self.rname=[]
        #self.script=[]
        
        if rtype == "2D":
            self.fqkey='/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
            self.tevkey='/Analyses/Basecall_1D_000/BaseCalled_template/Events'
            self.cevkey='/Analyses/Basecall_1D_000/BaseCalled_complement/Events'
        elif rtype == "1D":
            self.fqkey='/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
            self.tevkey='/Analyses/Basecall_1D_000/BaseCalled_template/Events'
            ##No complement for 1D
            self.cevkey='/Analyses/Basecall_1D_000/BaseCalled_template/Events'
        elif rtype == "BC":
            self.bcsummary = '/Analyses/Barcoding_000/Summary/barcoding'
            self.fqkey='/Analyses/Barcoding_000/Barcoding/Fastq'
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


        if self.rtype=="BC":
            self.bcid=hdf[self.bcsummary].attrs.get('barcode_arrangement').decode("utf-8")
            if self.bcid=="unclassified":
                self.bcnum=0
            else: 
                self.bcnum=int(self.bcid[2:4])
        else:
            self.bcnum=0


class runsum:
        
    def __init__(self, col, prefix, suffix, fqwrite=True, twrite=False):
        self.nreads=0
        self.qual=[]
        self.rname=[]
        self.rlength=[]
        self.stime=[]
        self.duration=[]
        self.etime=[]
        self.col=col
        self.fqwrite=fqwrite
        self.twrite=twrite
        if fqwrite:
            self.fout=gzip.open(prefix+'_'+suffix+'.fastq.gz', 'wb')
        if twrite:
            self.tout=gzip.open(prefix+'_'+suffix+'.tout.csv.gz', 'wt')
            delim=','
            self.tout.write(delim.join(['readname', 'length', 'avgqual', 'start.time', 'end.time', 'duration'])+'\n') 
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

        if self.fqwrite:             
            self.fout.write(read.fastq)
            if read.fastq.decode('ascii')[-1]!="\n":
                self.fout.write(bytes("\n", 'ascii'))
        if self.twrite:            
            delim=','
            self.tout.write(delim.join([read.rname,str(read.rlength),str(read.avgqual), str(read.stime), str(read.etime), str(read.duration)])+'\n')

    def finalstat(self):
        import statistics
        import numpy as np

        ##Empty check
        if self.rlength:
            self.totyield=sum(self.rlength)/1e6
            
            self.avglen=statistics.mean(self.rlength)
            self.minlen=np.percentile(self.rlength, 5)
            self.maxlen=np.percentile(self.rlength, 95)
            

            self.avgqual=statistics.mean(self.qual)
            self.minqual=np.percentile(self.qual, 5)
            self.maxqual=np.percentile(self.qual, 95)

        

    def close(self):
        self.fout.close()


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


parser = argparse.ArgumentParser( description='Extract fastq from fast5 and make some qc plots of the length/quality')
parser.add_argument('--input', '-i', type=str, required=True, help='input location of tarball containing fast5 files')
parser.add_argument('--barcode', '-b',  action='store_true', help='Sample is barcoded')
parser.add_argument('--oned', '-1', action='store_true', help='Put out 1D') 
parser.add_argument('--time', '-t', action='store_true', help='Put out 1D') 


args=parser.parse_args()


qscorethresh=9

##Random number dir location for tar
tmpdir=str(random.randint(1e6, 9e6))+'/'

os.mkdir(tmpdir)
prefix=os.path.basename(args.input)
prefix=prefix.replace('.tar.gz','')
print(prefix)
    
##open tarball to random dir
tarball=tarfile.open(name=args.input, mode='r')
    
for member in tarball.getmembers():
    if member.isreg():
        member.name=os.path.basename(member.name)
        tarball.extract(member,tmpdir)

if (args.barcode):
    bccolors=["black", "blue", "orange", "green", "red", "purple", "grey", "cyan", "fuchsia", "brown", "goldenrod", "darkgreen", "orangered"]
    allbc=[runsum(bccolors[i], prefix, "BC"+str(i).zfill(2)) for i in range(13)]
else:
    all2d=runsum("blue", prefix, "2d", twrite=args.time)
    pass2d=runsum("green", prefix, "2dhq", twrite=args.time)
    fail2d=runsum("red", prefix, "2dlq", twrite=args.time)

if (args.oned):
    alltemp=runsum("blue", prefix, "", fqwrite=False)
    has2dtemp=runsum("green", prefix, "", fqwrite=False)
    no2dtemp=runsum("red", prefix, "1d", twrite=args.time)

filelist=glob.glob(tmpdir+'*fast5')

nreads=0

chemver=""

        
for filename in filelist:
            
    try: 
        ##Load file
        hdf = h5py.File(filename, 'r')
        
        nreads += 1

        if (args.barcode):
            read2d=h5read('BC')
            read2d.check(hdf)
        else:
            read2d=h5read('2D')
            read2d.check(hdf)
            
        if (args.oned):
            read1d=h5read('1D')
            read1d.check(hdf)

            if (read1d.fqcheck):
                read1d.extract(hdf)
                #Extract 
                alltemp.addread(read1d)
                #addread
                if (read2d.fqcheck):
                    has2dtemp.addread(read1d)
                else:
                    no2dtemp.addread(read1d)
            
        if (read2d.fqcheck):
            read2d.extract(hdf)

            if (args.barcode):
                allbc[read2d.bcnum].addread(read2d)
            else:
                all2d.addread(read2d)
                if read2d.avgqual>qscorethresh:
                    pass2d.addread(read2d)
                else:
                    fail2d.addread(read2d)
        
        
                    
        hdf.close()        
    except Exception as e:
        print(e)
        print(filename+" failed to open!!")
        
                          

    


report=open((prefix+'.qc.csv'), 'w')
report.write('Report,Number of Reads,Total Yield(MB),Average Length,Average Quality\n')
report.write('Total Reads,'+str(nreads)+ ',,,\n')


if (args.oned):
    no2dtemp.close()
    fail1d=nreads-alltemp.nreads
    ##Plot 1D
    no2dtemp.finalstat()
    lq_scat(prefix+".1D", [alltemp, has2dtemp, no2dtemp])    
    report.write('1D Basecalling Failed,'+str(fail1d)+',,,\n')
    ##Calc stats
    report.write('1D only,'+str(no2dtemp.nreads)+','+ str(no2dtemp.totyield) + ',' + str(no2dtemp.avglen) + ',' + str(no2dtemp.avgqual) + '\n')

if (args.barcode):
    nfail=nreads
    for i in range(len(allbc)):
        ##Close out fastq
        allbc[i].close()
        ##Calc staats
        allbc[i].finalstat()
        nfail=nfail-allbc[i].nreads
        report.write('2D Basecalling BC'+str(i).zfill(2)+' Total,'+str(allbc[i].nreads)+','+ str(allbc[i].totyield) + ',' + str(allbc[i].avglen) + 
                     ',' + str(allbc[i].avgqual) + '\n')
    
    lq_scat(prefix+".BC", allbc)

else:
    all2d.close()
    pass2d.close()
    fail2d.close()

    ##Calc stats
    all2d.finalstat()
    pass2d.finalstat()
    fail2d.finalstat()


    ##Plot 2D 
    lq_scat(prefix+".2D", [all2d, pass2d, fail2d])
    report.write('2D Basecalling Failed,'+str(nreads-all2d.nreads)+',,,\n')

    ##output stats

    report.write('2D Basecalling Passed Total,'+str(all2d.nreads)+','+ str(all2d.totyield) + ',' + str(all2d.avglen) + ',' + str(all2d.avgqual) + '\n')
    report.write('2D Basecalling Passed HQ,'+str(pass2d.nreads)+','+ str(pass2d.totyield) + ',' + str(pass2d.avglen) + ',' + str(pass2d.avgqual) + '\n')
    report.write('2D Basecalling Passed LQ,'+str(fail2d.nreads)+','+ str(fail2d.totyield) + ',' + str(fail2d.avglen) + ',' + str(fail2d.avgqual) + '\n')




report.close()
    
    
##Clear tmpdir
shutil.rmtree(tmpdir)



