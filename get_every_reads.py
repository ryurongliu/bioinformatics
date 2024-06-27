import glob
import sys
from Bio import SeqIO
import gzip
import pyfastx
import os


"""
Custom subsampling tool for picking one sequence read every N reads.
"""

LINE_UP = '\033[1A'
LINE_CLEAR = '\x1b[2K'


size = 0 
if len(sys.argv) >= 2:
    print(sys.argv[1])
    size = int(sys.argv[1])
    
else:
    print("please input number of reads to skip")
    quit() 
    
    
print("Retrieving 1 read for every", size, "reads")

if not os.path.exists('every_' + str(size) + "_reads"):
   os.makedirs('every_' + str(size) + "_reads")

##get all zipped fastq files 
gzfiles = glob.glob("*fastq.gz")

for file in gzfiles:
    print("processing file", file)
    print()
    
    fastx = pyfastx.Fastq(file)
    #get length of fq
    num_reads = len(fastx)
    
    readfile = file

    writefile = "every_" + str(size) + "_reads/" + "every_" + str(size) + "_" + file
    handle_in = gzip.open(readfile, "rt")
    handle_out = gzip.open(writefile, "wt")
        
    fq = SeqIO.parse(handle_in, "fastq")
    
    i = 0
    for read in fq:
        if i%size == 0:
            print(LINE_UP, end=LINE_CLEAR)
            print('\t', "processing read", i+1, "/", num_reads)
            handle_out.write(read.format("fastq"))
        i+=1 

    handle_in.close()
    handle_out.close()
    
    