import pyfastx
import glob
import csv

"""

Counts the number of telomere repeats (TTAGGG or CCCTAA) in each read, and counts how many reads have telomere repeats as a % of the total number of reads in each sample. Outputs counts to a .csv file. 

"""


forward_tseq = "TTAGGG"
backward_tseq = "CCCTAA"

LINE_UP = '\033[1A'
LINE_CLEAR = '\x1b[2K'



#get all zipped fastq files 
gzfiles = glob.glob("*fastq.gz")


for k in range(len(gzfiles)):
    
    print("processing file", k+1, "/", len(gzfiles))
    file = gzfiles[k]
    print('\t', file)
    
    
    #get fastq data from zip 
    fq = pyfastx.Fastq(file)
    
    #get length of fq
    num_reads = len(fq)
    print('\t', "number of reads:", num_reads)
    
    #setup counts for this file 
    forward_repeats = {}
    for i in range(26):
        forward_repeats[i] = 0

    backward_repeats = {}
    for i in range(26):
        backward_repeats[i] = 0
        
    print()
    #go through all reads and look for repeats 
    for j in range(num_reads):
        print(LINE_UP, end=LINE_CLEAR)
        print('\t', "processing read", j+1, "/", num_reads)
        
        
        read = fq[j]
        seq = read.seq
        
        #look for forward repeats 
        if forward_tseq in seq:
            i = 1
            while forward_tseq*i in seq:
                i+=1
            num_repeats = i - 1
            forward_repeats[num_repeats] += 1
        
        else:
            forward_repeats[0] += 1
            
            
        #look for backward repeats
        if backward_tseq in seq:
            i = 1
            while backward_tseq*i in seq:
                i += 1
            
            num_repeats = i - 1
            backward_repeats[num_repeats] += 1
        
        else:
            backward_repeats[0] += 1
        
        
            
            
            
    #output to csv file 
    
    print('\t', "writing to csv...")
    csv_filename = file.replace(".gz", "").replace(".fastq", "") + ".csv"
    with open(csv_filename, "w") as csvfile:
        outwriter = csv.writer(csvfile)
        outwriter.writerow(["Total Reads", num_reads])
        outwriter.writerow([])
        outwriter.writerow(["Forward Repeat Counts (TTAGGG)"])
        outwriter.writerow(["number of repeats", "reads counted", "percentage of total reads"])
        for num, counts in forward_repeats.items():
            outwriter.writerow([num, counts, str(counts / num_reads * 100) + "%"])
        outwriter.writerow([])
        outwriter.writerow(["Backward Repeat Counts (CCCTAA)"])
        outwriter.writerow(["number of repeats", "reads counted", "percentage of total reads"])
        for num, counts in backward_repeats.items():
            outwriter.writerow([num, counts, str(counts / num_reads * 100) + "%"])
    
    print('\t', "saved to", csv_filename)
