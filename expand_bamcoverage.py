import csv
import os
import glob 
import sys

"""
Expands bamCoverage output to a uniform bin size.


Usage:

1. Make a folder and place the .bedgraph files you want to expand inside of it. 
2. Place this script inside that folder.
3. In terminal, go inside the folder. 
4. Run the script with the following command (no quotes): "python expand_bamcoverage.py"
5. To specify bin size (ex. 100), use the following command (no quotes): "python expand_bamcoverage.py 100"
    - default bin size is 50
    
Output files will be saved as .csv to a subfolder called "expanded", with the same filenames + the suffix "_expanded". 


"""

def expand_file(filename, bin_size):
    newrows = []

    with open(filename, 'r') as f:
        lines = f.readlines()

        #last entry needs checking for bin size 
        for row in lines[:-1]:
            #print(row)
            row = row.split("\t")
            chr_name = row[0]
            bin_start = int(row[1])
            bin_end = int(row[2])
            bin_amt = float(row[3])

            #figure out how many bins to split into
            diff = bin_end - bin_start
            num_bins = int(diff / bin_size)

            #make new rows for each new bin
            for i in range(num_bins):
                new_start = bin_start + i*bin_size
                new_end = bin_start + (i+1)*bin_size
                new_row = [chr_name, new_start, new_end, bin_amt]
                newrows.append(new_row)
                #print("\t", new_start, new_end)

        
        #deal with last row separately
        lastrow = lines[-1]
        #print(lastrow)
        chr_name = row[0]
        bin_start = int(row[1])
        bin_end = int(row[2])
        bin_amt = float(row[3])

        diff = bin_end - bin_start 
        num_bins = int(diff / bin_size)
        remainder = diff - num_bins * bin_size #last bin will be remainder-sized 

        #make new rows for each new full-sized bin
        for i in range(num_bins):
            new_start = bin_start + i*bin_size
            new_end = bin_start + (i+1)*bin_size
            new_row = [chr_name, new_start, new_end, bin_amt]
            newrows.append(new_row)
            #print("\t", new_start, new_end)

        new_lastrow = [chr_name, new_end, bin_end, bin_amt]
        newrows.append(new_lastrow)
        #print("\t", new_end, bin_end)
        
    return newrows



if __name__ == "__main__": 

    #set binsize as command line argument 
    bin_size = 50
    if len(sys.argv) >= 2:
        bin_size = int(sys.argv[1])

    print()
    print("expanding coverage tracks to bin size", bin_size)



    #make output folder
    outfolder = "expanded/"
    if not os.path.exists(outfolder):
        print("making output folder -->", outfolder) 
        os.mkdir(outfolder)

    else:
        print("output folder -->", outfolder)

    print()

    #grab all bedgraph files 
    files = glob.glob("*.bedgraph") 
    print(len(files), "file(s) found for expansion") 

    print()

    #expand all bedgraph files 
    for i in range(len(files)):
        file = files[i]
        print("expanding file", i+1, "of", len(files), "-->", file)
        expanded = expand_file(file, bin_size) 

        outfile = outfolder + file.split(".")[0] + "_expanded.csv"
        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            for row in expanded:
                writer.writerow(row)
        print("\t completed")



