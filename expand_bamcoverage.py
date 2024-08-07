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
    lines = []

    with open(filename, 'r') as f:
        lines = f.readlines()


    for row in lines:
        row = row.split("\t")
        chr_name = row[0]
        bin_start = int(row[1])
        bin_end = int(row[2])
        bin_amt = float(row[3])

        #figure out how many full bins to split into
        diff = bin_end - bin_start
        num_bins = int(diff / bin_size)

        #check for remainder
        remainder = diff - num_bins * bin_size 

        #make new rows for each new full-sized bin
        for i in range(num_bins):
            new_start = bin_start + i*bin_size
            new_end = bin_start + (i+1)*bin_size
            new_row = [chr_name, new_start, new_end, bin_amt]
            newrows.append(new_row)
            #print("\t", new_start, new_end)

        #if there is a remainder, make a remainder bin
        if remainder:
            #print(row)
            new_rem_bin = [chr_name, new_end, bin_end, bin_amt]
            newrows.append(new_rem_bin)
            #print(new_rem_bin)
        
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

        outfile = outfolder + file.split(".")[0] + "_expanded.txt"
        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            for row in expanded:
                writer.writerow(row)
        print("\t completed")



