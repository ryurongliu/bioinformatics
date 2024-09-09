from matplotlib import pyplot as plt
import pandas as pd
import numpy as np 
import os

import glob
import json

import argparse


#-------------------------------------------------------------------------------------------


#some constants

info_file = "ref/info.json"
subsample_nums = [1, 2, 3, 5]
control_names = ['SM-HA', 'SM-IgG']
types = ['avg', 'fold']

#-------------------------------------------------------------------------------------------

#for finding and saving initial constants of overall dataset

def name_subfolder(num):
    
    #returns subsample folder path by subsample number
    #(edit naming convention for future use)
    
    return '../May24_TRR_' + str(num) + '_mil/' #end folderpath with forward slash


def name_file(num, cntrl, type):
    
    #returns file name by control name, subsample number, type
    #(edit naming convention for future use)
    
    if type == 'avg':
        type_ps = 'averages'
    elif type == 'fold':
        type_ps = 'fold_enrichment'
    
    return cntrl + '_control_' + str(num) + "_mil_" + type_ps + ".txt"


def get_files(ss_nums):
    
    #returns files as dictionary
    #files[subnum][control group][type] 
    
    sub_folders = glob.glob("*mil") #subsample folders
    
    files = {
    }

    for i in ss_nums:
        
        sub_folder = name_subfolder(i)
        
        files[i] = {
        }
        
        for j in control_names:
            files[i][j] = {   
            }
            
            for k in types:
                files[i][j][k] = sub_folder + name_file(i, j, k)
                #print(os.path.exists(files[i][j][k]))
                
                
    return files

def get_sample_names(files): 
    
    #return samplenames listed by control 
    #sample_names[control]
    
    sample_names = {
    }
    
    for j in control_names:
        
        file = files[subsample_nums[0]][j][types[0]] #pick first file to get sample names from
        data = pd.read_csv(file, sep="\t") #read file
        sample_names[j] = [x[:-(len(types[0])+1)] for x in 
                           list(data) if types[0] in x] #parse names from headers
        
        
    return sample_names

def get_chr_info(files):
    
    #return chromosome names 
    file = files[subsample_nums[0]][control_names[0]][types[0]]
    data = pd.read_csv(file, sep="\t")
    
    chr_locs = dict.fromkeys(data["chromosome"])
    chr_names = list(chr_locs.keys())
    
    #iterate through to get index positions of each chromosome
    curr_chr = chr_names[0]
    chr_locs[curr_chr] = [0]
    for i in range(len(data['chromosome'])):
        this_chr = data['chromosome'][i]
        if this_chr != curr_chr:
            chr_locs[curr_chr].append(i-1)
            chr_locs[this_chr] = [i]
            curr_chr = this_chr

    chr_locs[curr_chr].append(i)
    
    return chr_names, chr_locs





def find_info():
    
    #find files, sample names, and chromosomes and save to info.json for easy access
    files = get_files(subsample_nums)
    sample_names = get_sample_names(files)
    chr_names, chr_locs = get_chr_info(files)
    
    info = {
        "files": files,
        "sample_names": sample_names,
        "chr_names": chr_names,
        "chr_locs": chr_locs
    }
    
    with open(info_file, "w") as f:
        json.dump(info, f)
    
    return 

def get_info(file):
    
    #return files, sample names, and chromosomes 
    with open(file, 'r') as f:
        info = json.load(f)
        
    files = info['files']
    sample_names = info['sample_names']
    chr_names = info['chr_names']
    chr_locs = info['chr_locs']
    
    return files, sample_names, chr_names, chr_locs

def return_info():
    
    #look for info and return files, sample names, chromosomes 
    
    #if info doesn't exist, find it
    if not os.path.exists(info_file):
        print('finding info about this dataset...')
        find_info()
        print('info saved to', info_file)
    
    print('using', info_file)
    files, sample_names, chr_names, chr_locs = get_info(info_file)
    
    return files, sample_names, chr_names, chr_locs

#-------------------------------------------------------------------------------------------

#for parsing config file

def parse_config(file, sample_names, chr_names, verbose=True):
    
    print("reading", file)
    
    #open config file and check fields for proper inputs 
    with open(file, "r") as f:
        config = json.load(f)
    
    data_config = config['data']
    plot_config = config['plot']
    
    parse_data_config(data_config, sample_names, chr_names, verbose=verbose)
    parse_plot_config(plot_config, data_config, verbose=verbose)
    
    return config

def parse_data_config(d_config, sample_names, chr_names,verbose=True):
    
    #parse data parameters of config 
    
    #type
    assert(d_config['type'] in types), \
    "Do not recognize type " + str(d_config['type']) + " in config, must be avg or fold"
    
    #control
    assert(d_config['control'] in control_names), \
    "do not recognize control parameter " + str(d_config['control']) + " in config"
    
    #sample
    assert(d_config['sample'] in sample_names[d_config['control']]), \
    "sample " + str(d_config['sample']) + " listed in config is unknown"
    
    #chr
    for chr_name in d_config['chromosomes']:
        assert(chr_name in chr_names or chr_name == "all"), \
        "chromosome " + str(chr_name) + " listed in config is unknown"
    
    #subsamples
    for s in d_config['subsamples']:
        assert(s in subsample_nums), "subsampling number " + str(s) + " in config not found"
        
    
    if verbose:
        #print data settings to console    
        print()
        print("-----------Data Settings------------")
        print()
        print("        type:", d_config['type'])
        print("     control:", d_config['control'])
        print("      sample:", d_config['sample'])
        print()
        print(" chromosomes:")
        print("\n".join("\t    " + str(x) for x in d_config['chromosomes']))
        print()
        print("  subsamples:")
        print("\n".join("\t    " + str(x) for x in d_config['subsamples']))
    
    return

def parse_plot_config(p_config, d_config, verbose=True):
    
    #parse plot parameters of config 
    
    #filename
    assert(type(p_config['filename']) is str), \
    "filename in config must be given as a string"
    
    assert(type(p_config['figure_size']) is tuple and \
           len(p_config['figure_size']) == 2 and \
           type(p_config['figure_size'][0]) is int or float), \
    "figure_size in config must be given as tuple (W, H) with int or float"
    
    assert(type(p_config['dpi']) is int), \
    "dpi in config must be given as int"
    
    #plotname
    assert(type(p_config['plot_title']) is str), \
    "plot_name in config must be given as a string"
    
    assert(type(p_config['plot_title_fontsize']) is int or float), \
    "plot_name_fontsize in config must be given as int or float"
    
    assert(type(p_config['spacing']) is int or float), \
    "spacing in config must be given as in or float"
    
    assert(type(p_config['cutoff_point']) is int or float), \
    "cutoff_point in config must be given as int or float"
    
    assert(type(p_config['cutoff_data']) is bool), \
    "cutoff_data in config must be given as True or False"
    
    assert(type(p_config['marker']) is str), \
    "marker in config must be given as str"
    
    assert(type(p_config['colors']) is list and \
           type(p_config['colors'][0]) is str and \
           len(p_config['colors']) == len(d_config['subsamples'])), \
    "colors in config must be given as a list of str, with one color for each subsample"
    
    assert(type(p_config['marker_size']) is int or float), \
    "marker_size in config must be given as int or float"
    
    assert(type(p_config['gene_arrows']) is bool), \
    "gene_arrows in config must be given as True or False"
    
    assert(type(p_config['chr_name_rotate']) is int or float), \
    "chr_name_rotate in config must be given as int or float, 0 = horizontal, 90 = vertical"
    
    assert(type(p_config['chr_name_fontsize']) is int or float), \
    "chr_name_fontsize in config must be given as int or float"
    
    assert(type(p_config['y_label_fontsize']) is int or float), \
    "y_label_fontsize in config must be given as int or float"
    
    
    if verbose:
        print()
        print("-----------Plot Settings------------")
        print()
        print("    filename:", p_config['filename'])
        print("  plot title:", p_config['plot_title'])
        print()
        print("cutoff point:", p_config['cutoff_point'])
        print("     y limit:", p_config['y_limit'])
        print()


    return

#-------------------------------------------------------------------------------------------

def get_wanted_files(config, files): 
    
    #return files needed for this plot 
    
    wanted_files = []
    for ss in config['subsamples']:
        wf = files[str(ss)][config['control']][config['type']]
        wanted_files.append(wf)
    
    return wanted_files

#-------------------------------------------------------------------------------------------

def make_plot(config, wanted_files, save=False, show=True, verbose=True):
    
    #make plot
    
    d_config = config['data']
    p_config = config['plot']
    
    #get list of chromosomes to plot
    chr_to_plot = d_config['chromosomes']
    if chr_to_plot == ['all']:
        chr_to_plot = chr_names
        
    #start plot
    fig, ax = plt.subplots(figsize=p_config['figure_size'])
    
    curr_max = 0
    
    #go through every subsample...
    for i in range(len(wanted_files)):
        
        file = wanted_files[i]
        color = p_config['colors'][i]
        ss_num = d_config['subsamples'][i] #matching color to subsample
        
        #read data from file
        data = pd.read_csv(file, sep="\t")
        
        #get sample data
        sample = d_config['sample']
        all_sample_data = list(data[sample + "_" + d_config['type']])
        
        
        chr_data = []     #extract data for each chromosome and save as y values
        chr_x = []        #match with incremental x values w/ spacing inbetween
        chr_label_x = []  #save midpoints for labelling chromosome names 
        chr_tick_x = []   #beginning and endpoints for ticks 
        curr_x = 0
        for chro in chr_to_plot:
            
            if verbose:
                print(chro, ': getting data')
            start = chr_locs[chro][0]
            end = chr_locs[chro][1]
            chr_data += all_sample_data[start:end + 1] #get data by start/end indices
            
            if p_config['cutoff_data']:
                if verbose:
                    print(chro, ': getting cutoff mask')
                chr_data = sliding_window_mask(chr_data, p_config['cutoff_point'], size=3)
            
            chr_size = end - start + 1 #size of chromosome (in bin units)
            
            x = list(range(curr_x, curr_x + chr_size)) #make xvals 
            chr_x += x #store xvals 
            
            chr_label_x.append(curr_x + chr_size/2) #center
            chr_tick_x.append(curr_x)
            chr_tick_x.append(curr_x + chr_size)
            
            curr_x = x[-1] + p_config['spacing'] #account for spacing 
          
        
        #plot data
        if verbose:
            print("plotting...")
        ax.scatter(chr_x, chr_data, marker=p_config['marker'], color=color, s=p_config['marker_size'], 
                  label=sample + ", subsampling = " + str(ss_num) + " million")
        
        maxi = max(chr_data)
        if maxi > curr_max:
            curr_max = maxi
            
    
    #configure remaining plot stuff
    
    ax.margins(x=0) #remove whitespace around x axis 
    ax.set_ylim(bottom=p_config['cutoff_point'])
   
    if verbose:
        print("max fold enrichment is", curr_max)
    if p_config['y_limit'] is not None:
        ax.set_ylim(top=p_config['y_limit'])
    
    if p_config['gene_arrows']:
        plot_gene_arrows(ax, chr_tick_x, chr_to_plot, p_config['cutoff_point'], p_config['figure_size'])
    
    #setting ticks
    if p_config['chr_start_end_lines']:
        ax.set_xticks(chr_tick_x, labels=[])
        ax.tick_params('x', length=40)
    else:
        ax.set_xticks([], labels=[])
    
    #setting chromosome name labels
    sec = ax.secondary_xaxis(location=0)
    sec.set_xticks(chr_label_x, labels=chr_to_plot, rotation=p_config['chr_name_rotate'], 
                  fontsize=p_config['chr_name_fontsize']) #label chromosomes
    if p_config['chr_start_end_lines']:
        sec.tick_params('x', length=0)
    sec.tick_params('x', pad=27)
    
    #y label
    if d_config['type'] == 'avg':
        y_label = 'Average'
    elif d_config['type'] == 'fold':
        y_label = "Fold Enrichment"
    ax.set_ylabel(y_label, fontsize=p_config['y_label_fontsize'])

    #plot name
    plt.title(p_config['plot_title'], fontsize=p_config['plot_title_fontsize'], y=1.02)
    ax.legend(fontsize=p_config['legend_fontsize'])
    
    #save plot inside of plots folder
    if save:
        if verbose:
            print('saving...')
        plt.savefig('plots/' + p_config['filename'], dpi=p_config['dpi'], bbox_inches='tight')
        #also save config
        with open('plots/' + p_config['filename'].strip(".png") + "_config.json", "w") as f:
            json.dump(config, f, indent=2)
    #show plot
    if show:
        plt.show()
        
    return

#-------------------------------------------------------------------------------------------

def sliding_window_mask(data, cutoff, size=3):
    
    #for window of size centered on datapoint - if all datapoints in window are > cutoff, retain datapoint
    #else, remove datapoint (set to 0)
    
    lookback = int(size/2) #number of points to look on either side of current datapoint
    
    masked = data.copy()
    for i in range(lookback): #set beginning and endpoints to 0
        masked[i] = 0
        masked[-i-1] = 0
    
    for i in range(lookback, len(data) - lookback):
        window = data[i-lookback: i+lookback+1]
        
        if not all(x >= cutoff for x in window):
            masked[i] = 0
        
    return masked

#-------------------------------------------------------------------------------------------

def plot_gene_arrows(ax, chr_ticks, chr_list, cutoff, figsize):
    
    #for scaling arrow size
    x_width = ax.get_xlim()[1]
    y_height = ax.get_ylim()[1] - ax.get_ylim()[0]
    x_scale = x_width / figsize[0]
    y_scale = y_height / figsize[1]
    
    #print(x_scale, y_scale)
    
    arrow_y = cutoff - y_scale / 5
    arrow_y_2 = cutoff - y_scale / 3
    arrow_head_width = y_scale / 10
    arrow_head_length = x_scale / 10
    arrow_width = y_scale / 50
    
    chr_start = chr_ticks[0::2]
    
    gene_info = get_gene_info(chr_list)
    
    for i in range(len(chr_list)):
        
        chr = chr_list[i]
        chr_start_x = chr_start[i]
        
        #plot based on file info 
        for gene in gene_info[chr]:
            
            start = chr_start_x + gene[0]
            length = gene[1]
            if gene[2] == 1: #forward
            
                ax.arrow(start, arrow_y, length, 0, color="k", clip_on=False, 
                         head_width=arrow_head_width, head_length = arrow_head_length, shape="right",
                         length_includes_head = True, width=arrow_width)
            elif gene[2] == -1: #backward
        
                ax.arrow(start, arrow_y_2, -length, 0, color="k", clip_on=False, 
                         head_width=arrow_head_width, head_length = arrow_head_length, shape="right",
                         length_includes_head = True, width=arrow_width)
    
    return

def get_gene_info(chr_list):
    
    #get gene info from file for each chr in list 
    
    #read files
    with open("ref/gene_locations.json", 'r') as f:
        gene_locs = json.load(f)
    with open("ref/pseudogene_locations.json", 'r') as f:
        pseudo_locs = json.load(f)
        
    gene_info = {}
    for chro in chr_list:
        gene_info[chro] = [] 
        all_genes = gene_locs[chro] + pseudo_locs[chro]
        for gene in all_genes:
        
            start_scaled = gene[0]/50
            end_scaled = gene[1]/50
            length = end_scaled - start_scaled
            strand = gene[2]
            
            #print(chro, start_scaled, end_scaled, strand)
            
            if strand == 1:
                gene_info[chro].append([start_scaled, length, strand])
            elif strand == -1:
                gene_info[chro].append([end_scaled, length, strand])
    
    return gene_info
        
#-------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                    prog='Plotter for Cut&Tag Sequencing, BamCoverage Fold Enrichment',
                    description='Reads config file and plots fold enrichment of specified chromosomes/samples.')
    
    parser.add_argument('-c', '--config_file', type=str, help="path to config file")
    parser.add_argument('-s', '--save', action='store_true', help='save plot and config file to plots folder')
    parser.add_argument('-q', '--quiet', action='store_true', help='suppress config parsing messages')
    parser.add_argument('-i', '--invisible', action='store_true', help='don\'t show plot')
    
    args = parser.parse_args()
    #print(args.config_file, args.save, args.quiet, args.invisible)
    
    if args.config_file is None:
        config_filename = 'config.json'
    else:
        config_filename = args.config_file
    save = args.save
    show = not (args.invisible)
    verbose = not (args.quiet)

    files, sample_names, chr_names, chr_locs = return_info() 
    config = parse_config(config_filename, sample_names, chr_names, verbose=verbose)
    wanted_files = get_wanted_files(config['data'], files)
    make_plot(config, wanted_files, save=save, show=show, verbose=verbose)
    