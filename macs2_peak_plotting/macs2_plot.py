from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
import json
import argparse

#-------------------------------------------------------------------------------------------

def return_info():
    
    #look for info and return files, sample names, chromosomes 
    info_file = 'ref/info.json'
    
    with open(info_file, 'r') as f:
        info = json.load(f)
        
    sample_names = info['sample_names']
    chr_names = info['chr_names']
    
    with open('ref/chr_lengths.json', 'r') as f:
        lengths = json.load(f)
    
    return sample_names, chr_names, lengths

#-------------------------------------------------------------------------------------------

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
    
    #sample
    assert(d_config['sample'] in sample_names), \
    "sample " + str(d_config['sample']) + " listed in config is unknown"
    
    #chr
    for chr_name in d_config['chromosomes']:
        assert(chr_name in chr_names or chr_name == "all"), \
        "chromosome " + str(chr_name) + " listed in config is unknown"

    
    if verbose:
        #print data settings to console    
        print()
        print("-----------Data Settings------------")
        print()
        print("      sample:", d_config['sample'])
        print()
        print(" chromosomes:")
        print("\n".join("\t    " + str(x) for x in d_config['chromosomes']))
        print()
    
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
    
    assert(type(p_config['color']) is str), \
    "colors in config must be given as a list of str, with one color for each subsample"
    
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
        print("     y limit:", p_config['y_limit'])
        print()

    return


#-------------------------------------------------------------------------------------------

def get_wanted_file(config): 
    
    #return files needed for this plot 
    
    files = glob.glob('data/*.bed')
    
    sample = config['data']['sample']
    
    sample_file = None
    i = 0
    
    while sample_file is None and i < len(files): 
        
        file = files[i]
        file_sample = file.split("[")[-1].split("]")[0].split("_")[-1]
        
        if file_sample == sample:
            sample_file = file
            
        i+=1
    
    
    return sample_file

#-------------------------------------------------------------------------------------------

def make_plot_dots(config, file, lengths, save=False, show=True, verbose=True):
    
    #make plot
    
    d_config = config['data']
    p_config = config['plot']
    
    #get list of chromosomes to plot
    chr_to_plot = d_config['chromosomes']
    if chr_to_plot == ['all']:
        chr_to_plot = chr_names
    
    #start plot
    fig, ax = plt.subplots(figsize=p_config['figure_size'])

    #get data from file
    peak_data = get_peak_data(file, chr_to_plot)


    #chr_data = []     #extract data for each chromosome and save as y values
    #chr_x = []        #match with incremental x values w/ spacing inbetween
    chr_label_x = []  #save midpoints for labelling chromosome names 
    chr_tick_x = []   #beginning and endpoints for ticks 
    curr_x = 0
    
    peak_vals = []
    x_vals = []
    
    for chro in chr_to_plot:
        
        chr_peak_data = peak_data[chro]
        chr_length = lengths[chro]

        if verbose:
            print(chro, ': getting peak data (locations', curr_x, "to", curr_x + chr_length, ")")
        
        for peak in chr_peak_data:
            
            peak_val = peak[2]
            peak_start = peak[0]
            peak_end = peak[1]
            
            peak_length = peak_end - peak_start
        
            x_val = curr_x + peak_start + peak_length / 2
            
            peak_vals.append(peak_val)
            x_vals.append(x_val)
            
            if verbose: 
                print("\t peak at location", x_val, "with value of", peak_val)
    

        chr_label_x.append(curr_x + chr_length/2) #center
        chr_tick_x.append(curr_x)
        chr_tick_x.append(curr_x + chr_length)

        curr_x += chr_length + p_config['spacing']*50 #account for spacing 


    #plot data
    if verbose:
        print("plotting...")
    ax.scatter(x_vals, peak_vals, color=p_config['color'], 
              label=d_config['sample'] + " narrow peaks", s = p_config['marker_size'])
    #ax.scatter(chr_x, chr_data, marker=p_config['marker'], color=color, s=p_config['marker_size'], 
              #label=sample + ", subsampling = " + str(ss_num) + " million")
            
    
    #configure remaining plot stuff
    
    ax.set_xlim(left=0, right=chr_tick_x[-1]) #remove whitespace around x axis 
    ax.set_ylim(bottom=0)
    
    if p_config['y_limit'] is not None:
        ax.set_ylim(top=p_config['y_limit'])
    
    if p_config['gene_arrows']:
        plot_gene_arrows(ax, chr_tick_x, chr_to_plot, 0, p_config['figure_size'])
    
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
    
    
    ax.set_ylabel("Narrow Peaks", fontsize=p_config['y_label_fontsize'])

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
            
            start = gene[0]
            end = gene[1]
            length = end - start
            strand = gene[2]
            
            if strand == 1:
                gene_info[chro].append([start, length, strand])
            elif strand == -1:
                gene_info[chro].append([end, length, strand])
    
    
    return gene_info
        
#-------------------------------------------------------------------------------------------

def get_peak_data(file, chrs_wanted):
    
    #return data as dictionary
    #peak_data['chromosome'] = [start, end, peak]
    
    data = pd.read_csv(file, sep='\t', header=None, usecols=[0, 1, 2, 4])
    
    peak_data = {}
    
    for chro in chrs_wanted:
        peak_data[chro] = []
     
    for i in range(len(data)):
        
        chro = data[0].iloc[i]
        if chro in chrs_wanted:
            
            start = data[1].iloc[i]
            end = data[2].iloc[i]
            peak = data[4].iloc[i]
            
            peak_data[chro].append([int(start), int(end), float(peak)])
            
    return peak_data
    
#-------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                    prog='Plotter for Cut&Tag Sequencing, MACS2 Narrow Peak',
                    description='Reads config file and plots narrow peaks of specified chromosomes/samples.')
    
    parser.add_argument('-c', '--config_file', type=str, help="path to config file")
    parser.add_argument('-s', '--save', action='store_true', help='save plot and config file to plots folder')
    parser.add_argument('-q', '--quiet', action='store_true', help='suppress config parsing messages')
    parser.add_argument('-i', '--invisible', action='store_true', help='don\'t show plot')
    
    args = parser.parse_args()
    #print(args.config_file, args.save, args.quiet, args.invisible)
    
    if args.config_file is None:
        config_filename = 'macs2_config.json'
    else:
        config_filename = args.config_file
    save = args.save
    show = not (args.invisible)
    verbose = not (args.quiet)
    
    sample_names, chr_names, lengths = return_info()
    config = parse_config(config_filename, sample_names, chr_names, verbose=verbose)
    file = get_wanted_file(config)
    make_plot_dots(config, file, lengths, save=save, verbose=verbose)
        
        