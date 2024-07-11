from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

import argparse
import pathlib


"""
Takes an input genome and creates three subgenomes:

1) genes extracted and separated by 200 bases of Ns
2) unspliced transcript fragments, +/- [fragment_size] nts around the beginning of the Open Reading Frame (ORF) of each gene, extracted and separated by 200 bases of Ns
3) spliced transcript fragments, [fragment-size] nts at the beginning of the ORF of each gene, extracted and separated by 200 bases of Ns






To use in command line, you must input arguments for the .gff file and .fasta file. 
    
    Ex.    python get_genes_unspliced_spliced.py Tbrucei_annotations.gff Tbrucei_genome.fasta
    


An optional argument --outfile allows you to name the output file prefix.

    Ex.    python get_genes_unspliced_spliced.py Tbrucei_annotations.gff Tbrucei_genome.fasta --outfile Tbrucei
    
With the above call, the output files will be named TBrucei_GENES.gff, TBrucei_UNSPLICED.gff, and TBrucei_SPLICED_147.gff (and similarly for .fasta). Otherwise, the output file will default to the input filenames: TBrucei_annotations_GENES.gff, TBrucei_genome_GENES.fasta, etc. 



A second optional argument --fsize allows you to define the fragment size. 

    Ex.    python get_genes_unspliced_spliced.py Tbrucei_annotations.gff Tbrucei_genome.fasta --outfile Tbrucei --fsize 150
    
If --fsize is not specified, the default fragment size is 147. 

"""


#------------------------------------------------------------------------------------------------------

def get_slice(seq, start, end):
    #returns slice specified by DNA-standard indexing 
    return seq[start-1:end]

#------------------------------------------------------------------------------------------------------

def get_all_info(gff_file, fasta_file, sequence_ids):
    
    #given input gff and fasta file, extract mRNA and pseudogenic_transcript entries
    #output gff and sequence information in dictionaries
    
    
    #gff info saved to dictionary wanted_entries = { sequence_id:seqRecord }
    wanted_entries = {}

    for id in sequence_ids:

        limit_info = dict(gff_id=[id], gff_type=["mRNA", "pseudogenic_transcript"])
        wanted_entries[id] = []
        in_handle = open(gff_file)
        for rec in GFF.parse(in_handle, limit_info=limit_info):
            #get rid of annotations...
            rec.annotations = {}
            wanted_entries[id] = rec

        in_handle.close()
        
        
    #get chromosome sequences from fasta
    chr_seqs_list = list(SeqIO.parse(fasta_file, "fasta"))

    #turn it into a dictionary for ease of use...
    chr_seqs = {}
    for seq in chr_seqs_list:
        seq.description = ""
        chr_seqs[seq.id] = seq
        
        
    return wanted_entries, chr_seqs


#------------------------------------------------------------------------------------------------------

def get_gene_sequence_slices(wanted_entries, chr_seqs, sequence_ids):
    
    #for each id, get exact gene sequence slices
    
    #store in a dictionary, matchingg list index to feature list index 
    
    
    #make dictionary
    exact_seq_slices = {}
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            exact_seq_slices[id] = [""]*len(wanted_entries[id].features)
            
    #for each id...
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            #for each feature of this chromosome...
            for i in range(len(wanted_entries[id].features)):
                feature = wanted_entries[id].features[i]
                start = feature.location.start
                end = feature.location.end
                
                seq = get_slice(chr_seqs[id].seq, start + 1, end)
                
                exact_seq_slices[id][i] = seq 
                
    return exact_seq_slices

#------------------------------------------------------------------------------------------------------

def get_unspliced_slices(wanted_entries, chr_seqs, sequence_ids, fragment_size):
    
    #get fragments around beginning of genes
    
    #make dictionary
    unspliced_slices = {}
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            unspliced_slices[id] = [""]*len(wanted_entries[id].features)
            
    #for each id...
    for id in sequence_ids:
        #for each feature...
        if type(wanted_entries[id]) != list:
            for i in range(len(wanted_entries[id].features)):

                feature = wanted_entries[id].features[i]
                start = feature.location.start
                end = feature.location.end 
                strand = feature.location.strand 

                
                #forward
                if strand == 1:
                    
                    if start >= fragment_size:
                        seq = get_slice(chr_seqs[id].seq, 
                                    start - (fragment_size - 1), start + fragment_size) #beginning shifted by 1
                    
                    #check for not enough sequence at beginning, pad with Ns
                    elif start < fragment_size: 
                        partial_seq = get_slice(chr_seqs[id].seq, 
                                               1, start + fragment_size)
                        
                        print(len(partial_seq))
                        pad = "N" * (fragment_size*2 - len(partial_seq))
                        seq = pad + partial_seq
                        print(len(seq), 'start padded')
                        print(id, i)
                    
                    
                    
                #backward
                elif strand == -1:
                    seq = get_slice(chr_seqs[id].seq, 
                                    end - (fragment_size - 1), end + fragment_size) #beginning shifted by 1
                    
                    #check for not enough sequence at end, pad with Ns
                    if len(seq) < fragment_size*2: 
                        print(len(seq))
                        pad = "N" * (fragment_size*2 - len(seq))
                        seq = seq + pad
                        print(len(seq), 'end padded')
                        print(id, i)
                
                unspliced_slices[id][i] = seq
                
    return unspliced_slices

#------------------------------------------------------------------------------------------------------

def get_spliced_slices(wanted_entries, chr_seqs, sequence_ids, fragment_size, 
                       splice_leader_forward, splice_leader_backward):
    
    #get fragment around gene beginning and add splice leader to beginning of fragment
    
    
    #make dictionary
    spliced_slices = {}
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            spliced_slices[id] = [""]*len(wanted_entries[id].features)
            
            
    #for each id...
    for id in sequence_ids:
        #for each feature...
        if type(wanted_entries[id]) != list:
            for i in range(len(wanted_entries[id].features)):

                feature = wanted_entries[id].features[i]
                start = feature.location.start
                end = feature.location.end 
                strand = feature.location.strand 

                if strand == 1:
                    seq = get_slice(chr_seqs[id].seq, start+1, start + fragment_size)
                    seq = splice_leader_forward + seq

                elif strand == -1:
                    seq = get_slice(chr_seqs[id].seq, end - (fragment_size-1), end)
                    seq = seq + splice_leader_backward

                spliced_slices[id][i] = seq
                
    return spliced_slices

#------------------------------------------------------------------------------------------------------

def make_200n_files_just_genes(chr_seqs, gene_info, piece_seqs, gff_file, fasta_file, outfile_prefix=None):
    
    #make gff and fasta file with gene sequences separated by spacers of 200N
    
    #chr_seqs is sequences from fasta 
    #gene_info should be wanted_entries from gff 
    #pieces should be dictionary of specific wanted parts of the sequence 
    
    #both indexed by chromosome name
    
    #gene_info and piece_seqs should be dicts of lists with matching lengths 
    
    ns = "N"*200
    spacer = Seq(ns)
    
    for id in sequence_ids:
        
        if type(wanted_entries[id]) != list:
            chr_seq = chr_seqs[id]
            gff_info = gene_info[id]
            pieces = piece_seqs[id]

            #first, combine all pieces w/ spacer of 200N
            combined_pieces = spacer.join(pieces)
            combined_pieces = combined_pieces.join((spacer, spacer))

            #then, replace sequence in fasta info 
            chr_seq.seq = combined_pieces

            #then, put sequence in gff info 
            gff_info.seq = combined_pieces

            #then, update indices in gff info
            sum = 200

            for i in range(len(gff_info.features)):
                feature = gff_info.features[i]
                size = feature.location.end - feature.location.start
                newstart = sum + 1 - 1 #????????
                newend = sum + size
                feature.location = FeatureLocation(newstart, newend, strand=feature.location.strand)

                sum += size + 200
            
            
    #save to fasta and gff
    
    if outfile_prefix is None:
        fasta_out_name = fasta_file.split(".fasta")[0] + "_GENES.fasta"
        gff3_out_name = gff_file.split(".gff")[0] + "_GENES.gff3"
        
    else:
        fasta_out_name = outfile_prefix + "_GENES.fasta"
        gff3_out_name = outfile_prefix + "_GENES.gff3"
        
    fasta_recs = list(chr_seqs.values())
    gff3_recs = list(gene_info.values())
    
    gff3_recs = [x for x in gff3_recs if type(x) != list]

    with open(fasta_out_name, "w") as out:
        SeqIO.write(fasta_recs, out, "fasta")
        
    with open(gff3_out_name, "w") as out:
        GFF.write(gff3_recs, out)
       
    
    
    return



#------------------------------------------------------------------------------------------------------

def make_200n_files_unspliced(chr_seqs, gene_info, piece_seqs, fragment_size, 
                              gff_file, fasta_file, outfile_prefix=None):
    
    #make gff and fasta file with unspliced sequences separated by 200N 
    
    #chr_seqs is sequences from fasta 
    #gene_info should be wanted_entries
    #pieces should be dictionary of specific wanted parts of the sequence 
    
    #both indexed by chromosome name
    
    #gene_info and piece_seqs should be dicts of lists with matching lengths 
    
    ns = "N"*200
    spacer = Seq(ns)
    
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            chr_seq = chr_seqs[id]
            gff_info = gene_info[id]
            pieces = piece_seqs[id]

            #first, combine all pieces w/ spacer of 200N
            combined_pieces = spacer.join(pieces)
            combined_pieces = combined_pieces.join((spacer, spacer))

            #then, replace sequence in fasta info 
            chr_seq.seq = combined_pieces

            #then, put sequence in gff info 
            gff_info.seq = combined_pieces

            sum = 200

            for i in range(len(gff_info.features)):
                feature = gff_info.features[i]
                newstart = sum + 1 - 1 
                newend = sum + (fragment_size*2)
                sum += 200 + (fragment_size*2)
                feature.location = FeatureLocation(newstart, newend, strand=feature.location.strand)
        
            
           
    #save to fasta and gff
    if outfile_prefix is None:
        fasta_out_name = fasta_file.split(".fasta")[0] + "_UNSPLICED_" + str(fragment_size) + ".fasta"
        gff3_out_name = gff_file.split(".gff")[0] + "_UNSPLICED_" + str(fragment_size) + ".gff3"
        
    else:
        fasta_out_name = outfile_prefix + "_UNSPLICED_" + str(fragment_size) + ".fasta"
        gff3_out_name = outfile_prefix + "_UNSPLICED_" + str(fragment_size) + ".gff3"
      
    
    fasta_recs = list(chr_seqs.values())
    gff3_recs = list(gene_info.values())
    
    gff3_recs = [x for x in gff3_recs if type(x) != list]

    with open(fasta_out_name, "w") as out:
        SeqIO.write(fasta_recs, out, "fasta")
        
    with open(gff3_out_name, "w") as out:
        GFF.write(gff3_recs, out)
        

    return


#------------------------------------------------------------------------------------------------------


def make_200n_files_spliced(chr_seqs, gene_info, piece_seqs, fragment_size, splice_leader_forward,
                             gff_file, fasta_file, outfile_prefix = None):
    
    
    #make fasta and gff of spliced fragments, separated by 200N
    
    
    #chr_seqs is sequences from fasta 
    #gene_info should be wanted_entries
    #pieces should be dictionary of specific wanted parts of the sequence 
    
    #both indexed by chromosome name
    
    #gene_info and piece_seqs should be dicts of lists with matching lengths 
    
    ns = "N"*200
    spacer = Seq(ns)
    
    for id in sequence_ids:
        if type(wanted_entries[id]) != list:
            chr_seq = chr_seqs[id]
            gff_info = gene_info[id]
            pieces = piece_seqs[id]

            #first, combine all pieces w/ spacer of 200N
            combined_pieces = spacer.join(pieces)
            combined_pieces = combined_pieces.join((spacer, spacer))

            #then, replace sequence in fasta info 
            chr_seq.seq = combined_pieces

            #then, put sequence in gff info 
            gff_info.seq = combined_pieces

            sum = 200

            for i in range(len(gff_info.features)):
                feature = gff_info.features[i]
                newstart = sum
                newend = sum + fragment_size + len(splice_leader_forward)
                sum += 200 + fragment_size + len(splice_leader_forward)
                feature.location = FeatureLocation(newstart, newend, strand=feature.location.strand)


           
    #save to fasta and gff
    if outfile_prefix is None:
        fasta_out_name = fasta_file.split(".fasta")[0] + "_SPLICED_" + str(fragment_size) + ".fasta"
        gff3_out_name = gff_file.split(".gff")[0] + "_SPLICED_" + str(fragment_size) + ".gff3"
        
    else:
        fasta_out_name = outfile_prefix + "_SPLICED_" + str(fragment_size) + ".fasta"
        gff3_out_name = outfile_prefix + "_SPLICED_" + str(fragment_size) + ".gff3"
    
    fasta_recs = list(chr_seqs.values())
    gff3_recs = list(gene_info.values())
    
    gff3_recs = [x for x in gff3_recs if type(x) != list]

    with open(fasta_out_name, "w") as out:
        SeqIO.write(fasta_recs, out, "fasta")
        
    with open(gff3_out_name, "w") as out:
        GFF.write(gff3_recs, out)
        
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

        
if __name__ == "__main__": 
    
    
    #parse command line arguments
    parser = argparse.ArgumentParser(description="Creates three files from a given genome: genes," + \
                                     " unspliced gene fragments, and spliced gene fragments.") 

    parser.add_argument('gff_file', type=str, help="input gff file")
    parser.add_argument('fasta_file', type=str, help="input fasta file")
    parser.add_argument('--outfile', type=str, help="optional outfile name, default labels appended to infile names") 
    parser.add_argument('--fsize', type=int, help="optional fragment size, default 147") 

    args = parser.parse_args()

    #check proper file extensions 
    gff_input_suff = pathlib.Path(args.gff_file).suffix
    if ".gff" not in gff_input_suff:
        parser.error("gff_file must be a .gff file") 

    fasta_input_suff = pathlib.Path(args.fasta_file).suffix
    if ".fasta" not in fasta_input_suff:
        parser.error("fasta_file must be a .fasta file") 

#------------------------------------------------------------------------------------------------------
            
    #set variables 
    
    gff_file = args.gff_file
    fasta_file = args.fasta_file

    if args.fsize is not None:
        frag_size = args.fsize
    else:
        frag_size = 147
        
    outfile_prefix = args.outfile
    
    
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

    #edit here if you want to change splice leader
    
    splice_leader_forward = "AACGCTATTATTAGAACAGTTTCTGTACTATATTG"
    splice_leader_backward = "CAATATAGTACAGAAACTGTTCTAATAATAGCGTT"
    
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
    
    #get sequence ID names from GFF3 file (from the header)
    examiner = GFFExaminer()
    in_handle = open(gff_file)
    limits = examiner.available_limits(in_handle)
    in_handle.close()

    sequence_ids = list(limits['gff_id'].keys())
    sequence_ids = [x[0] for x in sequence_ids]
    
    
    #get info from genome
    print("retrieving genome info...")
    wanted_entries, chr_seqs = get_all_info(gff_file, fasta_file, sequence_ids)
    print("info retrieved")
    
    #extract genes / unspliced fragments / spliced fragments
    
    print("extracting gene sequences...")
    exact_seq_slices = get_gene_sequence_slices(wanted_entries, chr_seqs, sequence_ids)
    print("extracting unspliced slices...") 
    unspliced_slices = get_unspliced_slices(wanted_entries, chr_seqs, sequence_ids, frag_size)
    print("extracting spliced slices...") 
    spliced_slices = get_spliced_slices(wanted_entries, chr_seqs, sequence_ids, frag_size,
                                   splice_leader_forward, splice_leader_backward)
    
    #export to files
    print("exporting to files...") 
    make_200n_files_just_genes(chr_seqs, wanted_entries, exact_seq_slices, 
                           gff_file, fasta_file, outfile_prefix=outfile_prefix)
    make_200n_files_unspliced(chr_seqs, wanted_entries, unspliced_slices, frag_size, 
                          gff_file, fasta_file, outfile_prefix=outfile_prefix)
    make_200n_files_spliced(chr_seqs, wanted_entries, spliced_slices, frag_size, splice_leader_forward,
                       gff_file, fasta_file, outfile_prefix=outfile_prefix)
    
    return