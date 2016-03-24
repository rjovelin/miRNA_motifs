# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:12:32 2016

@author: RJovelin
"""



# use this script to compare abundance of motifs upstream of TSS and 
# between mirna and TSS

from get_positions_from_mast import *
from mirna_data_for_meme import *


# use this function to get the TSS positions relative to mirna hairpin start
# the data uses TSS annotations for ncRNA from Kruesi et al 2013. Elife
def get_TSS_distances(TSS_file):
    '''
    (file) -> dict
    Return a dictionary with mirna name as key and the distance between TSS
    and mirna hairpin start
    '''
    
    # open file for reading
    infile = open(TSS_file)
    # skip header
    infile.readline()
    # create a dict of TSS positions
    TSS = {}
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get mirna name, distance pair
            if line[-1] != 'NA':
                TSS[line[2]] = int(line[-1])
    # close file
    infile.close()
    
    return TSS

# use this function to get the TSS positions in mirna upstream sequences    
def get_TSS_positions(TSS_file, upstream_file):
    '''
    (file, file) -> dict
    Take the file with TSS distances between the mirna hairpin and the TSS and 
    the file with mirna upstream sequences and return a dict with mirna name as
    key and a list with the position of the TSS in the upstream sequence and the
    length of the upstream sequence
    '''
    
    # get the distances between mirna and TSS
    TSS_dist = get_TSS_distances(TSS_file)
    
    # create a dict from fasta upstream sequences
    upstream_seqs = convert_fasta(upstream_file)
    
    # create a dict with mirna name and upstream sequence length
    seq_length = {}
    # loop over dict with mirna upstream sequences
    for mirna in upstream_seqs:
        # parse the name to get rid of the cel-suffix and the mirbase ID
        name = mirna[mirna.index('cel-')+4: mirna.index('_MI')]
        seq_length[name] = len(upstream_seqs[mirna])
    
    # create a dict to store the TSS positions
    TSS_pos = {}
    # loop over mirna with TSS
    for mirna in TSS_dist:
        # check that mirna has upstream seq
        if mirna in seq_length:
            # position = length - distance
            position = seq_length[mirna] - TSS_dist[mirna]
            # record only positions within the 1000 Kb upstream of mirna
            # mir-76 and mir-238 have TSS positions upstream the 1000 kb
            if position > 0:
                TSS_pos['cel-' + mirna] = [position, seq_length[mirna]]
            
    return TSS_pos
    
    
def sort_motifs_TSS(mast_output, TSS_file, upstream_file):
    '''
    (file, file, file) -> dict
    Return a dictionary with mirna name and list with motif
    counts upstream of TSS and with motif counts between TSS and mirna
    and the same counts normalized by sequence length
    '''
    
    # get the TSS positions
    TSS_pos = get_TSS_positions(TSS_file, upstream_file)
    
    # get motif positions
    motif_pos = motif_position_from_mast(mast_output)[1]

    # create a dict with motif counts
    motifs = {}
    # loop over mirnas with motifs
    for mirna in motif_pos:
        # check if mirna has TSS annotation
        name = mirna[:mirna.index('_MI')]
        if name in TSS_pos:
            # initialise counts [before TSS, after TSS, before_normalized, after_normalized]
            motifs[mirna] = [0, 0, 0, 0]
            for i in motif_pos[mirna]:
                if i <= TSS_pos[name][0]:
                    # motif upstream of TSS
                    motifs[mirna][0] += 1
                    # update normalized upstream count
                    motifs[mirna][2] += 1 / TSS_pos[name][0]
                elif i > TSS_pos[name][0]:
                    # motif downstream of TSS (between TSS and mirna)
                    motifs[mirna][1] += 1
                    # update normalized downstream count
                    motifs[mirna][3] += 1 / (TSS_pos[name][1] - TSS_pos[name][0])
    
    return motifs


if __name__ == '__main__':
    import os
    import numpy as np
    from scipy import stats
    # create a dict with motif correspondence
    motif_code = {'motif2': 'm1', 'motif5': 'm2', 'motif4': 'm3', 'motif3': 'm4', 'motif1': 'm5'}
    # open file for writing
    outputfile = input('enter the outputfile name:' )
    newfile = open(outputfile, 'w')
    # write header to file
    header = '\t'.join(['motif', 'N_mirnas', 'N_motifs_before_TSS', 'N_motifs_after_TSS',
                        'mean_before_TSS', 'mean_after_TSS',  'P_Wilcoxon_paired',
                        'mean_normalized_before_TSS', 'mean_normalized_after_TSS', 'P_Wilcoxon_paired'])
    newfile.write(header + '\n')                    
    files = [i for i in os.listdir() if '_mast.txt' in i]
    files.sort()
    # loop over file, get the number of motifs before TSS and after TSS
    for i in files:
        motifs = sort_motifs_TSS(i, 'Kruesi_media3.txt', 'celegans_mirna_upstream_1000.txt')
        # extract motif name
        motif_name = i[i.index('motif'): i.index('_', i.index('motif'))]
        after, before, after_norm, before_norm = [], [], [] ,[]
        for j in motifs:
            before.append(motifs[j][0])
            after.append(motifs[j][1])
            before_norm.append(motifs[j][2])
            after_norm.append(motifs[j][3])
        newfile.write('\t'.join([motif_code[motif_name], str(len(before)),  str(sum(before)), str(sum(after)),
                                 str(np.mean(before)), str(np.mean(after)),
                                 str(stats.wilcoxon(before, after)[1]),
                                 str(np.mean(before_norm)), str(np.mean(after_norm)),
                                 str(stats.wilcoxon(before_norm, after_norm)[1])]) + '\n')
    newfile.close()
    
