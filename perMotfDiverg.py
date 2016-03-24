from seqDiverg import *
import random
from get_positions_from_mast import motif_position_from_mast
from mirna_data_for_meme import convert_fasta

def randomization_K_motif(Length, N_repeats):
    '''
    (int, int) -> dict
    Compute divergence at a random motif of size Length taken from a random alignment file N_repeats times
    '''

        

    # dictionnary of filenames to draw from

    filenames = {0: 'mir-2225.txt', 1: 'mir-2226.txt', 2: 'mir-2227.txt', 3: 'mir-2231.txt', 4: 'mir-4813.txt', 5: 'mir-7579.txt',
                 6: 'mir-7580-1.txt', 7: 'mir-7581.txt', 8: 'mir-7582.txt', 9: 'mir-7584.txt', 10: 'mir-7585.txt',
                 11: 'mir-7586.txt', 12: 'mir-7587.txt', 13: 'mir-7588.txt', 14: 'mir-7589.txt', 15: 'mir-7591.txt', 16: 'mir-7593.txt',
                 17: 'mir-7594a.txt', 18: 'mir-7594b.txt', 19: 'mir-7596.txt', 20: 'mir-86.txt', 21: 'mir-87.txt', 22: 'mir-90a.txt',
                 23: 'mir-90b.txt', 24: 'mir-124.txt', 25: 'mir-124b.txt', 26: 'mir-124c.txt', 27: 'mir-228.txt', 28: 'mir-230.txt',
                 29: 'mir-231.txt', 30: 'mir-232-2.txt', 31: 'mir-233.txt', 32: 'mir-234.txt', 33: 'mir-235b.txt', 34: 'mir-236.txt',
                 35: 'mir-237.txt', 36: 'mir-238.txt', 37: 'mir-239a.txt', 38: 'mir-239b.txt', 39: 'mir-240.txt', 40: 'mir-241.txt',
                 41: 'mir-242.txt', 42: 'mir-244.txt', 43: 'mir-245.txt', 44: 'mir-246.txt', 45: 'mir-247.txt', 46: 'mir-82.txt',
                 47: 'mir-83.txt', 48: 'mir-84.txt', 49: 'mir-85.txt', 50: 'mir-52.txt', 51: 'mir-54d.txt', 52: 'mir-57.txt',
                  53: 'mir-58.txt', 54: 'mir-58b.txt', 55: 'mir-60.txt', 56: 'mir-61.txt', 57: 'mir-62.txt', 58: 'mir-70.txt',
                  59: 'mir-71.txt', 60: 'mir-72.txt', 61: 'mir-73a.txt', 62: 'mir-73b.txt', 63: 'mir-75.txt', 64: 'mir-76.txt',
                  65: 'mir-77-1.txt', 66: 'mir-79.txt', 67: 'mir-80.txt', 68: 'mir-81.txt', 69: 'mir-34.txt', 70: 'mir-35c-1.txt',
                  71: 'mir-36.txt', 72: 'mir-42.txt', 73: 'mir-45-1.txt', 74: 'mir-46.txt', 75: 'mir-47.txt', 76: 'mir-48.txt',
                  77: 'mir-49.txt', 78: 'mir-50.txt', 79: 'mir-51.txt', 80: 'lsy-6.txt', 81: 'mir-1.txt', 82: 'mir-2.txt', 83: 'let-7.txt',
                  84: 'lin-4.txt', 85: 'mir-784.txt', 86: 'mir-785b.txt', 87: 'mir-787.txt', 88: 'mir-788.txt', 89: 'mir-789.txt',
                  90: 'mir-790-2.txt', 91: 'mir-791.txt', 92: 'mir-792.txt', 93: 'mir-1817.txt', 94: 'mir-1822.txt', 95: 'mir-2214.txt',
                  96: 'mir-2222-2.txt', 97: 'mir-2223.txt', 98: 'mir-254.txt', 99: 'mir-255.txt', 100: 'mir-259.txt', 101: 'mir-268.txt',
                  102: 'mir-353.txt', 103: 'mir-354.txt', 104: 'mir-355.txt', 105: 'mir-356.txt', 106: 'mir-358.txt', 107: 'mir-359.txt',
                  108: 'mir-360.txt', 109: 'mir-392.txt', 110: 'mir-251.txt', 111: 'mir-252.txt', 112: 'mir-253.txt', 113: 'mir-248.txt', 114: 'mir-249.txt'}


    divergence = []
    

    # repeat N_repeats times
    while N_repeats > 0:
        # draw a random file from all upstream files
        file_number = random.randint(0, len(filenames)-1)
        myfile = open(filenames[file_number], 'r')
        # make a dictionnary with sequences in the file
        seq_ali = {}
        for line in myfile:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                seq_ali[line[1:]] = ""
                seq_name = line[1:]
            else:
                seq_ali[seq_name] += line
        myfile.close()
        genes = []
        for mirna in seq_ali:
            genes.append(mirna)
        # pick a random position in the sequence and grab a motif of size Length        
        random_position = random.randint(0, len(seq_ali[genes[0]])-1)
        motif1 = seq_ali[genes[0]][random_position: random_position + Length]
        motif2 = seq_ali[genes[1]][random_position: random_position + Length]

        K = compute_K(motif1, motif2)
        # accept K only if K is defined to make sure to get exactly N_repeats values
        # store the divergence values in a list
        if K != None:
            divergence.append(K)
            N_repeats -= 1
    
    return divergence
                        



def alignments_to_dictionnary():
    '''
    (list) -> dict
    Take a list of filenames of sequence aligments and return a dictionnary of dictionnaries
    wuth the cb name as key and dictionnary of gene name with their sequence as value: {cb_mirna: {cb_mirna: sequence, csp9_mirna: sequence}}
    '''

    filenames = ['mir-2225.txt', 'mir-2226.txt', 'mir-2227.txt', 'mir-2231.txt', 'mir-4813.txt',
                 'mir-7579.txt', 'mir-7580-1.txt', 'mir-7581.txt', 'mir-7582.txt', 'mir-7584.txt',
                 'mir-7585.txt', 'mir-7586.txt', 'mir-7587.txt', 'mir-7588.txt', 'mir-7589.txt', 'mir-7591.txt',
                 'mir-7593.txt', 'mir-7594a.txt', 'mir-7594b.txt', 'mir-7596.txt', 'mir-86.txt', 'mir-87.txt',
                 'mir-90a.txt', 'mir-90b.txt', 'mir-124.txt', 'mir-124b.txt', 'mir-124c.txt', 'mir-228.txt',
                 'mir-230.txt', 'mir-231.txt', 'mir-232-2.txt', 'mir-233.txt', 'mir-234.txt', 'mir-235b.txt',
                 'mir-236.txt', 'mir-237.txt', 'mir-238.txt', 'mir-239a.txt', 'mir-239b.txt', 'mir-240.txt',
                 'mir-241.txt', 'mir-242.txt', 'mir-244.txt', 'mir-245.txt', 'mir-246.txt', 'mir-247.txt', 'mir-82.txt',
                 'mir-83.txt', 'mir-84.txt', 'mir-85.txt', 'mir-52.txt', 'mir-54d.txt', 'mir-57.txt', 'mir-58.txt',
                 'mir-58b.txt', 'mir-60.txt', 'mir-61.txt', 'mir-62.txt', 'mir-70.txt', 'mir-71.txt', 'mir-72.txt',
                 'mir-73a.txt', 'mir-73b.txt', 'mir-75.txt', 'mir-76.txt', 'mir-77-1.txt', 'mir-79.txt', 'mir-80.txt',
                 'mir-81.txt', 'mir-34.txt', 'mir-35c-1.txt', 'mir-36.txt', 'mir-42.txt', 'mir-45-1.txt', 'mir-46.txt',
                 'mir-47.txt', 'mir-48.txt', 'mir-49.txt', 'mir-50.txt', 'mir-51.txt', 'lsy-6.txt', 'mir-1.txt', 'mir-2.txt',
                 'let-7.txt', 'lin-4.txt', 'mir-784.txt', 'mir-785b.txt', 'mir-787.txt', 'mir-788.txt', 'mir-789.txt', 'mir-790-2.txt',
                 'mir-791.txt', 'mir-792.txt', 'mir-1817.txt', 'mir-1822.txt', 'mir-2214.txt', 'mir-2222-2.txt', 'mir-2223.txt',
                 'mir-254.txt', 'mir-255.txt', 'mir-259.txt', 'mir-268.txt', 'mir-353.txt', 'mir-354.txt', 'mir-355.txt', 'mir-356.txt',
                 'mir-358.txt', 'mir-359.txt', 'mir-360.txt', 'mir-392.txt', 'mir-251.txt', 'mir-252.txt', 'mir-253.txt', 'mir-248.txt', 'mir-249.txt']

    alignments = {}
    for file in filenames:
        myfile = open(file, 'r')
        # make a dictionnary with sequences in the file
        seq_ali = {}
        for line in myfile:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                seq_ali[line[1:]] = ""
                seq_name = line[1:]
            else:
                seq_ali[seq_name] += line
        myfile.close()
        for mirna in seq_ali:
            if 'cbr' in mirna:
                alignments[mirna] = seq_ali

    return alignments
        


def check_seq_position(seq, seqali):
       
    gaps = 0

    for i in range(len(seqali)):
        if seqali[i] != '-':
            j = i - gaps
        else:
            gaps += 1
        print(i, j, seq[j], seqali[i], end = '\t')
        if seqali[i] == '-':
            print(seqali[i])
        else:
            print(seq[j] == seqali[i])    


def get_position_in_aligned_seq(position, seq, seqali):
    '''
    (int) -> int
    Take an position index in sequence seq and return the index corresponding
    to the same index in the aligned sequence seqali
    '''
    
    gaps = 0

    for i in range(len(seqali)):
        if seqali[i] != '-':
            j = i - gaps
        else:
            gaps += 1
        if j == position:
            assert seq[position] == seqali[i], 'no match'
            return i
            

def K_motif(mast_output, length, alignments, upstreams):
    '''
    (file, int, dict, file) -> 


    '''
    
    # get the positions from mast {mirna: [positions]}
    position_tuple = motif_position_from_mast(mast_output)
    motif_positions = position_tuple[1]
    
    # get the upstream sequences
    upstream_seqs = convert_fasta(upstreams)
     
    # create a dictionnary to store the K values of each motif for a given mirna
    mirna_K = {}
    
    # loop over mirnas with motifs
    for mirna in motif_positions:
        # check if mirna has ortholog 
        if mirna in alignments:
            # get cbr and csp9 mirna names
            mirnames = [i for i in alignments[mirna]]
            for i in mirnames:
                if i != mirna:
                    csp9 = i
            # loop over motif positions
            for position in motif_positions[mirna]:
                # get the corresponding position in the cbr aligned sequence
                j = get_position_in_aligned_seq(position, upstream_seqs[mirna], alignments[mirna][mirna])
                # get motif sequences in the aligned upstream sequences
                motif1 = alignments[mirna][mirna][j: j + length]
                motif2 = alignments[mirna][csp9][j: j + length]
                K = compute_K(motif1, motif2)
                if K != None:
                    if mirna in mirna_K:
                        mirna_K[mirna].append(K)
                    else:
                        mirna_K[mirna] = [K]
                
         
    return mirna_K



def count_substitutions_motif(mast_output, length, alignments, upstreams):
    '''
    (file, int, dict, file) -> dict
    '''
    
    # get the positions from mast {mirna: [positions]}
    position_tuple = motif_position_from_mast(mast_output)
    motif_positions = position_tuple[1]
    
    # get the upstream sequences
    upstream_seqs = convert_fasta(upstreams)
     
    # create a dictionnary to store the number of substitutions at each position
    # in the motifs
    mirna_subs = {}
    
    # loop over mirnas with motifs
    for mirna in motif_positions:
        # check if mirna has ortholog 
        if mirna in alignments:
            # get cbr and csp9 mirna names
            mirnames = [i for i in alignments[mirna]]
            for i in mirnames:
                if i != mirna:
                    csp9 = i
            # loop over motif positions
            for position in motif_positions[mirna]:
                # get the corresponding position in the cbr aligned sequence
                j = get_position_in_aligned_seq(position, upstream_seqs[mirna], alignments[mirna][mirna])
                motif1 = alignments[mirna][mirna][j: j + length].upper()
                motif2 = alignments[mirna][csp9][j: j + length].upper()
                # get motif positions in the aligned upstream sequences
                undefined = 0
                for k in range(len(motif1)):
                    # count gaps and Ns
                    if motif1[k] == '-' or motif1[k] == 'N' or motif2[k] == '-' or motif2[k] == 'N':
                        undefined += 1
                # consider only motifs with maximum 3 undefined positions
                if undefined <= 3:
                    for k in range(len(motif1)):
                        # check if nucleotide at given position are diferent
                        # do not count gaps and Ns
                        if motif1[k] != '-' and motif2[k] != '-':
                            if motif1[k] != 'N' and motif2[k] != 'N':
                                if motif1[k] != motif2[2]:
                                    # count 1 substitution
                                    if k in mirna_subs:
                                        mirna_subs[k] += 1
                                    else:
                                        mirna_subs[k] = 1
                    
    return mirna_subs
    
    

def randomization_to_file(random_K, outputfile):
    '''
    (dict, str) -> None
    Print the dictionnary of randomization of divergence K to a file outputfile
    '''


    myfile = open(outputfile, 'w')
    myfile.write('Repeat' + '\t' + 'K' +'\n')
    for i in range(len(random_K)):
        myfile.write(str(i) + '\t' + str(random_K[i]) + '\n')

    myfile.close()



def motif_K_to_file(motif_K, outputfile):
    '''
    (dict, str) -> None
    Print the dictionnary of divergence at a given motif to a file outputfile
    '''

    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'K' + '\n')
    for mirna in motif_K:
        for K in motif_K[mirna]:
            myfile.write(mirna + '\t' + str(K) + '\n')

    myfile.close()


def partition_mirna_withmotif_nomotif(motif_K_file, mirna_divergence, mirna_upstream, outputfile):
    '''
    (file, file, file, file) -> file
    Parition mirna divergence into 2 categories for mirnas in mirna_upstream that have or not a motif in motif_K_file 
    '''

    # make a set of mirnas that have a motif in motif_K_file
    motif_file = open(motif_K_file, 'r')
    mirna_with_motif = set()
    line = motif_file.readline()
    for line in motif_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0][:line[0].index('_')]
            mirna_with_motif.add(mirna)
    motif_file.close()

    # make a set of mirna that have no motif
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(mirna_upstream)
    mirna_no_motif = set()
    for mirna in upstreams:
        mirna_name = mirna[:mirna.index('_')]
        if not mirna_name in mirna_with_motif:
            mirna_no_motif.add(mirna_name)

    # make a dictionnary storing mirna divergence {mirna:[K_hairpin, K_mature, K_backbone]
    mirna_K_file = open(mirna_divergence, 'r')
    mirna_K = {}
    line = mirna_K_file.readline()
    for line in mirna_K_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            mirna_K[mirna] = line
            del mirna_K[mirna][0:2]
    mirna_K_file.close()

    # write the mirna divergence values into the output file
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'has_motif' + '\t' + 'K_hairpin' + '\t' + 'K_mature' + '\t' + 'K_backbone' + '\n')
    for mirna in mirna_with_motif:
        if mirna in mirna_K and len(mirna_K[mirna]) !=0:
            myfile.write(mirna + '\t' + 'motif' + '\t')
            for divergence in mirna_K[mirna][:-1]:
                myfile.write(divergence + '\t')
            myfile.write(mirna_K[mirna][-1] + '\n')
    for mirna in mirna_no_motif:
        if mirna in mirna_K and len(mirna_K[mirna]) !=0:
            myfile.write(mirna + '\t' + 'no_motif' + '\t')
            for divergence in mirna_K[mirna][:-1]:
                myfile.write(divergence + '\t')
            myfile.write(mirna_K[mirna][-1] + '\n')
    myfile.close()

    
       
        
def partition_mirna_with_without_all_motif(motif1_K_file, motif2_K_file, motif3_K_file, motif4_K_file, motif5_K_file, mirna_divergence, mirna_upstream, outputfile):
    '''
    (file, file, file, file, file, file, file, file) -> file
    Parition mirna divergence into 2 categories for mirnas in mirna_upstream that have or not a motif by combining all motifs in the motif_K_files 
    '''

    # make a set of mirnas that have a motif in motif_K_file
    mirna_with_motif = set()

    motif_filenames = [motif1_K_file, motif2_K_file, motif3_K_file, motif4_K_file, motif5_K_file]
    for motif_file in motif_filenames:
        motif_file = open(motif_file, 'r')
        line = motif_file.readline()
        for line in motif_file:
            line = line.rstrip()
            if line != '':
                line = line.split()
                mirna = line[0][:line[0].index('_')]
                mirna_with_motif.add(mirna)
        motif_file.close()

    # make a set of mirna that have no motif
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(mirna_upstream)
    mirna_no_motif = set()
    for mirna in upstreams:
        mirna_name = mirna[:mirna.index('_')]
        if not mirna_name in mirna_with_motif:
            mirna_no_motif.add(mirna_name)

    # make a dictionnary storing mirna divergence {mirna:[K_hairpin, K_mature, K_backbone]
    mirna_K_file = open(mirna_divergence, 'r')
    mirna_K = {}
    line = mirna_K_file.readline()
    for line in mirna_K_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            mirna_K[mirna] = line
            del mirna_K[mirna][0:2]
    mirna_K_file.close()

    # write the mirna divergence values into the output file
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'has_motif' + '\t' + 'K_hairpin' + '\t' + 'K_mature' + '\t' + 'K_backbone' + '\n')
    for mirna in mirna_with_motif:
        if mirna in mirna_K and len(mirna_K[mirna]) !=0:
            myfile.write(mirna + '\t' + 'motif' + '\t')
            for divergence in mirna_K[mirna][:-1]:
                myfile.write(divergence + '\t')
            myfile.write(mirna_K[mirna][-1] + '\n')
    for mirna in mirna_no_motif:
        if mirna in mirna_K and len(mirna_K[mirna]) !=0:
            myfile.write(mirna + '\t' + 'no_motif' + '\t')
            for divergence in mirna_K[mirna][:-1]:
                myfile.write(divergence + '\t')
            myfile.write(mirna_K[mirna][-1] + '\n')
    myfile.close()

        



def divergence_from_consensus(mast_output, mirna_upstream, consensus):
    '''
    (file, file, str) -> dict
    Compute divergence between each motif retrieved from mirna_upstream using the motif coordinates stored in mast_output and the motif consensus sequence
    '''
    from get_positions_from_mast import motif_position_from_mast
    from mirna_data_for_meme import convert_fasta
    from K_JC import compute_K_JC

    # get the motif positions
    
    mast_pos = motif_position_from_mast(mast_output)
    motif_positions = mast_pos[1]

    # convert the mirna upstream sequences into a dictionnary
    
    upstreams = convert_fasta(mirna_upstream)
    
    # make a dictionnary to store the motif sequences
    motif_seqs = {}
    for mirna in motif_positions:
        for position in motif_positions[mirna]:
            motif = upstreams[mirna][position: position + len(consensus)]
            if mirna in motif_seqs:
                motif_seqs[mirna].append(motif)
            else:
                motif_seqs[mirna] = [motif]
    
    # make a dictionary to store divergence for each motif in mirna upstream 
    motif_K = {}
    for mirna in motif_seqs:
        for motif in motif_seqs[mirna]:
            distance = compute_K_JC(consensus, motif)
            if mirna in motif_K:
                motif_K[mirna].append(distance)
            else:
                motif_K[mirna] = [distance]

    return motif_K

def match_distance_from_consensus_expression(mast_output, mirna_upstream, consensus, mirna_expression, outputfile):
    '''
    (file, file, str, file, file) -> file
    match mirna expression with divergence from the consensus sequence for each motif upstream the mirna
    '''

    # compute distance between consensus and each motif
    motif_K = divergence_from_consensus(mast_output, mirna_upstream, consensus)


    
    # make a dictionnary storing mirna expression {mirna:expression}
    expression_file = open(mirna_expression, 'r')
    expression = {}
    line = expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            expression[mirna] = line[1]
    expression_file.close()

     
    # write the mirna expression values and motif divergence for each motif and each mirna
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'motif divergence' + '\t' + 'expression' + '\n')
     
    for mirna in expression:
        for mirna_name in motif_K:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                for distance in motif_K[mirna_name]:
                    myfile.write(mirna_name + '\t' + str(distance) + '\t' + expression[mirna] + '\n')
               
    myfile.close()

   
def extract_motif_compute_K(mast_output_cbr, mast_output_csp9, Length, shift_cbr, shift_csp9, upstream_cbr, upstream_csp9):
    '''
    (file, file, int, int, int file, file) -> dict
    Extract the motifs of size Length in the upstream sequences of species 1 and 2 using their coordinates in the mast_output files,
    compute divergence between each motif and return a dictionnary storing the divergence values
    Eventually shift the motif so as to extract the sequence in common between the motifs of the 2 species
    '''

    from get_positions_from_mast import motif_position_from_mast
    from mirna_data_for_meme import convert_fasta
    from K_JC import compute_K_JC
    from get_positions_from_mast import get_mast_motif_shifted_right


    # get the mostif sequences into a dictionnary {mirna: [motif1, motif1]}
    motif_seq_cbr = get_mast_motif_shifted_right(mast_output_cbr, Length, shift_cbr, upstream_cbr)
    motif_seq_csp9 = get_mast_motif_shifted_right(mast_output_csp9, Length, shift_csp9, upstream_csp9)

    # create a dictionnary for sp1 and sp2 with each species motifs but with the same mirna names
    motifs_cbr = {}
    for mirna in motif_seq_cbr:
        mirnaname = mirna[mirna.index('-') + 1: mirna.index('_')]
        motifs_cbr[mirnaname] = motif_seq_cbr[mirna]

    motifs_csp9 = {}
    for mirna in motif_seq_csp9:
        mirnaname = mirna[mirna.index('_')+1:]
        motifs_csp9[mirnaname] = motif_seq_csp9[mirna]

    # make sets of mirnas present in one species but not the other
    mirna_remove_from_cbr = set()
    for mirna in motifs_cbr:
        if mirna not in motifs_csp9:
            mirna_remove_from_cbr.add(mirna)

    mirna_remove_from_csp9 = set()
    for mirna in motifs_csp9:
        if mirna not in motifs_cbr:
            mirna_remove_from_csp9.add(mirna)

    # remove the mirnas that are not shared between the 2 species
    for mirna in mirna_remove_from_cbr:
        del motifs_cbr[mirna]
    for mirna in mirna_remove_from_csp9:
        del motifs_csp9[mirna]

    # create a dictionnary to store the divergence values
    motif_K = {}
    for mirna in motifs_cbr:
        motif_K[mirna] = []

    # compute divergence between each motif in the order as they appear in the upstream sequences
    # if there are more motifs in one species ignore them
    for mirna in motifs_cbr:
        if len(motifs_cbr[mirna]) <= len(motifs_csp9[mirna]):
            for i in range(len(motifs_cbr[mirna])):
                gene1 = motifs_cbr[mirna][i]
                gene2 = motifs_csp9[mirna][i]
                divergence = compute_K_JC(gene1, gene2)
                if divergence != None:
                    motif_K[mirna].append(divergence)
        else:
            for i in range(len(motifs_csp9[mirna])):
                gene1 = motifs_csp9[mirna][i]
                gene2 = motifs_cbr[mirna][i]
                divergence = compute_K_JC(gene1, gene2)
                if divergence != None:
                    motif_K[mirna].append(divergence)
                    

    
    return motif_K
