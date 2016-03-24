def motif_position_from_mast(mast_output):
    '''
    (str) -> (dict, dict)
    Read the file output from MAST and extract the positions of the motif
    and return 2 dictionnaries of {mirna: [positions]},
    the first scaled to 1000, and the second with raw positions 
    '''


    # example of file output
    # mirna_name   position    scaled_position

    mastfile = open(mast_output, 'r')


    line = mastfile.readline()
    
    # read unil section II is reached
    while 'SECTION II' not in line:
        line = mastfile.readline()
    # skip over the first line of ****
    line = mastfile.readline()

    # keep readling until the line under the header is reached
    while '--' not in line:
        line = mastfile.readline()

    # store the positions and the scaled positions in 2 dictionnaries that have the same keys
    # if the length of the upstream is < 1000: scale the positions to 1000 to plot the histogramm
    # mirna_pos = {mirna: [positions]}
    # mirna_scaled = {mirna: [scaled positions]}
    mirna_pos = {}
    mirna_scaled = {}

    # extract position for each motif
    # stop reading when the line with *** is reached
    for line in mastfile:
        if '*' in line:
            break     
        else:
            # this is parsing the occurence information stored in the form d_[o]_d
            line = line.rstrip()
            if line!= '':
                line = line.split()
                line.remove(line[1])
                line[1] = line[1].split('_')
                # add the positions to the dictionnary of position
                count = line[1].count('[1]')
                i = 0
                mirna_pos[line[0]] = []
                position = 0
                while count > 0:
                    if i == 0:
                        position = int(line[1][i])
                    else:
                        position = position + int(line[1][i]) + 10
                    mirna_pos[line[0]].append(position)
                    i += 2
                    count -= 1

                # add the scaled positions to the dictionnary of scaled positions
                # scale the position to 1000 is the length of the upstream sequence is < 1000
                n_motif = line[1].count('[1]')
                total = line[1].count('[1]')
                j = 0
                length = 0
                mirna_scaled[line[0]] = []
                while n_motif > 0:
                    length += int(line[1][j])
                    j += 2
                    n_motif -= 1
                length = length + int(line[1][-1]) + total * 10
                for position in mirna_pos[line[0]]:
                    if length < 1000:
                        mirna_scaled[line[0]].append((position * 1000) / length)
                    else:
                        mirna_scaled[line[0]].append(position)
                                    
                    
    return mirna_scaled, mirna_pos    
        
    
def mast_motif_position_printer(mirna_pos, output_file):
    '''
    (dict) -> None
    Save the dictionnary mirna_scaled of positions for each mirna resulting from the MAST match for a given motif in the output_file
    '''

    positionfile = open(output_file, 'w')
    positionfile.write('mirna' + '\t' + 'position' + '\n')
    for key in mirna_scaled:
        positionfile.write(key + '\t')
        for item in mirna_scaled[key][:-1]:
            positionfile.write(str(item) + '\t')
        positionfile.write(str(mirna_scaled[key][-1]) + '\n')
        
    positionfile.close()



    
def plot_histogram_motif_positions(mirna_scaled, window, histo_file):
    '''
    (dict, str) -> None
    Build the histogram of motif position located within the window and save it into the histo_file
    '''

    # make a list of size window that contains only 0s:
    # each value in the list is the count of position for the range [0 - window[ etc
    range_counts = [0] * (1000 // window)

    # determine the index in the list range_count where the position should be added
    # count the number of times a position appear within the window
    for mirna in mirna_scaled:
        for position in mirna_scaled[mirna]:
            which_range = int(position) // window
            range_counts[which_range] += 1

    histogram = open(histo_file, 'w')
    histogram.write('Range' + '\t' + 'Count' + '\n')

    for i in range(len(range_counts)):
        histogram.write('[' + str(i * window) + '-' + str((i * window) + window -1) + ']' + '\t')
        histogram.write(str(range_counts[i]) + '\n')

    histogram.close()
        
        
        
def get_number_mast_motif_combined(mast_output1, mast_output2, storage_file, species, upstream_seq):
    '''
    (str, str, str, int. str) -> File
    Save the number of occurences of 2 MAST motifs in species from the mast_output1 and mast_output2 into the storage_file.
    Get the number of mirnas with 0 occurence by comparing mast_output and upstream_seqs
    '''

    outputfile = open(storage_file, 'a+')

    # get the motif positions from the mast_output file
    mast_positions1 = motif_position_from_mast(mast_output1)
    motif_positions1 = mast_positions1[1]

    mast_positions2 = motif_position_from_mast(mast_output2)
    motif_positions2 = mast_positions2[1]


    # make a dictionnary storing the sequence length
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(upstream_seq)

    seq_length = {}
    for mirna in upstreams:
        upstreams[mirna] = upstreams[mirna].upper()
        N_number = upstreams[mirna].count('N')
        total_length = len(upstreams[mirna]) - N_number
        seq_length[mirna] = total_length

    # make a set of mirna names for mirnas that have a motif
    mirna_with_motif = set()

    # normalize the motif count by the length of the upstream sequence
    for mirna in motif_positions1:
        if len(motif_positions1[mirna]) != 0:
            outputfile.write(species + '\t' + str(len(motif_positions1[mirna])/seq_length[mirna]) + '\n')
            mirna_with_motif.add(mirna)

    for mirna in motif_positions2:
        if len(motif_positions2[mirna]) != 0:
            outputfile.write(species + '\t' + str(len(motif_positions2[mirna])/seq_length[mirna]) + '\n')
            mirna_with_motif.add(mirna)

    for mirna in upstreams:
        if mirna not in mirna_with_motif:
            outputfile.write(species + '\t' + '0' + '\n')
       
    outputfile.close()




def get_number_mast_motif(mast_output, storage_file, species, upstream_seq):
    '''
    (str, str, str, int. str) -> File
    Save the number of occurences of a given MAST motif in species from the mast_output into the storage_file.
    Get the number of mirnas with 0 occurence by comparing mast_output and upstream_seqs
    '''

    outputfile = open(storage_file, 'a+')

    # get the motif positions from the mast_output file
    mast_positions = motif_position_from_mast(mast_output)
    motif_positions = mast_positions[1]

    # make a dictionnary storing the sequence length
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(upstream_seq)

    seq_length = {}
    for mirna in upstreams:
        upstreams[mirna] = upstreams[mirna].upper()
        N_number = upstreams[mirna].count('N')
        total_length = len(upstreams[mirna]) - N_number
        seq_length[mirna] = total_length

    # make a set of mirna names for mirnas that have a motif
    mirna_with_motif = set()

    # normalize the motif count by the length of the upstream sequence
    for mirna in motif_positions:
        if len(motif_positions[mirna]) != 0:
            outputfile.write(species + '\t' + str(len(motif_positions[mirna])/seq_length[mirna]) + '\n')
            mirna_with_motif.add(mirna)

    for mirna in upstreams:
        if mirna not in mirna_with_motif:
            outputfile.write(species + '\t' + '0' + '\n')
       
    outputfile.close()



def count_mirna_with_motif(mast_output, upstream_seq):
    '''
    (file, file) -> None
    Print the number of mirnas that have at least 1 motif in the mast_output file, the total of number of occurence of that motif,
    and the % of mirnas that a motif relative to the total number of mirnas in the upstream_seq file
    '''

    mast_positions = motif_position_from_mast(mast_output)
    motif_positions = mast_positions[1]

    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(upstream_seq)

    mirna_with_motif = len(motif_positions)
    total_mirna = len(upstreams)
    percent_mirna = mirna_with_motif/total_mirna * 100

    total_sites = 0

    for mirna in motif_positions:
        total_sites += len(motif_positions[mirna])
        
    print(mirna_with_motif)
    print(percent_mirna)
    print(total_sites)
 
    

def get_mast_motif(mast_output, Length, upstreams):
    '''
    (file, int, file) -> dict
    Extract the motifs of size Length from the upstream sequences
    using the coordinates stored in the mast_output file
    '''

    from mirna_data_for_meme import convert_fasta
    
    # get the positions from mast {mirna: [positions]}
    mast_positions = motif_position_from_mast(mast_output)
    motif_positions = mast_positions[1]

    # create a dictionnary to store the motifs {mirna:[]}
    motifs_seq = {}
    for mirna in motif_positions:
        motifs_seq[mirna] = []

    # get the motifs from the file of upstream sequense
    # store the motifs in the dictionary {mirna:[motifs]}
    upstream_seqs = convert_fasta(upstreams)
    for mirna in motif_positions:
        for position in motif_positions[mirna]:
            motif = upstream_seqs[mirna][position: position + Length]
            motifs_seq[mirna].append(motif)

    return motifs_seq


    
def get_mast_motif_shifted_right(mast_output, Length, shift_right, upstreams):
    '''
    (file, int, file) -> dict
    Extract the motifs of size Length from the upstream sequences
    using the coordinates stored in the mast_output file
    Shift the motif to the right by adding shift_right to the start index
    '''

    from mirna_data_for_meme import convert_fasta
    
    # get the positions from mast {mirna: [positions]}
    mast_positions = motif_position_from_mast(mast_output)
    motif_positions = mast_positions[1]

    # create a dictionnary to store the motifs {mirna:[]}
    motifs_seq = {}
    for mirna in motif_positions:
        motifs_seq[mirna] = []

    # get the motifs from the file of upstream sequense
    # store the motifs in the dictionary {mirna:[motifs]}
    upstream_seqs = convert_fasta(upstreams)
    for mirna in motif_positions:
        for position in motif_positions[mirna]:
            motif = upstream_seqs[mirna][position + shift_right: position + shift_right + Length]
            motifs_seq[mirna].append(motif)

    return motifs_seq
