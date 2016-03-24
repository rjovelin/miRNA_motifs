def cumulative_mirna_distribution(mast_output, inter_mirna_seqs, threshold, step, cumu_file):
    '''
    (file, file, int, int, file) -> None
    Build the cumulative frequency of miRNAs with at least one motif as a function of the inter_mirna length
    Starts the x-axis at threshold and increase by step. Store the values into the cumu_file
    '''

    from get_positions_from_mast import motif_position_from_mast
    from mirna_data_for_meme import convert_fasta

    # get the number of motif occurences
    motif_positions = motif_position_from_mast(mast_output)
    motifs = motif_positions[1]

    # store the inter_mirna sequences into a dictionnary
    inter_mirna = convert_fasta(inter_mirna_seqs)

    # make a dictionnary for each inter_mirna and store the number of motifs and the length of the sequence
    occurence = {}
    for mirna in inter_mirna:
        if mirna not in motifs:
            occurence[mirna] = [0, len(inter_mirna[mirna])]
        else:
            occurence[mirna] = [len(motifs[mirna]), len(inter_mirna[mirna])]

    # find the longest inter_mirna sequence
    lengths = [occurence[mirna][1] for mirna in occurence]
    longest = max(lengths)

    # open file to store the inter_mirna lengths and frequency of mirna with motif
    myfile = open(cumu_file, 'w')
    myfile.write('miRNA frequency' + '\t' + 'inter_mirna distance' + '\n')
            
    # count the frequency of mirnas with at least one motif starting at threshold
        
    while threshold < longest + step:
        total = 0
        for mirna in occurence:
            if occurence[mirna][1] <= threshold and occurence[mirna][0] != 0:
                total += 1
        myfile.write(str(total/len(occurence)) +'\t' + str(threshold) + '\n')
        threshold += step


    myfile.close()
    
