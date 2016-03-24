# get the inter-mirna sequences           
def get_inter_mirna_seq_from_clusters(mirna_coord_file, genome_file):   
    '''
    (file, file) -> dict
    Use the mirna coordinates stored in the mirna_coord file and return a dictionnary of strings corresponding to
    the sequences between the clustered mirnas using the genome sequence stored into the genome_file
    '''
    
    from reverse_complement import rev_compl
    from mirna_data_for_meme import convert_fasta
    from mirna_data_for_meme import mirna_indices_from_file
   
    # convert genome_file into a genome dictionnary {chromo: sequence}
    genome = convert_fasta(genome_file)

    # get the mirna coordinates from mirna_coord_file
    # mirna_coord is a dict with {'mirna':['chromosome', 'start', 'end', 'strand']
    mirna_coord = mirna_indices_from_file(mirna_coord_file)
    
    # make a dictionnary with chromosome as key
    # the mirna_coord dictionnary is unchanged
    # chromosome is a dictionnary in the form {chromo: [start, (end, mirna, strand)]
    chromosome = {}
    for mirna in mirna_coord:
        chromo = mirna_coord[mirna][0]
        if chromo in chromosome:
            chromosome[chromo].append([int(mirna_coord[mirna][1]), (int(mirna_coord[mirna][2]), mirna, mirna_coord[mirna][3])])
        else:
            chromosome[chromo] = [[int(mirna_coord[mirna][1]), (int(mirna_coord[mirna][2]), mirna, mirna_coord[mirna][3])]]


    # create a dictionnary to store the sequences
    inter = {}

    # go through each chromosome and store clustered mirnas into a list, except the first mirna of the cluster
    # sort the chromosome values according to the start value
    for chromo in chromosome:
        chromosome[chromo].sort()
        for i in range(len(chromosome[chromo])-1):
            d = int(chromosome[chromo][i][1][0]) - int(chromosome[chromo][i+1][0])
            # add strand condition keep first one if positif strand, keep last one if minus strand
            if abs(d) < 1000 and  chromosome[chromo][i][1][-1] == chromosome[chromo][i+1][1][-1]:
                start = int(chromosome[chromo][i][1][0])
                end = int(chromosome[chromo][i+1][0])
                sequence = genome[chromo][start:end]
                mirna = chromosome[chromo][i][1][1]
                # use the name of mirna upstream the inter-mirna sequence as key when sense = +
                # use the name of the mirna downstream the inter-mirna sequence as key when sense = - and take the reverse complement
                if chromosome[chromo][i][1][-1] == '+':
                    inter[mirna] = sequence
                if chromosome[chromo][i][1][-1] == '-':
                    inter[mirna] = rev_compl(sequence)

    return inter

def store_inter_mirna_to_file(inter, outputfile, threshold):
    '''
    (dict, file)
    Store the inter_mirna sequences from the inter dictionnary in the outputfile with a fasta format
    if the length of the sequence is greater or equal than the mininmum threshold
    '''

    myfile = open(outputfile, 'w')
    for mirna in inter:
        if len(inter[mirna]) >= threshold:
            myfile.write('>' + mirna + '\n')
            myfile.write(inter[mirna] + '\n')
    myfile.close()

    
