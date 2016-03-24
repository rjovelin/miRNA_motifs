# convert nematode genomes into a dictionnary
def convert_fasta(fasta):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    sequence ID as key and single string sequence as value
    '''
    
    genome = {}
    with open(fasta, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line
    return genome



# store nematode genome with chromosome/contig as sequence header
# and single string sequence into a file for future access
def store_genome(genome, storage):
    '''
    (dict) -> None
    Store the keys, values of dictionnary genome
    into a file with the key as sequence header
    and the value as single string fasta sequence
    '''

    newfile = open(storage, 'w')
        
    for chromosome in genome:
        newfile.write('>' + chromosome + '\n')
        newfile.write(genome[chromosome] + '\n')

    newfile.close()


# check which mirnas sequence are identical
# eventually adjust the mirna names in the table of mirna coordinates to match mirBase mirna names with identical sequences
def check_identical_mirnas(mirna):
    '''
    (str) -> list
    Lookup for identical hairpin sequences in the list of mirnas for a given species
    '''

    species = open(mirna, 'r')

    all_mirna = {}
    
    for line in species:
        line = line.strip()
        if line == '':
            continue
        elif line.startswith('>'):
            all_mirna[line[1:]] = ''
            mir = line[1:]
        else:
            all_mirna[mir] += line

    mir_seq = {}
    duplicates = []

    for mir, seq in all_mirna.items():
        if seq in mir_seq:
            mir_seq[seq].append(mir)
        else:
            mir_seq[seq] = [mir]

    for seq in mir_seq:
        if len(mir_seq[seq]) > 1:
            duplicates.append(mir_seq[seq])

    species.close()

    return duplicates

def remove_duplicates(duplicates, mirseq):
    '''
    (list, dict) -> dict
    Remove mirnas from the dictionnary of mirna sequences mirseq
    when they are identical to other mirnas stored in the list of duplicates
    until 1 unique sequence per mirna is kept
    '''

    for i in range(len(duplicates)):
        for j in range(1, len(duplicates[i])):
            del mirseq[duplicates[i][j]]

    return mirseq
    
# get the coordinates of each mirna for a given species
def find_mirna_indices(genome, mirna):
    '''
    (dict, str) -> dict
    Find the coordinates (start, end) of each gene in the mirna file within
    the dictionnary genome and store them in a dictionnary with mirna:[chromo, start, end, sense]
    '''


    # make a dictionnary storing mirna:seq for each mirna
    mirseq = convert_fasta(mirna)


    # check ientical sequences in file mirna
    duplicates = check_identical_mirnas(mirna)

    # remove duplicates from mirseq
    mirseq = remove_duplicates(duplicates, mirseq)


    # check that each chromo sequence is lower, and if not then concert it to lower
    for chromo in genome:
        if not genome[chromo].islower():
            genome[chromo] = genome[chromo].lower()

    # convert the sequence in mirseq in lower
    for mir in mirseq:
        mirseq[mir] = mirseq[mir].lower()
            
  
    from reverse_complement import rev_compl

    # make a dictionnary for each chromosome and store start for checking that a mirna is not already recorded
    start_check = {}
    for chromo in genome:
        start_check[chromo] = []

    # make a dictionnary to store information about each mirna location mir:[chromo, start, end, sense]
    mir_indices = {}
    for mir in mirseq:
        mir_indices[mir] = []


    for mir in mirseq:
        for chromo in genome:
            occurence = genome[chromo].count(mirseq[mir])     # count the number of times seq appears in the chromosome
            if occurence != 0:
                more = True
                i = 0
                start = 0
                while occurence > 0 and more:
                    start = genome[chromo].find(mirseq[mir], start + i)
                    if start not in start_check[chromo]:         # if start is not already recorded, record the mirna in mir_indices and record start in start_check
                        start_check[chromo].append(start)
                        end = start + len(mirseq[mir])
                        mir_indices[mir].append([chromo, start, end, '+'])
                        i += 1
                        occurence -= 1
                    else:
                        more = False
                        


    for mir in mirseq:
        for chromo in genome:
            rev_mir = rev_compl(mirseq[mir])                  # take the receverse complement of mir
            occurence = genome[chromo].count(rev_mir)     # count the number of times seq appears in the chromosome
            if occurence != 0:
                more = True
                i = 0
                start = 0
                while occurence > 0 and more:
                    start = genome[chromo].find(rev_mir, start + i)
                    if start not in start_check[chromo]:           # if start is not already recorded, record the mirna in mir_indices and record start in start_check
                        start_check[chromo].append(start)
                        end = start + len(mirseq[mir])
                        mir_indices[mir].append([chromo, start, end, '-'])
                        i += 1
                        occurence -= 1
                    else:
                        more = False
                        

    return mir_indices

# save information about mirna positions into a file
def store_mirna_coordinates(mir_indices, coordinates):
    '''
    (dict, str) -> None
    Store the dictionnary of mirna coordinates into a file
    '''
    # print the gene entries into a new file, starting with the header
    mircoord = open(coordinates, 'w')

    # write header to the updated file
    header = ['mirna', 'chromosome', 'start', 'end', 'strand']
    for item in header[:-1]:
        mircoord.write(item + '\t')
    mircoord.write(header[-1] + '\n')


    # write entries for each mirna
    for mir in mir_indices:
        if len(mir_indices[mir]) == 0:
            mircoord.write(mir + '\t')
            mircoord.write('not found\n')
        else:
            for i in range(len(mir_indices[mir])):
                mircoord.write(mir + '\t')
                mircoord.write(mir_indices[mir][i][0] +'\t')
                mircoord.write(str(mir_indices[mir][i][1] + 1) + '\t')
                mircoord.write(str(mir_indices[mir][i][2]) + '\t')
                mircoord.write(mir_indices[mir][i][-1] + '\n')
                
    mircoord.close()


# save mirna indices and positional information in a file
def mirna_indices_to_file(mir_indices, coordinates):
    '''
    (dict, str) -> None
    Store the dictionnary of mirna coordinates into a file
    '''
    # print the gene entries into a new file, starting with the header
    mircoord = open(coordinates, 'w')

    # write header to the updated file
    header = ['mirna', 'chromosome', 'start', 'end', 'strand']
    for item in header[:-1]:
        mircoord.write(item + '\t')
    mircoord.write(header[-1] + '\n')


    # write entries for each mirna
    for mir in mir_indices:
        for i in range(len(mir_indices[mir])):
            mircoord.write(mir + '\t')
            for j in range(len(mir_indices[mir][i])-1):
                mircoord.write(str(mir_indices[mir][i][j]) + '\t')
            mircoord.write(mir_indices[mir][i][-1] + '\n')

    mircoord.close()

# get back mirna indices stored in file
def mirna_indices_from_file(mir_indices):
    '''
    (str) -> dict
    Make a dictionnary with mirna positional information stored in filename
    in the form of {'mir': [chromo, start, end, strand]
    '''

    mir_coord = {}
    indice_file = open(mir_indices, 'r')
    indice_file.readline()
    for line in indice_file:
        line = line.strip()
        if line != '':
            entries = line.split()
            mir = entries[0]
            coordinates = entries[1:]
            mir_coord[mir] = coordinates

    return mir_coord
       

# identify the mirnas located within clusters            
def is_mirna_clustered(mirna_coord):   
    '''
    (dict) -> list
    Take the dictionnary of mirna indices and return a list of mirnas located in clusters, but not the first mirna of the cluster
    '''

    
    # mirna_coord is a dict with {'mirna':['chromosome', 'start', 'end', 'strand']

    chromosome = {}
    sense_clusters = []
    anti_sense_clusters = []
    
    # make a dictionnary with chromosome as key
    # the mirna_coord dictionnary is unchanged
    for mirna in mirna_coord:
        chromo = mirna_coord[mirna][0]
        if chromo in chromosome:
            chromosome[chromo].append([int(mirna_coord[mirna][1]), (int(mirna_coord[mirna][2]), mirna, mirna_coord[mirna][3])])
        else:
            chromosome[chromo] = [[int(mirna_coord[mirna][1]), (int(mirna_coord[mirna][2]), mirna, mirna_coord[mirna][3])]]

    # go through each chromosome and store clustered mirnas into a list, except the first mirna of the cluster
    for chromo in chromosome:
        chromosome[chromo].sort()
        for i in range(len(chromosome[chromo])-1):
            d = int(chromosome[chromo][i][1][0]) - int(chromosome[chromo][i+1][0])
            # add strand condition keep first one if positif strand, keep last one if minus strand
            if abs(d) < 1000 and  chromosome[chromo][i][1][-1] == chromosome[chromo][i+1][1][-1]:
                if chromosome[chromo][i][1][-1] == '+':
                    sense_clusters.append(chromosome[chromo][i+1][1][1])
                if chromosome[chromo][i][1][-1] == '-':
                    anti_sense_clusters.append(chromosome[chromo][i][1][1])
    
    sense_clusters.extend(anti_sense_clusters)

    return sense_clusters



def mirna_for_meme(mirna_indices):
    '''
    (str) -> dict
    Open the file mir_indices and remove the mirnas located in clusters, but keep the first mirna of the cluster
    Return a dictionnary of coordinates for the set of mirnas for MEME analysis
    '''

    # get back the indices stored in file
    mir_coord = mirna_indices_from_file(mirna_indices)
    
    # find the clustered mirnas
    mirnas_to_delete = is_mirna_clustered(mir_coord)


    # remove the clustered mirnas from the dictionnary of coordinates
    for mirna in mirnas_to_delete:
        del mir_coord[mirna]

    return mir_coord

       

def grab_upstream(mir_indices, genome, Length):
    '''
    (str, dict, int) -> dict
    Return a dictionnary of upstream sequence of given Length for each mirna in the file of mir_indices from the dictionnary genome
    '''


    import reverse_complement
    
    # build the dataset of mirnas for MEME analysis
    mir_coord = mirna_for_meme(mir_indices)

    upstream_seqs = {}

    for mir in mir_coord:
        if mir_coord[mir][-1] == '+':
            chromo = mir_coord[mir][0]
            start = int(mir_coord[mir][1])
            if start - Length >= 0:
                upstream = genome[chromo][(start - Length):start]
                upstream_seqs[mir] = upstream
            else:
                upstream = genome[chromo][:start]
                upstream_seqs[mir] = upstream
        elif mir_coord[mir][-1] == '-':
            chromo = mir_coord[mir][0]
            end = int(mir_coord[mir][2])
            if end + Length <= len(genome[chromo]):
                upstream = genome[chromo][end: (end + Length)]
                upstream =reverse_complement.rev_compl(upstream)
                upstream_seqs[mir] = upstream
            else:
                upstream = genome[chromo][end:]
                upstream =reverse_complement.rev_compl(upstream)
                upstream_seqs[mir] = upstream
                
    return upstream_seqs
    

def upstream_to_file(upstream_seqs, outputfile):
    '''
    (dict, str) -> None
    Save fasta upstream sequence of dictionnary upstream_seqs into an outputfile for use in MEME analysis
    '''

    fasta_file = open(outputfile, 'w')

    for mir in upstream_seqs:
        fasta_file.write('>' + mir + '\n')
        fasta_file.write(upstream_seqs[mir] + '\n')

    fasta_file.close()
