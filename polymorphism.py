def count_polymorphic_sites(fasta):
    '''
    (file) -> int
    Count the number of polymorphic sites in an multiple sequence alignment fasta file
    Precondition: all sequences are aligned and have the same length
    '''

    # convert the alignment file into a dictionnary
    from mirna_data_for_meme import convert_fasta
    alignment = convert_fasta(fasta)
    # fix the case of all sequences and make sure gaps are '-'
    for gene in alignment:
        alignment[gene] = alignment[gene].upper()
        alignment[gene] = alignment[gene].replace('~', '-')
       

    # make a list of keys
    genes = []
    for gene in alignment:
        genes.append(gene)

    polymorphisms = 0
    reference = genes[0]
    for i in range(len(alignment[reference])):
        nucleotides = set()
        for gene in alignment:
            nucleotides.add(alignment[gene][i])
        # exclude gap sites, Ns and sites, ambiguous sites and sites with more than 2 alleles
        for base in nucleotides:
            if base not in {'A', 'T', 'C', 'G'}:
                nucleotides = set()
        if len(nucleotides) == 2:
            polymorphisms += 1

    return polymorphisms
            
    
def compute_theta(fasta):
    '''
    (file) -> float
    Return the Waterson estimator theta per site of nuleotide diversity for the fasta alignment file
    '''

    polymorphisms = count_polymorphic_sites(fasta)

    # compute a: sum of 1/i from 1 to N-1 with N being the sample size
    from mirna_data_for_meme import convert_fasta
    alignment = convert_fasta(fasta)
    N = len(alignment)

    a = 0
    for i in range(1, N):
        a += 1/i

    # compute theta per sequence
    theta_sequence = polymorphisms / a

    # make a list of keys
    genes = []
    for gene in alignment:
        genes.append(gene)
    reference = genes[0]
    Total_Length = len(alignment[reference])
    
    # compute the number of sites excluding gaps, sites with Ns or ambiguous sites
    N_sites = Total_Length
    for i in range(Total_Length):
        nucleotides = set()
        excluded_sites = set()
        for gene in alignment:
            nucleotides.add(alignment[gene][i])
        for base in nucleotides:
            if base not in {'A', 'T', 'C', 'G'}:
                excluded_sites.add(base)
        if len(excluded_sites) != 0:
            N_sites -= 1

    # compute theta per site
    theta = theta_sequence / N_sites

    # round theta to the 6th digit
    return round(theta, 6)


def count_polymorphic_sites_motif(alignment, start, stop):
    '''
    (dict, str, int, int) -> int
    Count the number of polymorphic sites in the alignment between positions start and stop
    Precondition: all sequences are aligned and have the same length
    '''

    polymorphisms = 0
    for i in range(start, stop):
        nucleotides = set()
        for gene in alignment:
            nucleotides.add(alignment[gene][i])
        # exclude gap sites, Ns and sites, ambiguous sites and sites with more than 2 alleles
        for base in nucleotides:
            if base not in {'A', 'T', 'C', 'G'}:
                nucleotides = set()
        if len(nucleotides) == 2:
            polymorphisms += 1

    return polymorphisms

def compute_theta_motif(fasta, mast_output, Length, upstreams, outputfile):
    '''
    (file, file, int, file, file) -> file
    Return a dictionnary with theta for a given motif of size Length found in a mirna upstreams
    with positions stored in mast_output from an multiple sequence alignment fasta file
    Precondition: all sequences are aligned and have the same length
    '''

    from get_positions_from_mast import motif_position_from_mast
    from mirna_data_for_meme import convert_fasta
    import re

    # convert the alignment file into a dictionnary
    alignment = convert_fasta(fasta)
    # fix the case of all sequences and make sure gaps are '-'
    for gene in alignment:
        alignment[gene] = alignment[gene].upper()
        alignment[gene] = alignment[gene].replace('~', '-')

    # find the reference sequence containing the motif positions in alignment
    for gene in alignment:
        if 'crm' in gene:
            reference = gene
        
    # get the positions from mast {mirna: [positions]}
    position_tuple = motif_position_from_mast(mast_output)
    motif_positions = position_tuple[1]

    # create a dictionnary to store the motifs {mirna:[]}
    motifs_seq = {}
    for mirna in motif_positions:
        motifs_seq[mirna] = []

    # get the motifs from the file of upstream sequense
    # store the motifs in the dictionary {mirna:[motifs]}
    upstream_seqs = convert_fasta(upstreams)
    for mirna in motif_positions:
        for position in motif_positions[mirna]:
            motif = upstream_seqs[mirna][position: position + Length].upper()
            motifs_seq[mirna].append(motif)
    
    # create a dictionnary with reference as key and set of patterns as value
    # store the patterns in a set so as to avoid duplicate patterns
    # each pattern must be unique because each pattern will find all motifs using regex
    # if duplicate patterns are kept, divergence will be computed multiple times on the same motif
    mirna_pattern = {}

    if reference in motifs_seq:
        for motif in motifs_seq[reference]:
            pattern_list = []
            pattern = ''
            # a pattern include 3 gaps at most at each position
            for i in range(len(motif)-1):
                motif_with_gap = motif[0:i+1] +'-{0,3}' + motif[i+1:]
                pattern_list.append(motif_with_gap)
            for new_motif in pattern_list:
                pattern = pattern + new_motif + '|'
            pattern = pattern[0:-1]
            if reference in mirna_pattern:
                mirna_pattern[reference].add(pattern)
            else:
                mirna_pattern[reference] = {pattern,}
   
    # create a dictionnary to store the theta values of each motif for a given mirna
    mirna_polym = {}
    if reference in mirna_pattern:
        for pattern in mirna_pattern[reference]:
            for motif in re.finditer(pattern, alignment[reference]):
                start = motif.start()
                stop = motif.start() + len(motif.group())
                N_polym = count_polymorphic_sites_motif(alignment, start, stop)
                N_sites = len(motif.group())
                a = 0
                N = len(alignment)
                for i in range(1, N):
                    a += 1/i
                theta_sequence = N_polym / a
                for i in range(N_sites):
                    nucleotides = set()
                    excluded_sites = set()
                    for gene in alignment:
                        nucleotides.add(alignment[gene][i])
                    for base in nucleotides:
                        if base not in {'A', 'T', 'C', 'G'}:
                            excluded_sites.add(base)
                    if len(excluded_sites) != 0:
                        N_sites -= 1
                if N_sites != 0:
                    theta = theta_sequence / N_sites
                if reference in mirna_polym:
                    mirna_polym[reference].append(theta)
                else:
                    mirna_polym[reference] = [theta]
    
                
        # store the theta values into a file
        if mirna_polym != {}:
            with open(outputfile, 'a+') as myfile:
                for theta in mirna_polym[reference]:
                    myfile.write(reference + '\t' + str(theta) + '\n')
    

def randomization_theta_motif(Length, N_repeats, outputfile):
    '''
    (int, int, file) -> file
    Compute theta for N_repeats number of motifs of size Length chosen at random
    and save the values in the outputfile
    '''

    import random
    from mirna_data_for_meme import convert_fasta
    
    
    filenames = {0: 'mir-79.txt', 1: 'mir-80.txt', 2: 'mir-81.txt', 3: 'mir-82.txt', 4: 'mir-83.txt',
                 5: 'mir-85.txt', 6: 'mir-45.txt', 7: 'mir-46.txt', 8: 'mir-47.txt', 9: 'mir-48.txt',
                 10: 'mir-49.txt', 11: 'mir-50.txt', 12: 'mir-52-2.txt', 13: 'mir-57.txt', 14: 'mir-58b.txt',
                 15: 'mir-59.txt', 16: 'mir-60.txt', 17: 'mir-61.txt', 18: 'mir-62.txt', 19: 'mir-67.txt',
                 20: 'mir-70.txt', 21: 'mir-71.txt', 22: 'mir-72.txt', 23: 'mir-73.txt', 24: 'mir-75.txt',
                 25: 'mir-76.txt', 26: 'mir-77-2.txt', 27: 'mir-252.txt', 28: 'mir-253.txt', 29: 'mir-259.txt',
                 30: 'mir-356.txt', 31: 'mir-86.txt', 32: 'mir-87.txt', 33: 'mir-90.txt', 34: 'mir-124.txt',
                 35: 'mir-124b.txt', 36: 'mir-231.txt', 37: 'mir-232-1.txt', 38: 'mir-232-2.txt', 39: 'mir-233.txt',
                 40: 'mir-235.txt', 41: 'mir-236.txt', 42: 'mir-237.txt', 43: 'mir-238.txt', 44: 'mir-239a.txt',
                 45: 'mir-239b.txt', 46: 'mir-240.txt', 47: 'mir-241.txt', 48: 'mir-244.txt', 49: 'mir-245.txt',
                 50: 'mir-246-1.txt', 51: 'mir-246-2.txt', 52: 'mir-247.txt', 53: 'mir-248.txt', 54: 'mir-249.txt',
                 55: 'mir-2230.txt', 56: 'mir-2231.txt', 57: 'mir-2254.txt', 58: 'mir-7601.txt', 59: 'mir-7604.txt',
                 60: 'mir-7606.txt', 61: 'mir-7608.txt', 62: 'mir-784.txt', 63: 'mir-787.txt', 64: 'mir-788.txt',
                 65: 'mir-1817b.txt', 66: 'mir-2228.txt', 67: 'mir-35i.txt', 68: 'mir-38.txt', 69: 'mir-42.txt',
                 70: 'let-7.txt', 71: 'lin-4.txt', 72: 'mir-1.txt', 73: 'mir-2.txt', 74: 'mir-34.txt', 75: 'mir-84.txt'}

    thetas = []

    while N_repeats > 0:
        # draw a random file from all upstream files
        file_number = random.randint(0, len(filenames)-1)
        
        # make a dictionnary with sequences in the file
        alignment = convert_fasta(filenames[file_number])

        # fix the case of all sequences and make sure gaps are '-'
        for gene in alignment:
            alignment[gene] = alignment[gene].upper()
            alignment[gene] = alignment[gene].replace('~', '-')
        
        for mirna in alignment:
            if 'crm' in mirna:
                reference = mirna
             
        # pick a random position in the sequence and grab a motif of size Length        
        start = random.randint(0, len(alignment[reference])- Length)
        stop = start + Length


        N_polym = count_polymorphic_sites_motif(alignment, start, stop)
        N_sites = len(alignment[reference][start:stop])
        a = 0
        N = len(alignment)
        for i in range(1, N):
            a += 1/i
        theta_sequence = N_polym / a


        for i in range(N_sites):
            nucleotides = set()
            excluded_sites = set()
            for gene in alignment:
                nucleotides.add(alignment[gene][i])
            for base in nucleotides:
                if base not in {'A', 'T', 'C', 'G'}:
                    excluded_sites.add(base)
            if len(excluded_sites) != 0:
                N_sites -= 1
        if N_sites != 0:
            theta = theta_sequence / N_sites
            thetas.append(theta)
            N_repeats -= 1

    myfile = open(outputfile, 'w')
    for theta in thetas:
        myfile.write('random' + '\t' + str(theta) + '\n')
    myfile.close()
        

