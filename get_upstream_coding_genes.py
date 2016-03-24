def get_gene_coordinates_from_gff(species_genes):
    '''
    (file) -> dict
    Return a dictionary with the protein-coding genes as key
    amd a list of features as value
    Precondition: the gff file has been parsed to retain only information
    about protein coding genes on each line
    '''


    # make a dictionnary with the genes as key and list of feature as value
    # {gene: [chromosome, start, end, strand]}

    myfile = open(species_genes, 'r')

    gene_coord = {}

    for line in myfile:
        line = line.rstrip()
        line = line.split()
        gene_coord[line[8][8:line[8].index(';')]] = [line[0]]
        gene_coord[line[8][8:line[8].index(';')]].append(line[3])
        gene_coord[line[8][8:line[8].index(';')]].append(line[4])
        gene_coord[line[8][8:line[8].index(';')]].append(line[6])

    myfile.close()
    return gene_coord



def grab_upstream_random_coding_genes(species_genes, genome_file, gene_number, Length):
    '''
    (dict, file, int, int) -> dict
    Return a dictionnary of upstream sequence of given Length for each gene taken at random
    with upstream length  = 1000
    '''

    gene_coord = get_gene_coordinates_from_gff(species_genes)    

    import reverse_complement


    # convert the genome_file into a genome dictionnary
    # genome = {chromo: fasta}
    from mirna_data_for_meme import convert_fasta
    genome = convert_fasta(genome_file)


    # make a list of the keys in gene_coord
    genes = []
    for gene in gene_coord:
        genes.append(gene)

    # make a dictionnary with index as key and list of features as value
    # {i: [chromosome, start, end, strand, gene]}

    total_genes = {}
    for i in range(len(genes)):        
        total_genes[i] = gene_coord[genes[i]]
        total_genes[i].append(genes[i])



    # mke a dictionnary to store the upstream sequences
    upstream_gene_seqs = {}

    # select gene_number genes at random from the total_genes
    import random
    # make a set of genes already picked
    # because genes can be picked only once
    already_picked = set()
    while gene_number > 0:
        i = random.randint(0, len(total_genes)-1)
        if not i in already_picked:
            if total_genes[i][3] == '+':
                chromo = total_genes[i][0]
                start = int(total_genes[i][1])
                if start - Length >= 0:
                    upstream = genome[chromo][(start - Length):start]
                    count = upstream.lower().count('n')
                    if len(upstream) == Length and count == 0:
                        upstream_gene_seqs[total_genes[i][-1]] = upstream
                        already_picked.add(i)
                        gene_number -= 1
            if total_genes[i][3] == '-':
                chromo = total_genes[i][0]
                end = int(total_genes[i][2])
                if (end + Length) <= len(genome[chromo]):
                    upstream = genome[chromo][end: (end + Length)]
                    count = upstream.lower().count('n')
                    if len(upstream) == Length and count == 0:
                        upstream = reverse_complement.rev_compl(upstream)
                        upstream_gene_seqs[total_genes[i][-1]] = upstream
                        already_picked.add(i)
                        gene_number -= 1
            
         
    return upstream_gene_seqs





def count_genes_with_upstream_Length_no_N(genome_file, species_genes, Length):
    '''
    (file, file, int) -> None
    Count the number of genes in the genome that have upstream sequence of size Length using the gene coordinates
    stored in the species_genes file
    '''

    # convert the gene coordinates into a dictionary
    gene_coord = get_gene_coordinates_from_gff(species_genes)


    # convert the genome_file into a genome dictionnary
    from mirna_data_for_meme import convert_fasta
    genome = convert_fasta(genome_file)


    reverse = 0
    direct = 0

    for gene in gene_coord:
        if gene_coord[gene][3] == '-':
            chromo = gene_coord[gene][0]
            end = int(gene_coord[gene][2])
            if (end + Length) <= len(genome[chromo]):
                upstream = genome[chromo][end: (end + Length)]
                count = upstream.lower().count('n')
                if len(upstream) == Length and count == 0:
                    reverse += 1
        elif gene_coord[gene][3] == '+':
            chromo = gene_coord[gene][0]
            start = int(gene_coord[gene][1])
            if (start - Length) >= 0:
                upstream = genome[chromo][start - Length: start]
                count = upstream.lower().count('n')
                if len(upstream) == Length and count == 0:
                    direct += 1

    print('number of genes: ', direct + reverse)
    print('% of the total number of genes: ', (direct + reverse)/ len(gene_coord) * 100)


        
        
def random_genes_generator(gene_coord, gene_number):
    '''
    (dict, int) -> dict
    Pick gene_number of genes at random from the genome dictionnary
    '''

    # make a list of the keys in genome
    genes = []
    for gene in gene_coord:
        genes.append(gene)

    # make a dictionnary with index as key and list of features as value
    # {i: [chromosome, start, end, strand, gene]}

    total_genes = {}
    for i in range(len(genes)):        
        total_genes[i] = gene_coord[genes[i]]
        total_genes[i].append(genes[i])

    # select gene_number genes at random from the total_genes
    import random
    random_genes = {}
    # make a set of genes already picked
    # because genes can be picked only once
    already_picked = set()
    while gene_number > 0:
        i = random.randint(0, len(total_genes)-1)
        if i not in already_picked:
            random_genes[i] = total_genes[i]
            already_picked.add(i)
            gene_number -= 1

    return random_genes


    

def upstream_coding_to_file(upstream_gene_seqs, outputfile):
    '''
    (dict, str) -> file
    Save fasta upstream sequence of dictionnary upstream_gene_seqs into an outputfile for use in MAST analysis
    '''

    fasta_file = open(outputfile, 'w')

    for gene in upstream_gene_seqs:
        if len(upstream_gene_seqs[gene]) > 0:
            fasta_file.write('>' + gene + '\n')
            fasta_file.write(upstream_gene_seqs[gene] + '\n')

    fasta_file.close()



def is_mirna_intronic(species_cds_intron, mirna_coordinates):
    '''
    (file, file) -> dict
    Check if a mirna is located within an intron using the mirna and intron coordinates from files mirna_coordinates and species_cds_introns
    and return a dictionnary with this information for each mirna
    '''

    # make a dictionnary to store the positions of the cds and the introns {chromo:[[intron, start, end, strand], [cds, start, end, strand]]}

    intron_file = open(species_cds_intron, 'r')
    annotation = {}
    for line in intron_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[0] in annotation:
                annotation[line[0]].append([line[2], int(line[3]), int(line[4]), line[6]])
            else:
                annotation[line[0]] = [[line[2], int(line[3]), int(line[4]), line[6]]]
    intron_file.close()
    
    # make a dictionnary of mirna coordinates mir_coord = {'mir': [chromosome, start, end, strand]
    from mirna_data_for_meme import mirna_indices_from_file
    mirna_coord = mirna_indices_from_file(mirna_coordinates)

    # for each mirna, check if the mirna is located within an intron or not
    intronic_mirna = {}

    for mirna in mirna_coord:
        for chromo in annotation:
            if mirna_coord[mirna][0] == chromo and mirna_coord[mirna][1].isdigit():
                for i in range(len(annotation[chromo])):
                    if annotation[chromo][i][0] == 'intron':
                        if int(mirna_coord[mirna][1]) in range(annotation[chromo][i][1], annotation[chromo][i][2]) and mirna_coord[mirna][-1] == annotation[chromo][i][-1]:
                            intronic_mirna[mirna] = ['intronic', 'same strand']
                        elif int(mirna_coord[mirna][1]) in range(annotation[chromo][i][1], annotation[chromo][i][2]) and mirna_coord[mirna][-1] != annotation[chromo][i][-1]:
                            intronic_mirna[mirna] = ['intronic', 'opposite strand']

    for mirna in mirna_coord:
        if mirna not in intronic_mirna and mirna_coord[mirna][1].isdigit():
            intronic_mirna[mirna] = ['intergenic']
                
    return intronic_mirna


def count_intragenic_mirnas(species_cds_intron, mirna_coordinates, mirna_upstream, mast_output, intron = 'yes'):
    '''
    (file, file, file, file) -> None
    Print the number of intragenic and intergenic mirnas using mirna and gene coordinates files mirna_coordinates and parsed_gff_file
    for the mirnas using in the MEME analysis stored in the mirna_upstream file, separately for mirnas that have a motif in the mast_output and those that do not
    '''

    N_intragenic_same_strand_nomotif = 0
    N_intragenic_diff_strand_nomotif = 0
    N_intragenic_nomotif = 0
    N_intergenic_nomotif = 0

    N_intragenic_same_strand_motif = 0
    N_intragenic_diff_strand_motif = 0
    N_intragenic_motif = 0
    N_intergenic_motif = 0

    # make a dictionnary of mirnas that are used in the MEME analysis
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(mirna_upstream)

    # make a dictionnary with the mirna that have at least one motif
    from get_positions_from_mast import motif_position_from_mast
    mast_positions = motif_position_from_mast(mast_output)
    mirna_with_motif = mast_positions[1]

    # make a dictionnary storing the intragenic status of the mirna
    if intron == 'yes':
        intragenic_mirna = is_mirna_intronic(species_cds_intron, mirna_coordinates)
    elif intron == 'no':
        intragenic_mirna = is_intragenic(species_cds_intron, mirna_coordinates)

        
    
    for mirna in upstreams:
        if mirna in intragenic_mirna and mirna in mirna_with_motif:
            if intragenic_mirna[mirna][0] == 'intronic':
                N_intragenic_motif += 1
            if intragenic_mirna[mirna][0] == 'intronic' and intragenic_mirna[mirna][1] == 'same strand':
                N_intragenic_same_strand_motif += 1
            elif intragenic_mirna[mirna][0] == 'intronic' and intragenic_mirna[mirna][1] == 'opposite strand':
                N_intragenic_diff_strand_motif += 1
            else:
                N_intergenic_motif += 1
        elif mirna in intragenic_mirna and not mirna in mirna_with_motif:
            if intragenic_mirna[mirna][0] == 'intronic':
                N_intragenic_nomotif += 1
            if intragenic_mirna[mirna][0] == 'intronic' and intragenic_mirna[mirna][1] == 'same strand':
                N_intragenic_same_strand_nomotif += 1
            elif intragenic_mirna[mirna][0] == 'intronic' and intragenic_mirna[mirna][1] == 'opposite strand':
                N_intragenic_diff_strand_nomotif += 1
            else:
                N_intergenic_nomotif += 1

                
    print('N_intragenic_motif:\t{0}\tN_intergenic_motif:\t{1}'.format(N_intragenic_motif, N_intergenic_motif))
    print('N_intragenic_nomotif:\t{0}\tN_intergenic_nomotif:\t{1}'.format(N_intragenic_nomotif, N_intergenic_nomotif))
    
    print('N_intra_same_strand_motif:\t{0}\tN_intra_diff_strand_motif:\t{1}'.format(N_intragenic_same_strand_motif, N_intragenic_diff_strand_motif))
    print('N_intra_same_strand_nomotif:\t{0}\tN_intra_diff_strand_nomotif:\t{1}'.format(N_intragenic_same_strand_nomotif, N_intragenic_diff_strand_nomotif))

    print('N_upstream:\t{0}\tN_intra+inter:\t{1}'.format(len(upstreams), N_intragenic_motif + N_intergenic_motif + N_intragenic_nomotif + N_intergenic_nomotif))
    print(len(upstreams) == N_intragenic_motif + N_intergenic_motif + N_intragenic_nomotif + N_intergenic_nomotif)





# this includes code to find if a gene is within another gene
# use this code when the genome annotations doesn't have UTR annotations
# ie: gene coordinates = coordinates of the coding sequence + introns
def is_intragenic(species_cds_intron, mirna_coordinates):
    '''
    (file, file) -> dict
    Check if a mirna is located within a gene using the mirna and gene coordinates files parsed_gff_file and mirna_coordinates
    and return a dictionnary with this information for each mirna
    '''

    # make a dictionnary of gene coordinates gene_coord = {gene: [chromosome, start, end, strand]}
    gene_coord = get_gene_coordinates_from_gff(species_cds_intron)

    # make a dictionnary of mirna coordinates mir_coord = {'mir': [chromosome, start, end, strand]
    from mirna_data_for_meme import mirna_indices_from_file
    mirna_coord = mirna_indices_from_file(mirna_coordinates)

    # for each mirna, check if the mirna is located within a gene or not
    intragenic_mirna = {}

    for mirna in mirna_coord:
        for gene in gene_coord:
            if mirna_coord[mirna][0] == gene_coord[gene][0] and mirna_coord[mirna][1].isdigit():
                if int(mirna_coord[mirna][1]) in range(int(gene_coord[gene][1]), int(gene_coord[gene][2])) and mirna_coord[mirna][-1] == gene_coord[gene][-1]:
                    intragenic_mirna[mirna] = ['intronic', 'same strand']
                elif int(mirna_coord[mirna][1]) in range(int(gene_coord[gene][1]), int(gene_coord[gene][2])) and mirna_coord[mirna][-1] != gene_coord[gene][-1]:
                    intragenic_mirna[mirna] = ['intronic', 'opposite strand']


    for mirna in mirna_coord:
        if mirna not in intragenic_mirna and mirna_coord[mirna][1].isdigit():
            intragenic_mirna[mirna] = ['intergenic']
                
    return intragenic_mirna


# this function retrieves upstream sequence but doesn't filter upstream sequence < Length
##def grab_upstream_coding(random_genes, genome_file, Length):
##    '''
##    (dict, file, int) -> dict
##    Return a dictionnary of upstream sequence of given Length for each gene in the dictionnary random_genes
##    '''
##
##    import reverse_complement
##
##    # convert the genome_file into a genome dictionnary
##    # genome = {chromo: fasta}
##    from mirna_data_for_meme import convert_fasta
##    genome = convert_fasta(genome_file)
##
##    # use the coordinates stored in random_genes to fetch the upstream sequences
##    # random_genes = {i: [chromosome, start, end, strand, gene]}
##    upstream_gene_seqs = {}
##
##    for i in random_genes:
##        if random_genes[i][-2] == '+':
##            chromo = random_genes[i][0]
##            start = int(random_genes[i][1])
##            if start - Length >= 0:
##                upstream = genome[chromo][(start - Length):start]
##                upstream_gene_seqs[random_genes[i][-1]] = upstream
##            else:
##                upstream = genome[chromo][:start]
##                upstream_gene_seqs[random_genes[i][-1]] = upstream
##        elif random_genes[i][-2] == '-':
##            chromo = random_genes[i][0]
##            end = int(random_genes[i][2])
##            if end + Length <= len(genome[chromo]):
##                upstream = genome[chromo][end: (end + Length)]
##                upstream = reverse_complement.rev_compl(upstream)
##                upstream_gene_seqs[random_genes[i][-1]] = upstream
##            else:
##                upstream = genome[chromo][end:]
##                upstream = reverse_complement.rev_compl(upstream)
##                upstream_gene_seqs[random_genes[i][-1]] = upstream
##                
##    return upstream_gene_seqs

                                                                                                                                                 
