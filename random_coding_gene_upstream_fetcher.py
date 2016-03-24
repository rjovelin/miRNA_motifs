from get_upstream_coding_genes import get_gene_coordinates_from_gff
from get_upstream_coding_genes import upstream_coding_to_file



def grab_upstream_random_coding_genes(gene_coord, genome_file, gene_number, Length):
    '''
    (dict, file, int, int) -> dict

    Return a dictionnary of upstream sequence of given Length for each gene taken at random
    with upstream length  = 1000
    '''

    import reverse_complement


    # convert the genome_file into a genome dictionnary
    # genome = {chromo: fasta}
    from mirna_data_for_meme import convert_fasta
    genome = convert_fasta(genome_file)


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
