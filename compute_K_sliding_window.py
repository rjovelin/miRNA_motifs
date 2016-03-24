
import numpy as np
from scipy import stats
import math



def compute_K_sliding_window(gene1, gene2, window, step):
    '''
    (str, str, int, int) -> dict
    Compute pairwise divergence between 2 sequences with a sliding window
    of size window and incrementing by step
    '''

    # accept ~ or - as gaps
    # make sure the 2 sequences have same case
    Gene1 = gene1.replace('~', '-').upper()
    Gene2 = gene2.replace('~', '-').upper()

    windows_K = {}
    
    # go through the sequence, stop at length(seq) - window
    # get the windows in turn
    for i in range(0, len(Gene1)-window, step):
        w1 = Gene1[i: i + window]
        w2 = Gene2[i: i + window]
        D = 0
        gap = 0
        N = 0
        # compute divergence K between windows at same positions
        # round K to 6 decimals
        for j in range(len(w1)):
            if w1[j] == '-' or w2[j] == '-':
                gap += 1
            elif w1[j] == 'N' or w2[j] == 'N':
                N += 1
            elif w1[j] != w2[j]:
                D += 1
        L = len(w1) - gap - N
        # set a threshold of 30% gaps allowed within a given window
        # if more then don't compute divergence
        if L != 0 and (gap + N) <= (0.3 * L):
            K = D/L
            K = round(abs(K), 6)
            # store K in the dictionnary
            if i in windows_K:
                windows_K[i].append(K)
            else:
                windows_K[i] = [K]

    return windows_K


def reverse_sequence(dna):
    '''
    (str) -> str
    Reverse the order of sequence dna and return the reverse sequence

    >>> reverse_sequence('atgcata')
    atacgta
    '''

    rev_dna = ''
    for base in reversed(dna):
        rev_dna += base

    return rev_dna
        
        
# this function computes sliding windows specifically for upstream sequences
# by aligning windows to the right (reading from right to left)
# the index are reversed: ie: index 0 means end of the original sequence and
# region just upstream the focus sequence
def sliding_windows_K_all_genes_upstream(L_genes, window, step):
    '''
    (list) -> dict
    Take a list of filenames for files each containing 2 aligned sequences,
    compute sliding window divergence for each file and return a dictionnary
    storing all divergence values for each position
    '''

    # hint: assign os.listdir to a variable to get the list of filenames to use

    # make a dictionnary with identical keys to each windows_K dictionnary
    # create many more keys since each windows_K has different length
    # here length of each windows_K is < 2500

    all_K = {}
    for i in range(0, (2500 - window), step):
        all_K[i] = []



    # open each file
    # convert the 2 fasta seqs into 2 single string seqs
    for filename in L_genes:
        myfile = open(filename, 'r')
        genes = {}
        for line in myfile:
            line = line.rstrip()
            if line == '':
                continue
            elif line.startswith('>'):
                genes[line[1:]] = ""
                seq_name = line[1:]
            else:
                genes[seq_name] += line
        # assign the sequences in dictionnary genes to variables gene1 and gene2        
        gene_names = []
        for gene in genes:
            gene_names.append(gene)
        gene1 = genes[gene_names[0]]
        gene2 = genes[gene_names[1]]

        # need to take the reverse of the sequence
        # so that all sequences are aligned to the right
        gene1 = reverse_sequence(gene1)
        gene2 = reverse_sequence(gene2)

        myfile.close()

        # compute the silding window divergence
        sliding_K = compute_K_sliding_window(gene1, gene2, window, step)

        
        #update the all_K dictionnary with the sliding_K dictionnary
        for position in sliding_K:
            all_K[position].extend(sliding_K[position])

    # get rid of the keys with empty list we added at the beginning
    get_rid = []
    for position in all_K:
        if all_K[position] == []:
            get_rid.append(position)
    for item in get_rid:
        del all_K[item]

    # delete windows with small sample size: N < 10
    small = []
    for position in all_K:
        if len(all_K[position]) < 50:
            small.append(position)
    for item in small:
        del all_K[item]


    return all_K
   

def save_K_sliding_statistics_to_file(all_K, stats_output):
    '''
    Take a dictionnary all_K containing divergence for many genes in each window,
    compute the mean and its 95% confidence interval for every position using t_values
    and save it into a file stats_output.
    '''
    
    # create a lambda function to convert values to strings
    Gstr = lambda x: str(x)
      

    # make a dictionnary positions : list with mean K and 95%CI
    diverg = {}
    
    # loop over positions in dict    
    for i in all_K:
        # compute the mean K at that position
        mean_K = np.mean(all_K[i])
        # compute the standard error
        stderror = np.std(all_K[i]) / math.sqrt(len(all_K[i]))
        # compute the margin error (critical value = 1.96 for a 95% CI)
        margin = 1.96 * stderror
        lCI = mean_K - margin
        hCI = mean_K + margin
        diverg[i] = [i, len(all_K[i]), mean_K, stderror, np.std(all_K[i]), lCI, hCI]
        
    # write results to file
    outputfile = open(stats_output, 'w')
    outputfile.write('\t'.join(['position', 'N', 'mean', 'SEM', 'stdev', 'low_95%_CI', 'high_95%_CI']) + '\n')
    # make a list of positions
    positions = [i for i in diverg]
    # sort list
    positions.sort()
    # reverse list (index 0 means end of upstream sequence, position adjacent to premirna)
    positions.reverse()    
    
    # loop over indices in positions 
    for i in positions:
        outputfile.write('\t'.join(list(map(Gstr, diverg[i]))) + '\n')
        
    outputfile.close()
