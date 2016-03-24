def compute_K_JC(seq1, seq2):
    '''
    (str, str) -> float
    Return the Jukes-Cantor distance between sequences seq1 and seq
    Prerequisite: seq1 and seq2 have same length

    >>> compute_K_JC('CAGAAGGCCTCGCCGAATATGACACTGTTGCGTAGCATGACCAATGAGACTCTACAAA-AGAGTGAATGGATCTGACTACGCTCAGTGGAATACCCGGCGAGTGCCTTCT',
    'CAGAAGGCCTCGCCGAATATGACACTGTTGCGTAGCATGACCAATGAGACTCTCCAAAGAGAGTGAATGGATCTGACTACGCTCAGTGGAATACCCGGCGAGTGCCTTCT')
    0.009231

    >>> compute_K_JC('AGCCGACGGAACGAGTAAATCTCATCCTAATCTGGTTG-ACACAACA-CAAGGAAGTGCTACCGATTTGGCTTGGGATTGACTTGTGGAAATGGCT',
    'AGCCGACGGAACGAGTAAATCTCATCCTAATCTGGTAGCACACAACAACAAGGAAGTGCTACCGATTTGGCTTGGGATTGACTTGTGGAATTGGCG')
    0.032614

    >>> compute_K_JC('TCAACGATGTGCTGCCTAACAGCGGATTCCAGTAACTCTGAGCATATACATGAAGTTTCAAAACAATGTAAATGCCAGTCGTTGCTGGAATTCGATAAGCACTCGTGAC',
    'TCAACGCTGTGCTGCCTAACAGCGGATTCCAGTAACTCTGAGCATATACATGAAGTTTCAAAACAATGTAAATGCCAGTCGTTGCTGGAATTCGATAAGCACTCGTGAC')
    0.009231
    '''

    # test prerequisite of equal length
    assert len(seq1) == len(seq2), 'The sequences have different lengths, check alignment and/or missing end gaps'

    import math


    # accept - or ~ as gaps
    SEQ1 = seq1.replace('~', '-')
    SEQ2 = seq2.replace('~', '-')
    # make sure the sequences have the same case 
    SEQ1 = seq1.upper()
    SEQ2 = seq2.upper()

    D = 0
    gap = 0
    N = 0

    # count differences between seq1 and seq2, gapped positions are excluded, positions with N are excluded
    for i in range(len(SEQ1)):
        if SEQ1[i] == '-' or SEQ2[i] == '-':
            gap += 1
        elif SEQ1[i] == 'N' or SEQ2[i] == 'N':
            N += 1
        elif SEQ1[i] != SEQ2[i]:
            D += 1
    L = len(SEQ1) - gap - N
    # set a threshold of 30% gaps and Ns allowed within a given window
    # if more then don't compute divergence
    # K is not defined if p >= 0.75
    if L != 0 and (gap + N) <= (0.3 * L):
        p = D/ L
        if p < 0.75:
            K = (-3/4) * math.log(1-(4/3 * p))
            return round(abs(K), 6)
            

    
# use this function to compute K
def compute_K(seq1, seq2):
    '''
    (str, str) -> float
    Return the p distance between sequences seq1 and seq
    Prerequisite: seq1 and seq2 have same length
    '''

    # test prerequisite of equal length
    assert len(seq1) == len(seq2), 'The sequences have different lengths, check alignment and/or missing end gaps'

    # accept - or ~ as gaps
    SEQ1 = seq1.replace('~', '-')
    SEQ2 = seq2.replace('~', '-')
    # make sure the sequences have the same case 
    SEQ1 = seq1.upper()
    SEQ2 = seq2.upper()

    D = 0
    gap = 0
    N = 0

    # count differences between seq1 and seq2, gapped positions are excluded, positions with N are excluded
    for i in range(len(SEQ1)):
        if SEQ1[i] == '-' or SEQ2[i] == '-':
            gap += 1
        elif SEQ1[i] == 'N' or SEQ2[i] == 'N':
            N += 1
        elif SEQ1[i] != SEQ2[i]:
            D += 1
    L = len(SEQ1) - gap - N
    # set a threshold of 30% gaps and Ns allowed within a given window
    # if more then don't compute divergence
    # K is not defined if p >= 0.75
    if L != 0 and (gap + N) <= (0.3 * L):
        K = D/ L
        return K
    else:
        return None




    
