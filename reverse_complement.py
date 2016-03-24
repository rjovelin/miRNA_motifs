def rev_compl(dna):
    '''
    (str) -> (str)
    Return the reverse complementary sequence of string dna

    >>> rev_compl('atcg')
    'cgat'
    '''

    dna2 = dna.upper()
    dna_comp = ''
    for i in dna2:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'
        elif i == 'N':
            dna_comp += 'N'

    reverse_comp_dna = ''
    for i in reversed(dna_comp):
        reverse_comp_dna += i

    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
        
    return reverse_comp_dna
