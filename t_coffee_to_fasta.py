def convert_tcoffee_to_fasta(filename):
    '''
    (str) -> None
    Read the t-coffee alignment files and save the alignment in fasta format in a text file
    '''

    tcoffee = open(filename, 'r')
    tcoffee.readline()
    tcoffee.readline()

    ali = {}

    for line in tcoffee:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if '*' in line[0]:
                continue
            elif line[0] in ali:
                ali[line[0]] += line[1]
            else:
                ali[line[0]] = line[1]

    alignment = open(filename[:-4] + '.txt', 'w')
    for gene in ali:
        alignment.write('>' + gene + '\n')
        alignment.write(ali[gene] + '\n')

    tcoffee.close()
    alignment.close()



def get_cbrcsp9_fasta(L):
    '''
    (list) -> None
    Convert each tcoffee files in the list L into a fasta file
    '''

    # hint: assign os.listdir to a variable to get the list of files

    for file in L:
        convert_tcoffee_to_fasta(file)

        



            
        

    


    
